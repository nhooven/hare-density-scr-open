# PROJECT: Closed multi-session SCR
# SCRIPT: 05b - CBB-only model (sex ratio)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Feb 2026
# COMPLETED: 
# LAST MODIFIED: 13 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(nimble)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("for_model/CBB_test/constants.rds")
data.list <- readRDS("for_model/CBB_test/data.rds")
state.inits <- readRDS("for_model/CBB_test/state_inits.rds")

# ______________________________________________________________________________
# 3. Write model ----
# ______________________________________________________________________________

model.code <- nimbleCode({
  
  # OPEN SUB-MODEL
  # priors
  # phi - persistence
  phi ~ dunif(0, 1)
  
  # rho - per capita recruitment
  rho ~ dunif(0, 7)
  
  # gamma - probability of available individual being recruited [u, t]
  # loop over sites U
  for (u in 1:U) {
    
    # loop over years
    for (t in 2:YR) {
      
      # logit constraint is important
      logit(gamma[u, t - 1]) <- log(N[u, t - 1] + 1e-6) + 
                                log(rho) - 
                                log(n.avail[u, t - 1] + 1e-6)
      
    } # YR
    
  } # U
  
  # omega1 - initial states [1:3, U]
  for (u in 1:U) {
    
    # psi - initial inclusion [U]
    psi[u] ~ dunif(0, 1)
    
    omega1[1, u] <- 1 - psi[u]       # available with p = 1 - psi
    omega1[2, u] <- psi[u]           # recruited with p = psi
    omega1[3, u] <- 0.0              # cannot start dead
    
  } # U
  
  # omega - transition matrix [3 x 3, U, YR - 1]
  # rows: state at t
  # columns: state at t + 1
  for (u in 1:U) {
    
    for (t in 1:(YR - 1)) {
      
      # ENTERED
      omega[2, 1, u, t] <- 0.0
      omega[2, 2, u, t] <- phi
      omega[2, 3, u, t] <- 1 - phi
      
      # DIED
      omega[3, 1, u, t] <- 0.0
      omega[3, 2, u, t] <- 0.0
      omega[3, 3, u, t] <- 1.0
      
      # NOT ENTERED
      omega[1, 1, u, t] <- 1 - gamma[u, t]
      omega[1, 2, u, t] <- gamma[u, t]
      omega[1, 3, u, t] <- 0.0
      
    } # YR - 1
    
  } # U
  
  # open likelihood
  # loop over individuals
  for (i in 1:M) {
    
    # initial state
    z[i, 1] ~ dcat(omega1[1:3, site[i]])
    
    # subsequent states
    # loop over primary occasions [YR]
    for (t in 1:(YR - 1)) {
      
      z[i, t + 1] ~ dcat(omega[z[i, t], 1:3, site[i], t])
      
    } # YR
    
  } # M
  
  # CLOSED SCR SUB-MODEL
  # detection
  # alpha0 - baseline detection (log scale)
  alpha0_b0 ~ dnorm(0, sd = 1)  # intercept
  alpha0_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # alpha2 - trap response
  alpha2_b0 ~ dnorm(0, sd = 1)  # intercept
  alpha2_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # sigma - spatial scale of movement
  sigma_b0 ~ dnorm(log(45), sd = 0.5)     # intercept
  sigma_b1 ~ dnorm(0, sd = 1)             # effect of male
  
  # pooled sex probability
  psi.sex ~ dunif(0, 1)
  
  # calculate probabilities for closed likelihood
  # loop over individuals [M]
  for (i in 1:M) {
    
    # sex indicator (every latent individual gets one value)
    sex[i] ~ dbern(psi.sex)
    
    # detection
    alpha2[i] <- alpha2_b0 + alpha2_b1 * sex[i]
    log(sigma[i]) <- sigma_b0 + sigma_b1 * sex[i]
    alpha1[i] <- -1 / sigma[i]
    
    # loop over primary occasions [YR]
    for (t in 1:YR) {
      
      # z2 - indicator for state 2
      z2[i, t] <- step(z[i, t] - 1.5) * step(2.5 - z[i, t])
      
      # s - activity centers [M, 2, YR]
      s[i, 1, t] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
      s[i, 2, t] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
      
      # vectorized trap calculations
      # d - distances between s and each trap j
      d[i, 1:J, t] <- sqrt(pow(s[i, 1, t] - trap.coords[1:J, 1, site[i]], 2) + 
                           pow(s[i, 2, t] - trap.coords[1:J, 2, site[i]], 2))
      
      # g - distance kernel
      g[i, 1:J, t] <- exp(alpha1[i] * d[i, 1:J, t])
      
      # loop over secondary occasions [K] (indexed by site and year)
      for (k in 1:K[site[i], t]) {
        
        # alpha0 - baseline hazard of detection [i, k, t]
        log(alpha0[i, k, t]) <- alpha0_b0 + alpha0_b1 * sex[i] + 
                                alpha2[i] * prev.cap[i, k, t]
        
        # and loop over traps [J]
        for (j in 1:J) {
          
          # eta - un-normalized probability 
          eta[i, k, j, t] <- alpha0[i, k, t] * g[i, j, t] *
            
            # inclusion (1 if state == 2, 0 otherwise)
            z2[i, t] *
            
            # trap operation
            trap.op[j, k, t, site[i]] *
            
            # trap deaths
            trap.deaths[i, k, t]
          
          # p - normalized probabilities for categorical likelihood
          p[i, k, j, t] <- eta[i, k, j, t] / (1 + sum(eta[i, k, 1:J, t]))  # sum over all traps
          
        } # J
        
        # probability of not being captured as the complement of all trap-specific probs
        p[i, k, J + 1, t] <- 1 - sum(p[i, k, 1:J, t])
        
      } # K
      
    } # YR
    
  } # M
  
  # closed likelihood for observed individuals
  for (i in 1:n) {
    
    # loop over years
    for (t in 1:YR) {
      
      # loop over occasions
      for (k in 1:K[site[i], t]) {
        
        ch[i, k, t] ~ dcat(p[i, k, 1:(J + 1), t])
        
      } # K
      
    } # YR
    
  } # n
  
  # closed likelihood for non-detected individuals
  for (i in (n + 1):M) {
    
    # loop over years
    for (t in 1:YR) {
      
      zeroes[i, t] ~ dbern(1 - prod(1 - p[i, 1:K[site[i], t], 1:J, t]))
      
    } # YR
    
  } # (n + 1):M
  
  # DERIVED PARAMETERS
  # N
  # loop over sites
  for (u in 1:U) {
    
    # loop over years
    for (t in 1:YR) {
      
      # sex-specific N
      # total recruited individuals
      Nf[u, t] <- sum((z[1:M, t] == 2) * (site[1:M] == u) * (sex[1:M] == 0))
      Nm[u, t] <- sum((z[1:M, t] == 2) * (site[1:M] == u) * (sex[1:M] == 1))
      
      # total available individuals
      n.availf[u, t] <- sum((z[1:M, t] == 1) * (site[1:M] == u) * (sex[1:M] == 0))
      n.availm[u, t] <- sum((z[1:M, t] == 1) * (site[1:M] == u) * (sex[1:M] == 1))
      
      # pooled
      # total recruited individuals
      N[u, t] <- Nf[u, t] + Nm[u, t]
      
      # total available individuals
      n.avail[u, t] <- n.availf[u, t] + n.availm[u, t]
      
      # R - sex ratio (proportion male)
      R[u, t] <- Nm[u, t] / N[u, t]
      
    } # YR
    
  } # U
  
})

# ______________________________________________________________________________
# 4. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # OPEN
  phi = runif(1, 0, 1),
  rho = runif(1, 0, 7),
  psi = runif(3, 0, 1),
  z = state.inits[[1]],       # change based on which chain we're running
  
  # CLOSED
  s = array(c(cbind(runif(constant.list$M, -200, 200),
                    runif(constant.list$M, -200, 200)),
              cbind(runif(constant.list$M, -200, 200),
                    runif(constant.list$M, -200, 200)),
              cbind(runif(constant.list$M, -200, 200),
                    runif(constant.list$M, -200, 200))),
            dim = c(constant.list$M, 2, constant.list$YR)),
  alpha0_b0 = rnorm(1, 0, 1),
  alpha0_b1 = rnorm(1, 0, 1),
  alpha2_b0 = rnorm(1, 0, 1),
  alpha2_b1 = rnorm(1, 0, 1),
  sigma_b0 = rnorm(1, log(45), 0.5),
  sigma_b1 = rnorm(1, 0, 1),
  sex = ifelse(is.na(data.list$sex) == F, NA, 0)
  
)

# ______________________________________________________________________________
# 5. Parameters to monitor ----
# ______________________________________________________________________________

monitor <- c(
  
  "psi", "phi", "rho", "gamma",
  "alpha0_b0", "alpha0_b1", "alpha2_b0", "alpha2_b1", "sigma_b0", "sigma_b1",
  "N", "n.avail", "R"
  
)

# ______________________________________________________________________________
# 6. Set up and run model----
# ______________________________________________________________________________

# set up model (this takes forever!)
model.1 <- nimbleModel(
  
  code = model.code,
  constants = constant.list,
  data = data.list,
  inits = inits,
  calculate = F
  
)

model.1.conf <- configureMCMC(model.1, monitors = monitor)

# block sampling
# detection
model.1.conf$removeSamplers(c("alpha0_b0", "alpha2_b0", "sigma_b0"))
model.1.conf$addSampler(c("alpha0_b0", "alpha2_b0", "sigma_b0"), 
                        type = "RW_block",
                        control = list("propCov" = matrix(c(0.5, -0.3, -0.05,
                                                            -0.3, 0.3, 0.01,
                                                            -0.05, 0.01, 0.01),
                                                          nrow = 3),
                                       adaptScaleOnly = F))

# activity centers
model.1.conf$removeSamplers(c("s"))

for(i in 1:constant.list$M) {
  
  for (t in 1:constant.list$YR) {
    
    snew = paste0("s[", i, ", 1:2, ", t, "]")
    
    model.1.conf$addSampler(snew, type = "RW_block", silent = T)
    
  }
  
}

# build model
model.1.mcmc <- buildMCMC(conf = model.1.conf, monitors = monitor)

# compile model
compileNimble(model.1)
model.1.comp <- compileNimble(model.1.mcmc, project = model.1)

# run MCMC
model.1.run <- runMCMC(
  
  mcmc = model.1.comp,
  niter = 30000,
  nburnin = 15000,
  nchains = 1,
  thin = 10,
  samplesAsCodaMCMC = TRUE
  
)