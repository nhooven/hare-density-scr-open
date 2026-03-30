# PROJECT: Open multi-session SCR
# SCRIPT: 05a - Full open model (simplified)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Feb 2026
# COMPLETED: 
# LAST MODIFIED: 25 Mar 2026
# R VERSION: 4.4.3

# Here, we'll drop the cluster random intercepts for the detection parameters
# I just don't think these will be well-identified

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(nimble)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("for_model/constants.rds")
data.list <- readRDS("for_model/data.rds")
state.inits <- readRDS("for_model/state_inits.rds")

# which state inits? 
chain = 1

state.inits <- state.inits[[chain]]

# ______________________________________________________________________________
# 3. Write code ----
# ______________________________________________________________________________
# 3a. Functions ----
# ______________________________________________________________________________

# cat_p - calculate categorical probability vector
cat_p <- nimbleFunction(
  
  run = function (
    
    # detection parameters
    # stochastic from priors, all scalar
    # baseline hazard of detection
    lam0_b0 = double(0),
    lam0_b = double(1),
    
    # previous capture
    alpha2_b0 = double(0),
    alpha2_b = double(1),
    
    # exponential decay
    alpha1_b0 = double(0),
    alpha1_b = double(1),
    
    # individual data
    # scalars
    sex = double(0),          # sex[i]
    ret = double(0),          # det.ret[i, t]
    pil = double(0),          # det.pil[i, t]
    z = double(0),            # z[i, t]
    prev.cap = double(0),     # prev.cap[i, k, t]
    trap.deaths = double(0),  # trap.deaths[i, k ,t]
    
    # vectors
    s = double(1),            # [2] - s[i, 1:2, t]
    trap.op = double(1),      # [J] trap.op[1:J, k, t, site[i]]
    
    # scalar constants
    J = integer(0),            # J
    
    # trap coordinates - matrix [36, 2, 1]
    trap.coords = double(2) # trap.coords[1:J, 1:2, site[i]]
    
  ) {
    
    # declare data types
    p <- nimNumeric(length = J + 1, value = 0.0)
    eta <- nimNumeric(length = J, value = 0.0)
    eta.denom <- 1.0
    p.sum <- 0.0       # for subtracting, "trap 37"
    
    # detection
    # alpha2 - previous capture effect
    alpha2 <- alpha2_b0 + 
      
      alpha2_b[1] * sex + 
      alpha2_b[2] * ret + 
      alpha2_b[3] * pil
    
    # alpha1 - distance decay
    # we want this to be between zero and 1
    # we'll make it negative later
    logit_alpha1 <- alpha1_b0 + 
      
      alpha1_b[1] * sex + 
      alpha1_b[2] * ret + 
      alpha1_b[3] * pil
    
    alpha1 <- 1 / (1 + exp(-logit_alpha1))
    
    # lam0 - baseline hazard of detection [i, k, t]
    lam0 <- exp(
      
      lam0_b0 + 
        
        lam0_b[1] * sex + 
        lam0_b[2] * ret + 
        lam0_b[3] * pil + 
        
        alpha2 * prev.cap
      
    )
    
    # loop over traps [J]
    for (j in 1:J) {
      
      # trap calculations
      # d - distances between s and each trap j
      dx <- s[1] - trap.coords[j, 1]
      dy <- s[2] - trap.coords[j, 2]
      d <- sqrt(dx * dx + dy * dy)
      
      # eta - un-normalized probability 
      eta[j] <- lam0 * exp(-alpha1 * d) *
        
        # inclusion (1 if state == 2, 0 otherwise)
        z *
        
        # trap operation
        trap.op[j] *
        
        # trap deaths
        trap.deaths
      
      # add eta in the denominator to sum over all traps
      eta.denom <- eta.denom + eta[j]
      
    } # J
    
    # normalize probabilities
    for (j in 1:J) {
      
      # p - normalized probabilities for categorical likelihood
      p[j] <- eta[j] / eta.denom
      
      # increment p to the p.sum
      p.sum <- p.sum + p[j]
      
    }
    
    # probability of not being captured as the complement of all trap-specific probs
    p[J + 1] <- 1.0 - p.sum
    
    # return type (vector) - [J + 1]
    returnType(double(1))
    
    return(p)
    
  }
  
)

# bern_p - calculate Bernoulli probability for unobserved individuals
bern_p <- nimbleFunction(
  
  run = function (
    
    # detection parameters
    # baseline hazard of detection
    lam0_b0 = double(0),
    lam0_b = double(1),
    
    # previous capture
    alpha2_b0 = double(0),
    alpha2_b = double(1),
    
    # exponential decay
    alpha1_b0 = double(0),
    alpha1_b = double(1),
    
    # individual data
    # scalars
    sex = double(0),          # sex[i]
    ret = double(0),          # det.ret[i, t]
    pil = double(0),          # det.pil[i, t]
    z = double(0),            # z2[i, t]
    
    # vectors
    s = double(1),            # [2] - s[i, 1:2, t]
    prev.cap = double(1),     # [K] prev.cap[i, 1:K, t]
    trap.deaths = double(1),  # [K] trap.deaths[i, 1:K,t]
    
    # matrices
    trap.op = double(2),      # trap.op[1:J, 1:K, t, site[i]]
    
    # scalar constants
    J = integer(0),          # J
    K = integer(0),          # K[site[i]]
    
    # trap coordinates - matrix [36, 2, 1]
    trap.coords = double(2) # trap.coords[1:J, 1:2, site[i]]
    
  ) {
    
    # declare data types
    lam0 <- nimNumeric(K, 0.0)
    eta <- nimMatrix(0.0, nrow = K, ncol = J)
    p <- nimMatrix(0.0, nrow = K, ncol = J + 1)
    eta.denom <- nimNumeric(length = K, 0.0)
    p.sum <- nimNumeric(K, 0.0)       # for subtracting, "trap 37"
    
    # detection
    # alpha2 - previous capture effect
    alpha2 <- alpha2_b0 +
      
      alpha2_b[1] * sex + 
      alpha2_b[2] * ret + 
      alpha2_b[3] * pil
    
    # alpha1 - distance decay
    logit_alpha1 <- alpha1_b0 + 
      
      alpha1_b[1] * sex + 
      alpha1_b[2] * ret + 
      alpha1_b[3] * pil
    
    alpha1 <- 1 / (1 + exp(-logit_alpha1))
    
    # loop over secondary occasions K
    for (k in 1:K) {
      
      eta.denom[k] <- 1.0
      
      # lam0 - baseline hazard of detection [i, k, t]
      lam0[k] <- exp(
        
        lam0_b0 + 
          
          lam0_b[1] * sex + 
          lam0_b[2] * ret + 
          lam0_b[3] * pil + 
          
          alpha2 * prev.cap[k]
        
      )
      
      # and loop over traps [J]
      for (j in 1:J) {
        
        # trap calculations
        # d - distances between s and each trap j
        dx <- s[1] - trap.coords[j, 1]
        dy <- s[2] - trap.coords[j, 2]
        d <- sqrt(dx * dx + dy * dy)
        
        # eta - un-normalized probability 
        eta[k, j] <- lam0[k] * exp(-alpha1 * d) *
          
          # inclusion (1 if state == 2, 0 otherwise)
          z *
          
          # trap operation
          trap.op[j, k] *
          
          # trap deaths
          trap.deaths[k]
        
        # add eta in the denominator to sum over all traps
        eta.denom[k] <- eta.denom[k] + eta[k, j]
        
      } # J
      
      # normalize
      for (j in 1:J) {
        
        # p - normalized probabilities for categorical likelihood
        p[k, j] <- eta[k, j] / eta.denom[k]
        
        # increment p to the p.sum
        p.sum[k] <- p.sum[k] + p[k, j]
        
      } # J
      
      # probability of not being captured as the complement of all trap-specific probs
      p[k, J + 1] <- 1.0 - p.sum[k]
      
    } # K
    
    # multiplying probs is the same as adding log-probs
    # declare starting value
    log.prob <- 0.0
    
    # loop through K
    for (k in 1:K) {
      
      # loop through J
      for (j in 1:J) {
        
        log.prob <- log.prob + log(1 - p[k, j])
        
      }
      
    }
    
    # exponentiate and subtract from 1
    prob <- 1 - exp(log.prob)
    
    # return
    returnType(double(0))
    return(prob)
    
  }
  
)

# ______________________________________________________________________________
# 3b. Model ----
# ______________________________________________________________________________

model.code <- nimbleCode({
  
  # ___________________________
  # OPEN STATE-SPACE SUB-MODEL
  # ___________________________
  
  # t1 -> t2 -> t3 -> t4
  #   tr1   tr2   tr3
  
  # demographic parameter priors
  # phi - persistence (logit scale)
  phi_b0 ~ dlogis(0, 1)
  phi_sd ~ T(dt(0, sigma = 1, df = 1), 0, )    # SD - half-Cauchy
  
  # rho - per capita recruitment (log scale)
  rho_b0 ~ dunif(log(0.001), log(7))
  rho_sd ~ T(dt(0, sigma = 1, df = 1), 0, )    # SD - half-Cauchy
  
  # linear coefficients on phi  (n = 7)
    # b1 - male effect
    # b2 - post1 effect
    # b3 - post2 effect
    # b4 - ret x post1 effect
    # b5 - ret x post2 effect
    # b6 - pil x post1 effect
    # b7 - pil x post2 effect
  for (x in 1:7) {
    
    phi_b[x] ~ dlogis(0, 1)        # logistic prior to avoid logit woes
    
  }
  
  # linear coefficients on rho parameters (n = 7)
    # b1 - post1 effect
    # b2 - post2 effect
    # b3 - ret x post1 effect
    # b4 - ret x post2 effect
    # b5 - pil x post1 effect
    # b6 - pil x post2 effect
  for (x in 1:6) {
    
    rho_b[x] ~ dnorm(0, sd = 1)
    
  }
  
  # cluster-specific random scaling factors
  for (c in 1:4) {
    
    phi_c[c] ~ dnorm(0, sd = 1)
    rho_c[c] ~ dnorm(0, sd = 1)
    
  }
  
  # calculate phi and rho by individual
  for (i in 1:M) {
    
    for (t in first.year[site[i]]:YR) {
      
      # phi
      phi[i, t] <- ilogit(
        
        (phi_b0 + phi_sd * phi_c[cluster[i]]) +
          
          phi_b[1] * sex[i] +
          phi_b[2] * op.post1[i, t] +
          phi_b[3] * op.post2[i, t] +
          phi_b[4] * op.ret[i] * op.post1[i, t] +
          phi_b[5] * op.ret[i] * op.post2[i, t] +
          phi_b[6] * op.pil[i] * op.post1[i, t] +
          phi_b[7] * op.pil[i] * op.post2[i, t]
        
      )
      
      # rho
      rho[i, t] <- exp(
        
        (rho_b0 + rho_sd * rho_c[cluster[i]]) +
          
          rho_b[1] * op.post1[i, t] +
          rho_b[2] * op.post2[i, t] +
          rho_b[3] * op.ret[i] * op.post1[i, t] +
          rho_b[4] * op.ret[i] * op.post2[i, t] +
          rho_b[5] * op.pil[i] * op.post1[i, t] +
          rho_b[6] * op.pil[i] * op.post2[i, t]
        
      )
      
    } # YR
    
  } # M
  
  # gamma - probability of entry
  # conditional on how many individuals were in states 2 and 1 in t - 1
  # calculate by individual-year
  for (i in 1:M) {
    
    for (t in (first.year[site[i]] + 1):YR) {
      
      logit(gamma[i, t - 1]) <- log(N[site[i], t - 1] + 1e-6) +
                                log(rho[i, t - 1]) -
                                log(N.avail[site[i], t - 1] + 1e-6)
      
    } # YR
    
  } # M
  
  # omega1 - initial state matrix
  for (u in 1:U) {
    
    # psi - initial inclusion [U]
    psi[u] ~ dunif(0, 1)
    
    # omega1 - initial states [1:3, U]
    # this will be for the first year BY SITE
    omega1[1, u] <- 1 - psi[u]          # available with p = 1 - psi
    omega1[2, u] <- psi[u]              # recruited with p = psi
    omega1[3, u] <- 0.0                 # cannot start dead
    
  } # U
  
  # open likelihood [M]
  for (i in 1:M) {
    
    # initial state (first.year by site)
    z[i, first.year[site[i]]] ~ dcat(omega1[1:3, site[i]])
    
    # subsequent states (YR > first.year)
    for (t in (first.year[site[i]]):(YR - 1)) {
      
      # omega - open state probabilities [3, 3, M, t]
      # rows: state at t
      # columns: state at t + 1
      # ENTERED
      omega[2, 1, i, t] <- 0.0
      omega[2, 2, i, t] <- phi[i, t]
      omega[2, 3, i, t] <- 1 - phi[i, t]
      
      # DIED
      omega[3, 1, i, t] <- 0.0
      omega[3, 2, i, t] <- 0.0
      omega[3, 3, i, t] <- 1.0
      
      # NOT ENTERED
      omega[1, 1, i, t] <- 1 - gamma[i, t]
      omega[1, 2, i, t] <- gamma[i, t]
      omega[1, 3, i, t] <- 0.0
      
      z[i, t + 1] ~ dcat(omega[z[i, t], 1:3, i, t])
      
    } # YR - 1
    
  } # M
  
  # ___________________________
  # CLOSED SCR SUB-MODEL
  # ___________________________
  
  # detection parameter priors
  # lam0 - baseline hazard detection (log scale)
  lam0_b0 ~ dnorm(0, sd = 1)  # mean 
  
  # alpha2 - trap response
  alpha2_b0 ~ dnorm(0, sd = 1)  # mean 
  
  # alpha1 - distance decay 
  alpha1_b0 ~ dlogis(0, 1)     # mean
  
  # detection linear coefficients (n = 3)
  # b1 - male effect
  # b2 - ret effect
  # b3 - pil effect
  for (x in 1:3) {
    
    lam0_b[x] ~ dnorm(0, sd = 1)  
    alpha2_b[x] ~ dnorm(0, sd = 1)  
    alpha1_b[x] ~ dlogis(0, 1)
    
  }
  
  # pooled sex probability
  psi.sex ~ dunif(0, 1)
  
  # calculations and likelihood for observed individuals [n]
  for (i in 1:n) {
    
    # sex indicator
    sex[i] ~ dbern(psi.sex)
    
    for (t in (first.year[site[i]]):YR) {
      
      # z indicator
      z2[i, t] <- step(z[i, t] - 1.5) * step(2.5 - z[i, t])
      
      # s - activity centers [M, 2, YR]
      s[i, 1, t] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
      s[i, 2, t] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
      
      for (k in 1:K[site[i], t]) {
        
        # calculate categorical probability vector 
        p[i, 1:(J + 1), k, t] <- cat_p(
          
          # baseline hazard of detection
          lam0_b0 = lam0_b0,
          lam0_b = lam0_b[1:3],
          
          # previous capture effect
          alpha2_b0 = alpha2_b0,
          alpha2_b = alpha2_b[1:3],
          
          # spatial scale of movement
          alpha1_b0 = alpha1_b0,
          alpha1_b = alpha1_b[1:3],
          
          # constants
          sex = sex[i],
          ret = det.ret[i, t],
          pil = det.pil[i, t],
          z = z2[i, t],
          s = s[i, 1:2, t],
          prev.cap = prev.cap[i, k, t],
          trap.deaths = trap.deaths[i, k, t],
          trap.op = trap.op[ , k, t, site[i]],
          J = J,
          trap.coords = trap.coords[ , 1:2, site[i]]
          
        )
        
        # categorical likelihood
        ch[i, k, t] ~ dcat(p[i, 1:(J + 1), k, t])
        
      } # K
      
    } # YR
    
  } # n
  
  # calculations and likelihood for unobserved individuals ("zeroes" trick) [(n + 1):M]
  for (i in (n + 1):M) {
    
    # sex indicator
    sex[i] ~ dbern(psi.sex)
    
    for (t in (first.year[site[i]]):YR) {
      
      # z indicator
      z2[i, t] <- step(z[i, t] - 1.5) * step(2.5 - z[i, t])
      
      # s - activity centers [M, 2, YR]
      s[i, 1, t] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
      s[i, 2, t] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
      
      # calculate Bernoulli probability
      p.zero[i, t] <- bern_p(
        
        # baseline hazard of detection
        lam0_b0 = lam0_b0,
        lam0_b = lam0_b[1:3],
        
        # previous capture effect
        alpha2_b0 = alpha2_b0,
        alpha2_b = alpha2_b[1:3],
        
        # spatial scale of movement
        alpha1_b0 = alpha1_b0,
        alpha1_b = alpha1_b[1:3],
        
        # constants
        sex = sex[i],
        ret = det.ret[i, t],
        pil = det.pil[i, t],
        z = z2[i, t],
        s = s[i, 1:2, t],
        prev.cap = prev.cap[i, 1:8, t],
        trap.deaths = trap.deaths[i, 1:8, t],
        trap.op = trap.op[ , 1:8, t, site[i]],
        J = J,
        K = K[site[i], t],
        trap.coords = trap.coords[ , 1:2, site[i]]
        
      )
      
      # Bernoulli likelihood
      zeroes[i, t] ~ dbern(p.zero[i, t])
      
    } # YR
    
  } # M
  
  # ___________________________
  # DERIVED QUANTITIES
  # ___________________________
  
  # N - counts of state 2 individuals [U, YR]
  # N.avail - counts of state 1 individuals [U, YR]
  for (u in 1:U) {
    
    for (t in 1:YR) {
      
      N[u, t] <- inprod(step(z[1:M, t] - 1.5) * step(2.5 - z[1:M, t]), which.site[1:M, u])
      N.avail[u, t] <- inprod(step(z[1:M, t] - 0.5) * step(1.5 - z[1:M, t]), which.site[1:M, u])
      
    } # YR
    
  } # U
  
})

# ______________________________________________________________________________
# 4. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # STOCHASTIC
  # OPEN
  phi_b0 = rlogis(1, 0, 1),
  phi_sd = rexp(1, 1),
  phi_b = rlogis(7, 0, 1), 
  
  # rho - per capita recruitment (log scale)
  rho_b0 = runif(1, log(0.001), log(7)),
  rho_sd = rexp(1, 1),
  rho_b = rnorm(6, 0, sd = 1),
  
  # cluster-specific random scaling factors
  phi_c = rnorm(4, 0, sd = 1),
  rho_c = rnorm(4, 0, sd = 1),
  
  # initial inclusion
  psi = runif(constant.list$U, 0, 1),
  psi.sex = runif(1, 0, 1),
  
  # states
  z = state.inits,       
  
  # CLOSED
  # ACs - we'll keep these pretty close to the trap array to start
  s = array(runif(constant.list$M * constant.list$YR, -50, 50),   
            dim = c(constant.list$M, 2, constant.list$YR)),
  
  # detection
  # baseline hazard
  lam0_b0 = rnorm(1, 0, 1),
  lam0_b = rnorm(3, 0, 1),
  
  # previous capture
  alpha2_b0 = rnorm(1, 0, 1),
  alpha2_b = rnorm(3, 0, 1),
  
  # spatial scale
  alpha1_b0 = rnorm(1, 0, 1),
  alpha1_b = rnorm(3, 0, 1),
  
  # latent covariates
  sex = ifelse(is.na(data.list$sex) == F, NA, 0)
  
)

# ______________________________________________________________________________
# 5. Parameters to monitor ----
# ______________________________________________________________________________

# for PPC: z, sex, AC coordinates
# can probably drop gamma, and psi after initial testing
# N.all can be split, and we won't need the available after tuning the augmented guys

monitor <- c(
  
  # open
  "phi_b0", "phi_sd", "phi_c", "phi_b", 
  "rho_b0", "rho_sd", "rho_c", "rho_b",
  "psi", 
  
  # detection
  "lam0_b0", "lam0_b", 
  "alpha2_b0", "alpha2_b", 
  "alpha1_b0", "alpha1_b", 
  
  # states and counts
  "N", "N.avail", "z"
  
)

# ______________________________________________________________________________
# 6. Set up model----
# ______________________________________________________________________________

model.1 <- nimbleModel(
  
  code = model.code,
  constants = constant.list,
  data = data.list,
  inits = inits,
  calculate = T
  
)  

# ______________________________________________________________________________
# 7. Compile functions and model ----
# ______________________________________________________________________________

# functions
cat_p_compiled <- compileNimble(cat_p)
bern_p_compiled <- compileNimble(bern_p)

# model
model.1.comp <- compileNimble(model.1)

# ______________________________________________________________________________
# 8. Configure on the uncompiled model ----
# ______________________________________________________________________________

# MCMC configuration
model.1.conf <- configureMCMC(model.1, monitors = monitor)

#______________________________________________________________________________
# 9. Add block samplers ----
# ______________________________________________________________________________
# 9a. Open parameters ----
# ______________________________________________________________________________

# remove samplers
model.1.conf$removeSamplers(
  
  c("phi_b0", "phi_sd", "phi_c[1]", "phi_c[2]", "phi_c[3]", "phi_c[4]",
    "rho_b0", "rho_sd", "rho_c[1]", "rho_c[2]", "rho_c[3]", "rho_c[4]")
)

# add samplers
# phi
model.1.conf$addSampler(
  
  c("phi_b0", "phi_sd", "phi_c[1]", "phi_c[2]", "phi_c[3]", "phi_c[4]"), 
  
  type = "RW_block",
  control = list(
    
    adaptInterval = 50,
    adaptScaleOnly = F
    
  )
  
)

# rho
model.1.conf$addSampler(
  
  c("rho_b0", "rho_sd", "rho_c[1]", "rho_c[2]", "rho_c[3]", "rho_c[4]"), 
  
  type = "RW_block",
  control = list(
    
    adaptInterval = 50,
    adaptScaleOnly = F
    
  )
  
)

# ______________________________________________________________________________
# 9b. Detection parameters ----
# ______________________________________________________________________________

# remove samplers
model.1.conf$removeSamplers(
  
  c("lam0_b0", "lam0_b[1]", 
    "alpha2_b0", "alpha2_b[1]",
    "alpha1_b0", "alpha1_b[1]")
)

# add samplers
model.1.conf$addSampler(
  
  c("lam0_b0", "lam0_b[1]", 
    "alpha2_b0", "alpha2_b[1]",
    "alpha1_b0", "alpha1_b[1]"), 
  
  type = "RW_block",
  control = list(
    
    adaptInterval = 50,
    adaptScaleOnly = F
    
  )
  
)

# ______________________________________________________________________________
# 10. Build MCMC from configuration ----
# ______________________________________________________________________________

# build MCMC
mcmc.1 <- buildMCMC(model.1.conf)

# compile again with the MCMC
mcmc.1.comp <- compileNimble(mcmc.1, project = model.1)

# ______________________________________________________________________________
# 9. Run sampling ----
# ______________________________________________________________________________

model.1.run <- runMCMC(
  
  mcmc = mcmc.1.comp,
  niter = 100,
  nburnin = 50,
  nchains = 1,
  samplesAsCodaMCMC = TRUE
  
)
