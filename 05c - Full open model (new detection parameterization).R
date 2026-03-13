# PROJECT: Open multi-session SCR
# SCRIPT: 05c - Full open model (new detection parameterization)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Feb 2026
# COMPLETED: 
# LAST MODIFIED: 13 Mar 2026
# R VERSION: 4.4.3

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
    lam0_sd = double(0),
    lam0_c = double(0),
    lam0_b = double(1),
    
    # previous capture
    alpha2_b0 = double(0),
    alpha2_sd = double(0),
    alpha2_c = double(0),
    alpha2_b = double(1),
    
    # spatial scale of movement
    sigma_b0 = double(0),
    sigma_sd = double(0),
    sigma_c = double(0),
    sigma_b = double(1),
    
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
    alpha2 <- (alpha2_b0 + alpha2_sd * alpha2_c) + 
      
      alpha2_b[1] * sex + 
      alpha2_b[2] * ret + 
      alpha2_b[3] * pil
    
    # sigma and alpha1 - distance decay
    sigma <- exp(
      
      (sigma_b0 + sigma_sd * sigma_c) + 
                   
        sigma_b[1] * sex + 
        sigma_b[2] * ret + 
        sigma_b[3] * pil
      
      )
    
    alpha1 <- -1 / sigma
    
    # lam0 - baseline hazard of detection [i, k, t]
    lam0 <- exp(
      
      (lam0_b0 + lam0_sd * lam0_c) + 
        
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
      eta[j] <- lam0 * exp(alpha1 * d) *
        
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
    lam0_sd = double(0),
    lam0_c = double(0),
    lam0_b = double(1),
    
    # previous capture
    alpha2_b0 = double(0),
    alpha2_sd = double(0),
    alpha2_c = double(0),
    alpha2_b = double(1),
    
    # spatial scale of movement
    sigma_b0 = double(0),
    sigma_sd = double(0),
    sigma_c = double(0),
    sigma_b = double(1),
    
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
    alpha2 <- (alpha2_b0 + alpha2_sd * alpha2_c) + 
      
      alpha2_b[1] * sex + 
      alpha2_b[2] * ret + 
      alpha2_b[3] * pil
    
    # sigma and alpha1 - distance decay
    sigma <- exp(
      
      (sigma_b0 + sigma_sd * sigma_c) + 
        
        sigma_b[1] * sex + 
        sigma_b[2] * ret + 
        sigma_b[3] * pil
      
      )
    
    alpha1 <- -1 / sigma
    
    # loop over secondary occasions K
    for (k in 1:K) {
      
      eta.denom[k] <- 1.0
      
      # lam0 - baseline hazard of detection [i, k, t]
      lam0[k] <- exp(
        
        (lam0_b0 + lam0_sd * lam0_c) + 
          
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
        eta[k, j] <- lam0[k] * exp(alpha1 * d) *
          
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
  
  # demographic parameter priors
  # phi - persistence
  phi ~ dunif(0, 1)
  
  # rho - per capita recruitment
  rho ~ dunif(0, 7)
  
  # gamma - probability of available individual being recruited [u, t]
  for (u in 1:U) {
    
    # 2nd year (by site) to year 4
    for (t in (first.year[u] + 1):YR) {
      
      # logit constraint is important
      logit(gamma[u, t - 1]) <- log(N[u, t - 1] + 1e-6) + 
                                log(rho) - 
                                log(N.avail[u, t - 1] + 1e-6)
      
    } # YR
    
  } # U
  
  # omegas - transition matrices
  for (u in 1:U) {
    
    # psi - initial inclusion [U]
    psi[u] ~ dunif(0, 1)
    
    # omega1 - initial states [1:3, U]
    # this will be for the first year BY SITE
    omega1[1, u] <- 1 - psi[u]          # available with p = 1 - psi
    omega1[2, u] <- psi[u]              # recruited with p = psi
    omega1[3, u] <- 0.0                 # cannot start dead
    
    # omega - transition matrix [3 x 3, U, YR - 1]
    # rows: state at t
    # columns: state at t + 1
    for (t in (first.year[u]):(YR - 1)) {
      
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
  
  # open likelihood - by individual [M]
  for (i in 1:M) {
    
    # initial state (first.year by site)
    z[i, first.year[site[i]]] ~ dcat(omega1[1:3, site[i]])
    
    # subsequent states (YR > first.year)
    for (t in (first.year[site[i]]):(YR - 1)) {
      
      z[i, t + 1] ~ dcat(omega[z[i, t], 1:3, site[i], t])
      
    } # YR - 1
    
  } # M
  
  # ___________________________
  # CLOSED SCR SUB-MODEL
  # ___________________________
  
  # detection parameter priors
  # lam0 - baseline hazard detection (log scale)
  lam0_b0 ~ dnorm(0, sd = 1)  # mean 
  lam0_sd ~ T(dt(0, sigma = 1, df = 1), 0, )    # SD - half-Cauchy
  
  # alpha2 - trap response
  alpha2_b0 ~ dnorm(0, sd = 1)  # mean 
  alpha2_sd ~ T(dt(0, sigma = 0.5, df = 1), 0, )    # SD - half-Cauchy
  
  # sigma - spatial scale of movement
  sigma_b0 ~ dnorm(log(45), sd = 0.5)     # mean
  sigma_sd ~ T(dt(0, sigma = 1, df = 1), 0, )   # SD - half-Cauchy
  
  # detection linear coefficients (n = 3)
    # b1 - male effect
    # b2 - ret effect
    # b3 - pil effect
  for (x in 1:3) {
    
    lam0_b[x] ~ dnorm(0, sd = 1)  
    alpha2_b[x] ~ dnorm(0, sd = 1)  
    sigma_b[x] ~ dnorm(0, sd = 1)  
    
  }
  
  # c - cluster-specific random scaling factors
  for (c in 1:4) {
    
    lam0_c[c] ~ dnorm(0, sd = 1)
    alpha2_c[c] ~ dnorm(0, sd = 1)
    sigma_c[c] ~ dnorm(0, sd = 1)
    
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
          lam0_sd = lam0_sd,
          lam0_c = lam0_c[cluster[i]],
          lam0_b = lam0_b[1:3],
          
          # previous capture effect
          alpha2_b0 = alpha2_b0,
          alpha2_sd = alpha2_sd,
          alpha2_c = alpha2_c[cluster[i]],
          alpha2_b = alpha2_b[1:3],
          
          # spatial scale of movement
          sigma_b0 = sigma_b0,
          sigma_sd = sigma_sd,
          sigma_c = sigma_c[cluster[i]],
          sigma_b = sigma_b[1:3],
          
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
        lam0_sd = lam0_sd,
        lam0_c = lam0_c[cluster[i]],
        lam0_b = lam0_b[1:3],
        
        # previous capture effect
        alpha2_b0 = alpha2_b0,
        alpha2_sd = alpha2_sd,
        alpha2_c = alpha2_c[cluster[i]],
        alpha2_b = alpha2_b[1:3],
        
        # spatial scale of movement
        sigma_b0 = sigma_b0,
        sigma_sd = sigma_sd,
        sigma_c = sigma_c[cluster[i]],
        sigma_b = sigma_b[1:3],
        
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
  phi = runif(1, 0, 1),
  rho = runif(1, 0, 7),
  psi = runif(constant.list$U, 0, 1),
  psi.sex = runif(1, 0, 1),
  z = state.inits,       
  
  # CLOSED
  # ACs - we'll keep these pretty close to the trap array to start
  s = array(runif(constant.list$M * constant.list$YR, -50, 50),   
            dim = c(constant.list$M, 2, constant.list$YR)),
  
  # detection
  # baseline hazard
  lam0_b0 = rnorm(1, 0, 1),
  lam0_sd = runif(1, 0.1, 1),
  lam0_b = rnorm(3, 0, 1),
  
  # previous capture
  alpha2_b0 = rnorm(1, 0, 1),
  alpha2_sd = runif(1, 0.1, 1),
  alpha2_b = rnorm(3, 0, 1),
  
  # spatial scale
  sigma_b0 = rnorm(1, log(45), 0.5),
  sigma_sd = runif(1, 0.1, 1),
  sigma_b = rnorm(3, 0, 1),
  
  # c - random effect scaling parameters (keep close to zero to start)
  lam0_c = rnorm(4, 0, 0.25),
  alpha2_c = rnorm(4, 0, 0.25),
  sigma_c = rnorm(4, 0, 0.25),
  
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
  "psi", "phi", "rho",
  
  # detection
  "lam0_b0", "lam0_sd", "lam0_c", "lam0_b", 
  "alpha2_b0", "alpha2_sd", "alpha2_c", "alpha2_b", 
  "sigma_b0", "sigma_sd", "sigma_c", "sigma_b", 
  
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
# 9a. Detection parameters ----
# ______________________________________________________________________________

# matrices
# lam0 (lam0_b0, lam0_sd, lam0_c, lam0_b)
propCov.lam0 <- matrix(
  
  data = c(0.022,	-0.004,	-0.010,	-0.007,	-0.007,	-0.007,	-0.018,	-0.003,	-0.005,
           -0.004,	0.125,	0.040,	0.023,	0.035,	0.035,	0.011,	0.015,	0.006,
           -0.010,	0.040,	0.021,	0.012,	0.017,	0.015,	0.007,	0.004,	-0.001,
           -0.007,	0.023,	0.012,	0.012,	0.011,	0.011,	0.000,	-0.004,	-0.003,
           -0.007,	0.035,	0.017,	0.011,	0.018,	0.013,	0.002,	0.000,	-0.001,
           -0.007,	0.035,	0.015,	0.011,	0.013,	0.019,	0.004,	-0.001,	-0.002,
           -0.018,	0.011,	0.007,	0.000,	0.002,	0.004,	0.046,	0.008,	0.007,
           -0.003,	0.015,	0.004,	-0.004,	0.000,	-0.001,	0.008,	0.045,	0.009,
           -0.005,	0.006,	-0.001,	-0.003,	-0.001,	-0.002,	0.007,	0.009,	0.060
  ),
  
  nrow = 9,
  ncol = 9,
  byrow = T
  
)

# alpha2 (alpha2_b0, alpha2_c, alpha2_b[1-2])
propCov.alpha2 <- matrix(
  
  data = c(0.023,		-0.055,	-0.022,	-0.027,	-0.024,	-0.016,	-0.010,	
           -0.055,		0.638,	0.175,	0.125,	0.033,	0.024,	0.017,	
           -0.022,		0.175,	0.581,	0.048,	0.030,	-0.010,	0.001,	
           -0.027,		0.125,	0.048,	0.637,	0.070,	0.009,	-0.011,	
           -0.024,		0.033,	0.030,	0.070,	0.746,	-0.017,	-0.010,	
           -0.016,		0.024,	-0.010,	0.009,	-0.017,	0.035,	-0.001,	
           -0.010,		0.017,	0.001,	-0.011,	-0.010,	-0.001,	0.042
  ),
  
  nrow = 7,
  ncol = 7,
  byrow = T
  
)

# sigma (sigma_b0, sigma_sd, sigma_c)
propCov.sigma <- matrix(
  
  data = c(0.018,	0.006,	-0.052,	-0.032,	-0.049,	-0.032,
           0.006,	0.017,	-0.020,	0.019,	-0.030,	0.030,
           -0.052,	-0.020,	0.210,	0.097,	0.166,	0.100,
           -0.032,	0.019,	0.097,	0.158,	0.065,	0.146,
           -0.049,	-0.030,	0.166,	0.065,	0.269,	0.057,
           -0.032,	0.030,	0.100,	0.146,	0.057,	0.235
  ),
  
  nrow = 6,
  ncol = 6,
  byrow = T
  
)

# check positive-definiteness
eigen(propCov.lam0)$values
eigen(propCov.alpha2)$values
eigen(propCov.sigma)$values

# check condition
kappa(propCov.lam0)
kappa(propCov.alpha2)
kappa(propCov.sigma)

# remove samplers
model.1.conf$removeSamplers(
  
  c("lam0_b0", "lam0_sd", "lam0_c[1]", "lam0_c[2]", "lam0_c[3]", "lam0_c[4]", "lam0_b[1]", "lam0_b[2]", "lam0_b[3]",
    "alpha2_b0", "alpha2_c[1]", "alpha2_c[2]", "alpha2_c[3]", "alpha2_c[4]", "alpha2_b[1]", "alpha2_b[2]",
    "sigma_b0", "sigma_sd", "sigma_c[1]", "sigma_c[2]", "sigma_c[3]", "sigma_c[4]")
  
)

# add samplers
# lam0
model.1.conf$addSampler(
  
  c("lam0_b0", "lam0_sd", "lam0_c[1]", "lam0_c[2]", "lam0_c[3]", "lam0_c[4]", "lam0_b[1]", "lam0_b[2]", "lam0_b[3]"), 
  
  type = "RW_block",
  control = list(
    
    "propCov" = propCov.lam0,
    adaptInterval = 50,
    adaptScaleOnly = F
    
  )
  
)

# alpha2
model.1.conf$addSampler(
  
  c("alpha2_b0", "alpha2_c[1]", "alpha2_c[2]", "alpha2_c[3]", "alpha2_c[4]", "alpha2_b[1]", "alpha2_b[2]"), 
  
  type = "RW_block",
  control = list(
    
    "propCov" = propCov.alpha2,
    adaptInterval = 50,
    adaptScaleOnly = F
    
  )
  
)

# sigma
model.1.conf$addSampler(
  
  c("sigma_b0", "sigma_sd", "sigma_c[1]", "sigma_c[2]", "sigma_c[3]", "sigma_c[4]"), 
  
  type = "RW_block",
  control = list(
    
    "propCov" = propCov.sigma,
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
  niter = 10000,
  nburnin = 5000,
  nchains = 1,
  samplesAsCodaMCMC = TRUE
  
)