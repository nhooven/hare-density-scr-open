# PROJECT: Closed multi-session SCR
# SCRIPT: 06 - Re-factorization tests
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Feb 2026
# COMPLETED: 
# LAST MODIFIED: 14 Feb 2026
# R VERSION: 4.4.3

# Here we want to see if we can speed up the closed part of the model
# by re-factoring

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(nimble)
library(tictoc)
library(MCMCvis)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("for_model/CBB_test/constants.rds")
data.list <- readRDS("for_model/CBB_test/data.rds")

# ______________________________________________________________________________
# 3. Subset data ----

# here we only want individuals we caught in 2022 from 2A and 2B
# we'll do our own augmentation afterward

# ______________________________________________________________________________

# site indices (1 and 2)
site.indices <- which(constant.list$site %in% c(1, 2))

# which individuals were caught in year 1?
z.1 <- which(data.list$z[ ,1] == 2)

# which individuals were in 2A/2B in year 1?
indivs <- z.1[which(z.1 %in% site.indices)]

# augment
# indivs per site
n.site <- c(sum(constant.list$site[indivs] == 1),
            sum(constant.list$site[indivs] == 2))

n.aug <- n.site * 5

# ______________________________________________________________________________
# 4. New lists ----
# ______________________________________________________________________________

constant.list.1 <- list(
  
  n = length(indivs),
  M = length(indivs) + sum(n.aug),
  J = 36,
  U = 2,
  # no need for YR
  # no need for S.areas,
  K = constant.list$K[1:2, 1],
  site = c(constant.list$site[indivs], rep(1, n.aug[1]), rep(2, n.aug[2]))  
  
)

data.list.1 <- list(
  
  z = c(rep(1, length(indivs)), rep(NA, sum(n.aug))),
  ch = data.list$ch[indivs, , 1],
  prev.cap = rbind(data.list$prev.cap[indivs, , 1],
                   matrix(data = 0,
                          nrow = sum(n.aug),
                          ncol = 8)),
  trap.op = data.list$trap.op[ , , 1, 1:2],
  S.lim = data.list$S.lim[1:2, 1:2, 1:2],
  trap.coords = data.list$trap.coords[ , 1:2, 1:2],
  sex = c(data.list$sex[indivs], rep(NA, sum(n.aug))),
  zeroes = c(rep(NA, length(indivs)), rep(0, sum(n.aug)))
  
)

# ______________________________________________________________________________
# 4. Model code ----
# ______________________________________________________________________________
# 4a. MODEL 1 - BASELINE ----

# this model use the basic parameterization that I've been using
# we'll keep linear predictors for sex

# ______________________________________________________________________________

model.1.code <- nimbleCode({
  
  # detection
  # alpha0 - baseline detection
  alpha0_b0 ~ dunif(-10, 1)      # intercept
  alpha0_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # alpha2 - trap response
  alpha2_b0 ~ dnorm(0, sd = 1)  # intercept
  alpha2_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # sigma - spatial scale of movement
  sigma_b0 ~ dunif(log(10), log(70))      # intercept
  sigma_b1 ~ dnorm(0, sd = 1)             # effect of male
  
  for (u in 1:U) {
    
    psi[u] ~ dunif(0, 1)
    
  }
  
  psi.sex ~ dunif(0, 1)
  
  # calculate probabilities for closed likelihood
  # loop over individuals [M]
  for (i in 1:M) {
    
    # inclusion
    z[i] ~ dbern(psi[site[i]])
    
    # sex indicator (every latent individual gets one value)
    sex[i] ~ dbern(psi.sex)
    
    # detection
    alpha2[i] <- alpha2_b0 + alpha2_b1 * sex[i]
    log(sigma[i]) <- sigma_b0 + sigma_b1 * sex[i]
    alpha1[i] <- -1 / sigma[i]
    
    # s - activity centers [M, 2]
    s[i, 1] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
    s[i, 2] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
    
    # vectorized trap calculations
    # d - distances between s and each trap j
    d[i, 1:J] <- sqrt(pow(s[i, 1] - trap.coords[1:J, 1, site[i]], 2) + 
                      pow(s[i, 2] - trap.coords[1:J, 2, site[i]], 2))
      
    # g - distance kernel
    g[i, 1:J] <- exp(alpha1[i] * d[i, 1:J])
    
    # loop over secondary occasions [K] (indexed by site)
    for (k in 1:K[site[i]]) {
      
      # alpha0 - un-normalized baseline detection probability [i, k]
      logit(alpha0[i, k]) <- alpha0_b0 + alpha0_b1 * sex[i] +
                             alpha2[i] * prev.cap[i, k]
      
      # and loop over traps [J]
      for (j in 1:J) {
        
        # eta - un-normalized probability 
        eta[i, k, j] <- alpha0[i, k] * g[i, j] *
          
          # inclusion
          z[i] *
          
          # trap operation
          trap.op[j, k, site[i]]
        
        # p - normalized probabilities for categorical likelihood
        p[i, k, j] <- eta[i, k, j] / (1 + sum(eta[i, k, 1:J]))  # sum over all traps
        
      } # J
      
      # probability of not being captured as the complement of all trap-specific probs
      p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
      
    } # K
    
  } # M
  
  # closed likelihood for observed individuals
  for (i in 1:n) {
    
    # loop over occasions
    for (k in 1:K[site[i]]) {
      
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
      
    } # K
    
  } # n
  
  # closed likelihood for non-detected individuals
  for (i in (n + 1):M) {
    
    zeroes[i] ~ dbern(1 - prod(1 - p[i, 1:K[site[i]], 1:J]))
    
  } # (n + 1):M
  
  # DERIVED PARAMETERS
  # N
  # loop over sites
  for (u in 1:U) {
    
    # sex-specific N
    # total recruited individuals
    Nf[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 0))
    Nm[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 1))
    
    # total available individuals
    n.availf[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 0))
    n.availm[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 1))
    
    # pooled
    # total recruited individuals
    N[u] <- Nf[u] + Nm[u]
    
    # total available individuals
    n.avail[u] <- n.availf[u] + n.availm[u]
    
    # R - sex ratio (proportion male)
    R[u] <- Nm[u] / N[u]
    
  } # U
  
})

# ______________________________________________________________________________
# 4b. MODEL 2 - exp(alpha0) ----

# this avoids the logit transform
# really the p0 isn't a probability anyway

# ______________________________________________________________________________

model.2.code <- nimbleCode({
  
  # detection
  # alpha0 - baseline detection
  alpha0_b0 ~ dunif(-10, 1)      # intercept
  alpha0_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # alpha2 - trap response
  alpha2_b0 ~ dnorm(0, sd = 1)  # intercept
  alpha2_b1 ~ dnorm(0, sd = 1)  # effect of male
  
  # sigma - spatial scale of movement
  sigma_b0 ~ dunif(log(10), log(70))      # intercept
  sigma_b1 ~ dnorm(0, sd = 1)             # effect of male
  
  psi ~ dunif(0, 1)
  psi.sex ~ dunif(0, 1)
  
  # calculate probabilities for closed likelihood
  # loop over individuals [M]
  for (i in 1:M) {
    
    # inclusion
    z[i] ~ dbern(psi)
    
    # sex indicator (every latent individual gets one value)
    sex[i] ~ dbern(psi.sex)
    
    # detection
    alpha2[i] <- alpha2_b0 + alpha2_b1 * sex[i]
    log(sigma[i]) <- sigma_b0 + sigma_b1 * sex[i]
    alpha1[i] <- -1 / sigma[i]
      
    # s - activity centers [M, 2, YR]
    s[i, 1] ~ dunif(S.lim[site[i], 1, 1], S.lim[site[i], 2, 1])
    s[i, 2] ~ dunif(S.lim[site[i], 1, 2], S.lim[site[i], 2, 2])
      
    # vectorized trap calculations
    # d - distances between s and each trap j
    d[i, 1:J] <- sqrt(pow(s[i, 1] - trap.coords[1:J, 1, site[i]], 2) + 
                      pow(s[i, 2] - trap.coords[1:J, 2, site[i]], 2))
      
    # g - distance kernel
    g[i, 1:J] <- exp(alpha1[i] * d[i, 1:J])
      
      # loop over secondary occasions [K] (indexed by site)
      for (k in 1:K[site[i]]) {
        
      # alpha0 - un-normalized baseline detection probability [i, k, t]
      log(alpha0[i, k]) <- alpha0_b0 + alpha0_b1 * sex[i] + 
                           alpha2[i] * prev.cap[i, k]
        
        # and loop over traps [J]
        for (j in 1:J) {
          
          # eta - un-normalized probability 
          eta[i, k, j] <- alpha0[i, k] * g[i, j] *
            
            # inclusion
            z[i] *
            
            # trap operation
            trap.op[j, k, site[i]]
          
          # p - normalized probabilities for categorical likelihood
          p[i, k, j] <- eta[i, k, j] / (1 + sum(eta[i, k, 1:J]))  # sum over all traps
          
        } # J
        
        # probability of not being captured as the complement of all trap-specific probs
        p[i, k, J + 1] <- 1 - sum(p[i, k, 1:J])
        
      } # K
    
  } # M
  
  # closed likelihood for observed individuals
  for (i in 1:n) {
      
    # loop over occasions
    for (k in 1:K[site[i]]) {
        
      ch[i, k] ~ dcat(p[i, k, 1:(J + 1)])
        
    } # K
    
  } # n
  
  # closed likelihood for non-detected individuals
  for (i in (n + 1):M) {
      
      zeroes[i] ~ dbern(1 - prod(1 - p[i, 1:K[site[i]], 1:J]))
    
  } # (n + 1):M
  
  # DERIVED PARAMETERS
  # N
  # loop over sites
  for (u in 1:U) {
      
      # sex-specific N
      # total recruited individuals
      Nf[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 0))
      Nm[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 1))
      
      # total available individuals
      n.availf[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 0))
      n.availm[u] <- sum((z[1:M] == 1) * (site[1:M] == u) * (sex[1:M] == 1))
      
      # pooled
      # total recruited individuals
      N[u] <- Nf[u] + Nm[u]
      
      # total available individuals
      n.avail[u] <- n.availf[u] + n.availm[u]
      
      # R - sex ratio (proportion male)
      R[u] <- Nm[u] / N[u]
    
  } # U
  
})

# ______________________________________________________________________________
# 5. Shared inputs ----
# ______________________________________________________________________________
# 5a. Inits ----
# ______________________________________________________________________________

inits <- list(
  
  z = c(rep(NA, constant.list.1$n), rep(0, constant.list.1$M - constant.list.1$n)),
  
  # CLOSED
  psi = runif(2, 0, 1),
  s = cbind(runif(constant.list.1$M, -200, 200),
            runif(constant.list.1$M, -200, 200)),
  alpha0_b0 = runif(1, -10, 1),
  alpha0_b1 = rnorm(1, 0, 1),
  alpha2_b0 = rnorm(1, 0, 1),
  alpha2_b1 = rnorm(1, 0, 1),
  sigma_b0 = runif(1, log(10), log(70)),
  sigma_b1 = rnorm(1, 0, 1),
  sex = ifelse(is.na(data.list.1$sex) == F, NA, 0)
  
)

# ______________________________________________________________________________
# 5b. Monitors ----
# ______________________________________________________________________________

monitor <- c(
  
  "alpha0_b0", "alpha0_b1", "alpha2_b0", "alpha2_b1", "sigma_b0", "sigma_b1",
  "psi", "N"
  
)

# ______________________________________________________________________________
# 6. Model setup ----
# ______________________________________________________________________________
# 6a. MODEL 1 ----
# ______________________________________________________________________________

model.1 <- nimbleModel(
  
  code = model.1.code,
  constants = constant.list.1,
  data = data.list.1,
  inits = inits,
  calculate = F
  
)

model.1.conf <- configureMCMC(model.1, monitors = monitor)

# build model
model.1.mcmc <- buildMCMC(conf = model.1.conf)

# compile model
compileNimble(model.1)
model.1.comp <- compileNimble(model.1.mcmc, project = model.1)

# ______________________________________________________________________________
# 6b. MODEL 2 ----
# ______________________________________________________________________________

model.2 <- nimbleModel(
  
  code = model.2.code,
  constants = constant.list.1,
  data = data.list.1,
  inits = inits,
  calculate = F
  
)

model.2.conf <- configureMCMC(model.2, monitors = monitor)

# build model
model.2.mcmc <- buildMCMC(conf = model.2.conf)

# compile model
compileNimble(model.2)
model.2.comp <- compileNimble(model.2.mcmc, project = model.2)

# ______________________________________________________________________________
# 7. Time models ----

# we could time the whole process, but let's just do the MCMC now

# ______________________________________________________________________________

# model 1
tic()
model.1.run <- runMCMC(
  
  mcmc = model.1.comp,
  niter = 10000,
  nburnin = 5000,
  nchains = 1,
  samplesAsCodaMCMC = TRUE
  
)
toc()

# 724 s

# model 2
tic()
model.2.run <- runMCMC(
  
  mcmc = model.2.comp,
  niter = 10000,
  nburnin = 5000,
  nchains = 1,
  samplesAsCodaMCMC = TRUE
  
)
toc()

# 694 s

# ______________________________________________________________________________
# 8. Examine sampling ----
# ______________________________________________________________________________
# 8a. alpha0_b0 ----
# ______________________________________________________________________________

MCMCtrace(model.1.run, pdf = F, params = c("alpha0_b0"))    # both sampled poorly
MCMCtrace(model.2.run, pdf = F, params = c("alpha0_b0")) 

# ______________________________________________________________________________
# 8b. alpha2_b0 ----
# ______________________________________________________________________________

MCMCtrace(model.1.run, pdf = F, params = c("alpha2_b0"))    # both sampled poorly
MCMCtrace(model.2.run, pdf = F, params = c("alpha2_b0"))    

# ______________________________________________________________________________
# 8c. sigma_b0 ----
# ______________________________________________________________________________

MCMCtrace(model.1.run, pdf = F, params = c("sigma_b0"))    # both sampled poorly
MCMCtrace(model.2.run, pdf = F, params = c("sigma_b0"))    

# ______________________________________________________________________________
# 8d. N ----

# maxed these out. should add more augmented

# ______________________________________________________________________________

MCMCtrace(model.1.run, pdf = F, params = c("N"))    
MCMCtrace(model.2.run, pdf = F, params = c("N"))    # sampled much worse

MCMCsummary(model.1.run, params = c("n.avail"))    
MCMCtrace(model.2.run, pdf = F, params = c("n.avail"))

# 02-14-2026 - solved the alpha2 issue
# now I can actually use these models to improve sampling speed and efficiency