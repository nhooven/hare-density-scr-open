# PROJECT: Closed multi-session SCR
# SCRIPT: 04 - Prepare all data for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Jan 2026
# LAST MODIFIED: 02 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# I would like all of this to be done in a separate script so all I need to 
# feed to Kamiak is the modeling script and the constant and data lists

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# capture histories
open.ch <- readRDS("for_model/open_ch.rds")
closed.ch <- readRDS("for_model/closed_ch.rds")

# previous capture
prev.cap <- readRDS("for_model/prev_cap.rds")

# covariates
indiv.covs <- readRDS("for_model/indiv_covs.rds")

# occasions by session
occ.sess <- readRDS("for_model/occ_sess.rds")

# trap operation matrix
trap.op <- readRDS("for_model/trap_op.rds")

# trap coords
trap.coords <- readRDS("for_model/trap_coords.rds")

# S limits
S.lim <- readRDS("for_model/S_lim.rds")

# S areas
S.area <- readRDS("for_model/S_area.rds")

# ______________________________________________________________________________
# 3. Data cleaning ----
# ______________________________________________________________________________
# 3a. Trap operation ----

# this is a list of arrays
# can we have a 4-D array?
# apparently we can...
# I believe nimble can only handle 4

# ______________________________________________________________________________
  
trap.op.4D <- array(data = NA, dim = c(nrow(trap.op[[1]]),
                                       ncol(trap.op[[1]]),
                                       4,
                                       12))

for (i in 1:12) {
  
  trap.op.4D[ , , , i] <- trap.op[[i]]
  
}

# ______________________________________________________________________________
# 3b. Open CH ----

# we should replace any 2* with NA and ensure this is an integer matrix

# ______________________________________________________________________________

open.ch[open.ch == "2*"] <- NA

open.ch.1 <- matrix(as.integer(open.ch),
                    nrow = nrow(closed.ch),
                    ncol = 4)

# ______________________________________________________________________________
# 4. Define required quantities ----
# ______________________________________________________________________________

# number of sites U
U = 12

# number of captured individuals 
n.u <- table(indiv.covs[[1]])     # by site
n = sum(n.u)                      # total

# number of traps J
J = 36

# S limits [U, 2, 2]
S.lim

# trap coords [J, 2, U]
trap.coords

# S areas (in ha) [U]
S.areas <- S.area$area

# ______________________________________________________________________________
# 5. Data augmentation ----
# ______________________________________________________________________________

# define how many individuals to add by site
# ultimately this should be related to the total number we ever saw in a site
# maybe start with 5 x n.u?

# totals by site
n.aug.u <- n.u * 5

n.aug = sum(n.aug.u)

# total individuals in the dataset
M = n + n.aug

# open CH states to add
open.ch.aug <- matrix(NA, nrow = n.aug, ncol = 4)

open.ch.all <- rbind(open.ch.1, open.ch.aug)

# zeroes matrix (Chandler method) [M, yr]
# this can be by open CH for simplicity
zeroes <- matrix(data = cbind(c(rep(NA, times = n), rep(0, times = n.aug)),
                              c(rep(NA, times = n), rep(0, times = n.aug)),
                              c(rep(NA, times = n), rep(0, times = n.aug)),
                              c(rep(NA, times = n), rep(0, times = n.aug))),
                 ncol = 4)

# add zero rows to prev.cap
prev.cap.1 <- array(data = NA, dim = c(M, 8, 4))

# loop through years
for (y in 1:4) {
  
  # bind together
  prev.cap.1[ , , y] <- rbind(prev.cap[ , , y],
                              matrix(data = 0,
                                     nrow = n.aug,
                                     ncol = 8))
  
}  

# covariate assignment
# siteID
aug.siteID <- rep(1:12, times = n.aug.u)
  
# clusterID
aug.clusterID <- case_when(aug.siteID %in% c(1:3) ~ 1,
                           aug.siteID %in% c(4:6) ~ 2,
                           aug.siteID %in% c(7:9) ~ 3,
                           aug.siteID %in% c(10:12) ~ 4)

# sex
aug.sex <- rep(NA, times = n.aug)

# treatment
aug.ret <- matrix(NA, nrow = n.aug, ncol = 4)
aug.pil <- matrix(NA, nrow = n.aug, ncol = 4)

# add zeroes for pre-treatment
aug.ret[ , c(1:2)] <- 0
aug.pil[ , c(1:2)] <- 0

# post
aug.ret[which(aug.siteID %in% c(1, 5, 8, 10)), c(3:4)] <- 1
aug.ret[which(aug.siteID %notin% c(1, 5, 8, 10)), c(3:4)] <- 0

aug.pil[which(aug.siteID %in% c(2, 4, 7, 11)), c(3:4)] <- 1
aug.pil[which(aug.siteID %notin% c(2, 4, 7, 11)), c(3:4)] <- 0

# bind all indiv covs together
indiv.covs.M <- list(
  
  # site
  c(indiv.covs[[1]], aug.siteID),
  
  # cluster
  c(indiv.covs[[2]], aug.clusterID),
  
  # ret
  rbind(indiv.covs[[3]], aug.ret),
  
  # pil
  rbind(indiv.covs[[4]], aug.pil),
  
  # sex
  c(indiv.covs[[5]], aug.sex),
  
  # indivID
  c(indiv.covs[[6]], (length(indiv.covs[[6]]) + 1):M)
  
)

# ______________________________________________________________________________
# 6. Additional covariates / constants ----
# ______________________________________________________________________________
# 6a. Forest type covariate ----

# 0 for SFL, 1 for XMC

# ______________________________________________________________________________

indiv.covs.M[[7]] <- ifelse(indiv.covs.M[[2]] == 4, 1, 0)

# ______________________________________________________________________________
# 6b. "First year" index for each site ----
# ______________________________________________________________________________

first.year <- c(2, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2)

# ______________________________________________________________________________
# 7. Build lists for model ----
# ______________________________________________________________________________
# 7a. Constants ----
# ______________________________________________________________________________  

constant.list <- list(
  
  # scalar constants (for loop indices)
  n = n,
  M = M,
  J = J,
  U = 12,
  YR = 4,
  
  # site-specific constants [U]
  S.areas = S.areas,
  
  # K/session lookup [U, YR]
  K = occ.sess,
  
  # "first" year by site [U]
  first.year = first.year,
  
  # indices for each individual [M]
  site = indiv.covs.M[[1]],
  cluster = indiv.covs.M[[2]],
  ft = indiv.covs.M[[7]]
  
)

# ______________________________________________________________________________
# 7b. Data ----
# ______________________________________________________________________________  
  
data.list <- list(
  
  # open capture histories (states - z) [M, YR]
  z = open.ch.all,
  
  # closed capture histories [n, max(K), YR]
  ch = closed.ch,
  
  # previous capture [M, max(K), YR]
  prev.cap = prev.cap.1,
  
  # trap operation [J, max(K), YR, U]
  trap.op = trap.op.4D,
  
  # state space limits [U, 2, 2]
  S.lim = S.lim,
  
  # trap coordinates [J, 2, U]
  trap.coords = trap.coords,
  
  # individual data [M]
  ret = indiv.covs.M[[3]],
  pil = indiv.covs.M[[4]],
  sex = indiv.covs.M[[5]],
  
  # data augmentation [M, YR]
  zeroes = zeroes
  
)

# ______________________________________________________________________________
# 7c. Sensible state inits ----

# we need to pass state inits to nimble so it can freely estimate latent states
# these should be flexible so each chain can get its own starts
# but NOT flexible enough to propose impossible states

# this function will "march" along each indivs open CH, referencing a simple
# transition matrix

# ______________________________________________________________________________  

# function
make_init_states <- function (x) {
  
  # transition matrix 3 x 3
  trans.mat <- matrix(
    
    data = c(0.5, 0.5, 0,
             0, 0.5, 0.5,
             0, 0, 1),
    nrow = 3,
    ncol = 3,
    byrow = T
    
  )
  
  na.indices <- which(is.na(x))
  
  # how many NAs?
  n.na <- length(na.indices)
  
  # if there are no integer states
  if (n.na == 4) {
    
    # vector to hold all integer states
    integer.states <- vector()
    
    integer.states[1] <- rbinom(1, 1, 0.5) + 1
    
    # for next values
    for (k in 2:4) {
      
      # only proceed if k is in na.indices
      if (k %in% na.indices) {
        
        # look at k - 1
        trans.probs.k <- trans.mat[integer.states[k - 1], ]
        
        integer.states[k] <-  which(rmultinom(1, 1, trans.probs.k) == 1)
        
      }
      
    }
    
  }
  
  # else, assign x to integer.states
  else {
    
    integer.states <- x
    
    # first index is NA and second IS NOT 3
    if (k == 1 & 
        k %in% na.indices &
        integer.states[2] != 3) {
      
      integer.states[1] <- rbinom(1, 1, 0.5) + 1
      
    }
    
    # if first index is NA and second IS 3
    if (k == 1 & 
        k %in% na.indices &
        integer.states[2] == 3) {
      
      integer.states[1] <- 2
      
    }
    
    # for next values
    for (k in 2:4) {
      
      # only proceed if k is in na.indices
      if (k %in% na.indices) {
        
        # look at k - 1
        trans.probs.k <- trans.mat[integer.states[k - 1], ]
        
        integer.states[k] <-  which(rmultinom(1, 1, trans.probs.k) == 1)
        
      }
      
    }
    
  }
  
  # convert true states to NA
  y <- integer.states
  
  y[which(is.na(x) == F)] <- NA
  
  return(y)
  
}

# apply function thrice (one for each chain)
state.inits <- list()

state.inits[[1]] <- t(apply(open.ch.all, 1, make_init_states))
state.inits[[2]] <- t(apply(open.ch.all, 1, make_init_states))
state.inits[[3]] <- t(apply(open.ch.all, 1, make_init_states))

# ______________________________________________________________________________
# 8. Save to file ----
# ______________________________________________________________________________    

saveRDS(constant.list, "for_model/constants.rds")
saveRDS(data.list, "for_model/data.rds")
saveRDS(state.inits, "for_model/state_inits.rds")
