# PROJECT: Open multi-session SCR
# SCRIPT: 04 - Prepare all data for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Jan 2026
# LAST MODIFIED: 25 Feb 2026
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

# trap deaths
trap.deaths <- readRDS("for_model/trap_deaths.rds")

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
S.areas <- S.area$area.rect

# ______________________________________________________________________________
# 5. Data augmentation ----
# ______________________________________________________________________________

# define how many individuals to add by site
# ultimately this should be related to the total number we ever saw in a site
# from the CBB-only model, we still had plenty to play with by the end
# let's do x 5

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

# add zero rows to prev.cap and one rows to trap.deaths
prev.cap.1 <- array(data = NA, dim = c(M, 8, 4))
trap.deaths.1 <- array(data = NA, dim = c(M, 8, 4))

# loop through years
for (y in 1:4) {
  
  # bind together
  prev.cap.1[ , , y] <- rbind(prev.cap[ , , y],
                              matrix(data = 0,
                                     nrow = n.aug,
                                     ncol = 8))
  
  trap.deaths.1[ , , y] <- rbind(trap.deaths[ , , y],
                                 matrix(data = 1,
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
# 6c. Year and site treatment variables ----
# ______________________________________________________________________________

# year treatment [YR]
post1 <- c(0, 0, 1, 0)
post2 <- c(0, 0, 0, 1)

# site ret [U]
site.ret <- c(1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0)

# site pil [U]
site.pil <- c(0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0)

# site cluster
site.clust <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3))

# ______________________________________________________________________________
# 6d. Site indicator for N counting function ----
# ______________________________________________________________________________

# which site [M, U]
which.site <- matrix(0, nrow = M, ncol = U)

for (i in 1:M) {
  
  which.site[i, indiv.covs.M[[1]][i]] <- 1
  
}

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
  
  # year-specific constants [YR]
  post1 = post1,
  post2 = post2,
  
  # site-specific constants [U]
  S.areas = S.areas,
  site.ret = site.ret,
  site.pil = site.pil,
  site.clust = site.clust,
  
  # K/session lookup [U, YR]
  K = occ.sess,
  
  # "first" year by site [U]
  first.year = first.year,
  
  # indices for each individual [M]
  site = indiv.covs.M[[1]],
  cluster = indiv.covs.M[[2]],
  ft = indiv.covs.M[[7]],
  
  # previous capture [M, max(K), YR]
  prev.cap = prev.cap.1,
  
  # trap deaths [M, max(K), YR]
  trap.deaths = trap.deaths.1,
  
  # trap operation [J, max(K), YR, U]
  trap.op = trap.op.4D,
  
  # state space limits [U, 2, 2]
  S.lim = S.lim,
  
  # trap coordinates [J, 2, U]
  trap.coords = trap.coords,
  
  # site indicator
  which.site = which.site
  
)

# ______________________________________________________________________________
# 7b. Data ----
# ______________________________________________________________________________  
  
data.list <- list(
  
  # open capture histories (states - z) [M, YR]
  z = open.ch.all,
  
  # closed capture histories [n, max(K), YR]
  ch = closed.ch,
  
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
  
  # which indices to fill?
  na.indices <- which(is.na(x))
  
  # how many NAs?
  n.na = length(na.indices)
  
  # ALL STATES ARE OBSERVED
  # z inits get NAs
  if (n.na == 0) {
    
    y <- rep(NA, 4)
    
  }
  
  # ALL STATES ARE UNOBSERVED
  # march along, filling in possible states
  if (n.na == 4) {
    
    # first one can be either a 1 or 2
    x.1 <- vector()  # vector to hold all integer states
    
    x.1[1] <- rbinom(1, 1, 0.5) + 1
    
    # loop through the rest
    for (t in 2:4) {
      
      if (t == 2) {
        
        x.1[2] <- case_when(
          
          # t1 == 1
          x.1[1] == 1 & x.1[3] == 1 ~ 1,                             # 1
          x.1[1] == 1 & x.1[3] == 2 ~ rbinom(1, 1, 0.5) + 1,         # 1 or 2
          x.1[1] == 1 & x.1[3] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
          x.1[1] == 1 & is.na(x.1[3]) == T ~ rbinom(1, 1, 0.5) + 1,  # 1 or 2
          
          # t1 == 2
          x.1[1] == 2 & x.1[3] == 2 ~ 2,                             # 2
          x.1[1] == 2 & x.1[3] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
          x.1[1] == 2 & is.na(x.1[3]) == T ~ rbinom(1, 1, 0.5) + 2   # 2 or 3
          
        )
        
        # t == 3
      } else if (t == 3) {
        
        x.1[3] <- case_when(
          
          # t2 == 1
          x.1[2] == 1 & x.1[4] == 2 ~ rbinom(1, 1, 0.5) + 1,         # 1 or 2
          x.1[2] == 1 & x.1[4] == 3 ~ 2,                             # 2
          x.1[2] == 1 & is.na(x.1[4]) == T ~ rbinom(1, 1, 0.5) + 1,  # 1 or 2
          
          # t2 == 2
          x.1[2] == 2 & x.1[4] == 2 ~ 2,                             # 2
          x.1[2] == 2 & x.1[4] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
          x.1[2] == 2 & is.na(x.1[4]) == T ~ rbinom(1, 1, 0.5) + 2,  # 2 or 3
          
          # t2 == 3
          x.1[2] == 3 ~ 3                                            # 3
          
        )
        
        # t == 4
      } else {
        
        # only need to check x.1[3]
        x.1[4] <- case_when(
          
          x.1[3] == 1 ~ 2,
          x.1[3] == 2 ~ rbinom(1, 1, 0.5) + 2,
          x.1[3] == 3 ~ 3
          
        )
        
      } 
      
    } # t
    
    #change x.1 to y by adding NAs
    y <- x.1
    
    y[which(is.na(x) == F)] <- NA
    
  }
  
  # SOME STATES ARE OBSERVED, SOME ARE UNOBSERVED
  if (n.na %in% c(1:3)) {
    
    x.1 <- x
    
    # loop through years, checking previous AND next state
    for (t in 1:4) {
      
      # proceed if the current state is NA
      if (is.na(x.1[t]) == T) {
        
        # t == 1
        if (t == 1) {
          
          # this state can only be 1 or 2
          x.1[1] <- rbinom(1, 1, 0.5) + 1
          
          # t == 2
        } else if (t == 2) {
          
          x.1[2] <- case_when(
            
            # t1 == 1
            x.1[1] == 1 & x.1[3] == 1 ~ 1,                             # 1
            x.1[1] == 1 & x.1[3] == 2 ~ rbinom(1, 1, 0.5) + 1,         # 1 or 2
            x.1[1] == 1 & x.1[3] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
            x.1[1] == 1 & is.na(x.1[3]) == T ~ rbinom(1, 1, 0.5) + 1,  # 1 or 2
            
            # t1 == 2
            x.1[1] == 2 & x.1[3] == 2 ~ 2,                             # 2
            x.1[1] == 2 & x.1[3] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
            x.1[1] == 2 & is.na(x.1[3]) == T & x.1[4] != 2 ~ rbinom(1, 1, 0.5) + 2,   # 2 or 3
            x.1[1] == 2 & is.na(x.1[3]) == T & x.1[4] == 2 ~ 2        # 2
            
          )
          
          # t == 3
        } else if (t == 3) {
          
          x.1[3] <- case_when(
            
            # t2 == 1
            x.1[2] == 1 & x.1[4] == 2 ~ rbinom(1, 1, 0.5) + 1,         # 1 or 2
            x.1[2] == 1 & x.1[4] == 3 ~ 2,                             # 2
            x.1[2] == 1 & is.na(x.1[4]) == T ~ rbinom(1, 1, 0.5) + 1,  # 1 or 2
            
            # t2 == 2
            x.1[2] == 2 & x.1[4] == 2 ~ 2,                             # 2
            x.1[2] == 2 & x.1[4] == 3 ~ rbinom(1, 1, 0.5) + 2,         # 2 or 3
            x.1[2] == 2 & is.na(x.1[4]) == T ~ rbinom(1, 1, 0.5) + 2,  # 2 or 3
            
            # t2 == 3
            x.1[2] == 3 ~ 3                                            # 3
            
          )
          
          # t == 4
        } else {
          
          # only need to check x.1[3]
          x.1[4] <- case_when(
            
            x.1[3] == 1 ~ 2,
            x.1[3] == 2 ~ rbinom(1, 1, 0.5) + 2,
            x.1[3] == 3 ~ 3
            
          )
          
        } 
        
      }
      
    } # t
    
    # change x.1 to y by adding NAs
    y <- x.1
    
    y[which(is.na(x) == F)] <- NA
    
  }
  
  return(y)
  
}

# apply function thrice (one for each chain)
state.inits <- list()

state.inits[[1]] <- t(apply(open.ch.all, 1, make_init_states))
state.inits[[2]] <- t(apply(open.ch.all, 1, make_init_states))
state.inits[[3]] <- t(apply(open.ch.all, 1, make_init_states))

# check validity
# merge states function
merge_states <- function (open.ch.all, state.inits) {
  
  all.states <- data.frame()
  
  for (i in 1:nrow(open.ch.all)) {
    
    focal.ch <- open.ch.all[i, ]
    focal.inits <- state.inits[i, ]
    
    # bind together
    both <- rbind(focal.ch, focal.inits)
    
    # extract (I think this works)
    real.states <- both[which(is.na(both) == F)]
    
    # bind in
    all.states <- rbind(all.states, real.states)
    
  }
  
  return(all.states)
  
}

merge.1 <- merge_states(open.ch.all, state.inits = state.inits[[1]])
merge.2 <- merge_states(open.ch.all, state.inits = state.inits[[2]])
merge.3 <- merge_states(open.ch.all, state.inits = state.inits[[3]])

sum(apply(merge.1, 1, is.unsorted))
sum(apply(merge.2, 1, is.unsorted))
sum(apply(merge.3, 1, is.unsorted))

# ______________________________________________________________________________
# 8. Check likelihood-breaking in the closed CH ----

# We cannot have individuals captured in a trap that has its categorical
# probability "zeroed out" by the trap operation/trap death variable

# I'll write a function that checks this for every closed CH
# and returns a list of problematic individuals

# ______________________________________________________________________________    

check_chs <- function (data.list, constant.list) {
  
  ch <- data.list$ch
  trap.op <- data.list$trap.op
  trap.deaths <- data.list$trap.deaths
  site <- constant.list$site
  
  # subset to only captured individuals
  trap.deaths.1 <- trap.deaths[1:nrow(ch), ,]
  
  # blank df
  problem.ch.all <- data.frame()
  
  # loop through indivs i
  for (i in 1:nrow(ch)) {
    
    # subset trap.op to correct site
    trap.op.focal <- trap.op[ , , , site[i]]
    
    # loop through year t
    for (t in 1:4) {
      
      # focal closed session 1:8
      ch.focal <- ch[i, , t]
      
      # focal trap deaths 1:8
      trap.deaths.focal <- trap.deaths.1[i, , t]
      
      # focal trap op 1:36, 1:8
      trap.op.focal <- trap.op[ , , t, site[i]]
      
      # check ch against trap deaths
      # IF hare was caught at all (ch.focal != 37)
      if (any(ch.focal != 37) & any(trap.deaths.focal == 0)) {
        
        #  then trap.deaths must be 1
        if (any(which(ch.focal != 37) %in% which(trap.deaths.focal == 0))) {
          
          # bind in
          problem.ch <- data.frame(i = i,
                                   t = t,
                                   prob = "trap.deaths")
          
          problem.ch.all <- rbind(problem.ch.all, problem.ch)
          
        }
        
      }
      
      # check ch against trap op
      if (any(ch.focal != 37, na.rm = T)) {
        
        # trap(s) and occasions
        k.all <- which(ch.focal != 37)
        
        # loop through
        for (k in k.all) {
          
          j <- ch.focal[k]
          
          if (0 %in% trap.op.focal[j, k]) {
            
            # bind in
            problem.ch <- data.frame(i = i,
                                     t = t,
                                     prob = "trap.op")
            
            problem.ch.all <- rbind(problem.ch.all, problem.ch)
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  return(problem.ch.all)
  
}

# use function
(check.chs <- check_chs(data.list, constant.list))

# ______________________________________________________________________________
# 9. Save to file ----
# ______________________________________________________________________________    

saveRDS(constant.list, "for_model/constants.rds")
saveRDS(data.list, "for_model/data.rds")
saveRDS(state.inits, "for_model/state_inits.rds")
