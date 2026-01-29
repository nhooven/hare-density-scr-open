# PROJECT: Open multi-session SCR
# SCRIPT: 01 - Data checking and formatting
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 12 Jan 2026
# COMPLETED: 29 Jan 2026
# LAST MODIFIED: 29 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________
# 2a. Open capture histories ----

# this dataset includes individual ID information, sites trapped, and
# known open population states for 3-4 years

# if individuals dispersed to different sites, they get multiple observations
# to preserve the multi-session structure (since this is rare)

# asterisks for state 2 during PRE1 are for trapped animals that won't be given
# states in the large model because we can't include the winter 2022-2023 trapping

# from this dataset, we want to keep MRID, site, sex, and each state

# ______________________________________________________________________________

open.ch <- read.csv("trapping_data/open_ch.csv")

# keep needed columns and rename to make everything easier
open.ch <- open.ch %>%
  
  dplyr::select(MRID,
                Site,
                Sex,
                PRE1,
                PRE2,
                POST1,
                POST2) %>%
  
  rename(y1 = PRE1,
         y2 = PRE2,
         y3 = POST1,
         y4 = POST2) %>%
  
  # replace "" with NA in y1
  mutate(y1 = ifelse(y1 == "", NA, y1))

# ______________________________________________________________________________
# 2b. Mark-recapture data ----

# these are J x K frames of MRIDs, that I will need to manipulate for SCR
# we'll extract both individual capture histories and a trap operation matrix from this
# we can also pull session-specific Ks from here

# ______________________________________________________________________________

mr.data.2022 <- read.csv("trapping_data/trapping_2022.csv")
mr.data.2023 <- read.csv("trapping_data/trapping_2023.csv")
mr.data.2024 <- read.csv("trapping_data/trapping_2024.csv")
mr.data.2025 <- read.csv("trapping_data/trapping_2025.csv")

# ______________________________________________________________________________
# 3. Site / session / treatment lookup table ----
# ______________________________________________________________________________

site.lookup <- data.frame(
  
  # integer session ID
  sessionID = 1:41,
  
  # number of secondary occasions
  K = c(8, 8, 6, 5, 7,
        6, 5, 6, 6, 6, 6, 6, 4, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6),
  
  # character site name
  SiteName = c(c("2A", "2B", "2C", "3A", "3C"),
               rep(c("1A", "1B", "1C", "2A", "2B", "2C",
                     "3A", "3B", "3C", "4A", "4B", "4C"),
                   times = 3)),
  
  # integer site ID
  siteID = c(c(4, 5, 6, 7, 9),
             rep(1:12,
                 times = 3)),
  
  # integer year YR
  year = c(rep(1, times = 5),
           rep(2, times = 12),
           rep(3, times = 12),
           rep(4, times = 12)),
  
  # treatment variable
  trt = c(rep("untreated", times = 17),
          rep(c("ret", "pil", "untreated",
                "pil", "ret", "untreated",
                "pil", "ret", "untreated",
                "ret", "pil", "untreated"),
              times = 2))

)

# write to .csv for reference
write.csv(site.lookup, "for_model/site_lookup.csv", row.names = F)

# ______________________________________________________________________________
# 4. Bind MR data together with unique session IDs ----
# ______________________________________________________________________________

# 2022
mr.data.2022.1 <- mr.data.2022 %>%
  
  # subset only sessions we want
  filter(Site %in% c("2A", "2B", "2C", "3A", "3C")) %>%
  
  # remove the last two columns
  dplyr::select(-c("D9", "D10")) %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 1) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D8)

# 2023
mr.data.2023.1 <- mr.data.2023 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 2) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D6) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# 2024
mr.data.2024.1 <- mr.data.2024 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 3) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D8) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# 2025
mr.data.2025.1 <- mr.data.2025 %>%
  
  # replace missing data explictly with NAs
  mutate_at(vars(matches("D")), ~ifelse(. == "", NA, .)) %>%
  
  # add year variable
  mutate(year = 4) %>%
  
  # rename Site
  rename(SiteName = Site) %>%
  
  # join in other index variables from site.lookip
  left_join(site.lookup, by = c("SiteName", "year")) %>%
  
  # rearrange sensibly
  dplyr::select(sessionID, siteID, year, SiteName, trt, Trap, D1:D6) %>%
  
  # arrange by sessionID
  arrange(sessionID)

# bind all together (should have 8 total occasion columns)
mr.data <- bind_rows(mr.data.2022.1, mr.data.2023.1, mr.data.2024.1, mr.data.2025.1)

# ______________________________________________________________________________
# 5. Check that MRIDs match exactly ----
# ______________________________________________________________________________
# 5a. Check function ----
# ______________________________________________________________________________

check_mrid <- function (focalID) {
  
  # DATA SETUP
  
  # extract correct indices
  focal.indices <- site.lookup %>% filter(sessionID == focalID)
  
  # keep correct column (based on primary occasion)
  if (focal.indices$year == 1) {open.ch.1 <- open.ch[ , c(1:3, 4)]}
  if (focal.indices$year == 2) {open.ch.1 <- open.ch[ , c(1:3, 5)]}
  if (focal.indices$year == 3) {open.ch.1 <- open.ch[ , c(1:3, 6)]}
  if (focal.indices$year == 4) {open.ch.1 <- open.ch[ , c(1:3, 7)]}
  
  # change name to "y" for generality
  names(open.ch.1)[4] <- "y"
  
  # subset focal indivs from open cH
  focal.indivs <- open.ch.1 %>%
    
    filter(Site == focal.indices$SiteName,
           y == 2) %>%
    
    # drop y column
    dplyr::select(-y)
  
  # subset just the secondary session CHs
  mr <- mr.data %>%
    
    filter(sessionID == focalID) %>%
    
    dplyr::select(D1:D8)
  
  # CHECK IDs
  
  # generate list of MRIDs from entire trapping data
  # bind together into one vector
  all.outcomes <- as.vector(as.matrix(mr))
  
  # remove trap indicators
  all.outcomes.unique <- unique(all.outcomes)
  
  trap.mrid <- all.outcomes.unique[! all.outcomes.unique %in% c("O", "B", "X", "C", "E", NA)]
  
  # compare to focal.indivs
  # which MRIDs are IN MR and NOTIN indivs?
  iMR.niIndivs <- trap.mrid[which(trap.mrid %notin% focal.indivs$MRID)]
  
  # which MRIDs are NOTIN MR and IN indivs?
  niMR.iIndivs <- focal.indivs$MRID[which(focal.indivs$MRID %notin% trap.mrid)]
  
  # dfs to hold information
  # in MR and not in Indivs
  if (length(iMR.niIndivs) > 0) {
    
    iMR.niIndivs.df <- data.frame(
      
      MRID = iMR.niIndivs,
      which.list = "iMR.niIndivs",
      where.else = ifelse(length(open.ch.1$y[which(open.ch.1$MRID %in% iMR.niIndivs)]) > 0,
                          open.ch.1$y[which(open.ch.1$MRID %in% iMR.niIndivs)],
                          NA)
      
    )
    
  } else {
    
    iMR.niIndivs.df <- 
      
      data.frame(
        
        MRID = NA,
        which.list = "iMR.niIndivs",
        where.else = NA
        
      )
    
  }
  
  # not in MR and in Indivs
  if (length(niMR.iIndivs) > 0) {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = niMR.iIndivs,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  } else {
    
    niMR.iIndivs.df <- data.frame(
      
      MRID = NA,
      which.list = "niMR.iIndivs",
      where.else = NA
      
    )
    
  }
  
  check.df <- rbind(iMR.niIndivs.df, niMR.iIndivs.df)
  
  # return
  return(check.df)
  
}

# ______________________________________________________________________________
# 5b. Change conditions ----

# if we have individuals captured originally in (an)other unit(s),
# change their MR entry(ies) to "C"
# this is to account for the trap being partially closed during the night.
# We'll keep those individuals with their original unit and assume any
# their AC is somewhere between the units.

# an alternative is to merge the grids (for 2A/2B, for example)
# but this complicates any inference on treatment effects
# so for simplicity, since grid-switching (within year) is rare,
# we'll ignore it for the most part

# I think it'll be better to do this manually
# we'll keep code for any we change

# NOTE: some individuals will show up in the open CH for time t because they were 
# in state 2 for t - 1 and t + 1

# ______________________________________________________________________________

# session 1
check_mrid(1)

mr.data[4, 12] <- "C"
mr.data[20, 10] <- "C"

# session 8
check_mrid(8)

mr.data[265, 10] <- "C"

# session 10
mr.data[339, 10] <- "C"

# session 33
mr.data[1159, 10] <- "C"
mr.data[1166, 11] <- "C"

# session 38
mr.data[1333, 10] <- "C"

# ______________________________________________________________________________
# 5c. Sanity check ----
# ______________________________________________________________________________

all.checks <- data.frame()

for (i in 1:41) {
  
  focal.check <- check_mrid(i)
  
  all.checks <- rbind(all.checks, focal.check)
  
}

# ______________________________________________________________________________
# 6. Reformat MR data ----

# this will result in an array with dimensions [i, max(K), t]
# we can no longer index everything by session, so the site.lookup will be our friend

# split open CH by MRID x Site
open.ch.1 <- open.ch %>% mutate(MRID.Site = paste0(MRID, ".", Site))
open.ch.split <- split(open.ch.1, f = ~ MRID.Site)

# ______________________________________________________________________________
# 6a. Function to find closed CHs from MR data ----
# ______________________________________________________________________________

extract_closed_CH <- function (x) {
  
  # this function accepts a list of open CHs split by MRID AND Site
  
  # the function will return a matrix of dimension [4, 8] for each indiv
  # we need this to be easily bindable into an array after the fact
  
  # primary occasion column names
  yr.names <- c("y1", "y2", "y3", "y4")
  
  # loop through all primary occasions, checking whether indiv has a closed CH
  # pack into a length 4 list
  indiv.ch <- matrix(NA, nrow = 4, ncol = max(site.lookup$K))
  
  for (i in 1:4) {
    
    # K occasions (must be a 0-8 integer)
    if (length(site.lookup$K[site.lookup$SiteName == x$Site &
                             site.lookup$year == i]) > 0) {
      
      focal.K <- site.lookup$K[site.lookup$SiteName == x$Site &
                                 site.lookup$year == i]
      
    } else {
      
      focal.K <- 0
      
    }
    
    # if the state is 1, 3, or NA, fill the matrix with J + 1
    # for non-recruited states, the CH won't matter because the lp will be zeroed
    # for latent states, the animal wasn't captured at all
    if (x[ , yr.names[i]] %in% c(1, 3, NA)) {
      
      # ensure that any extra occasions get NAs, for consistency
      indiv.ch[i, ] <- cbind(matrix(data = 37, 
                                    nrow = 1, 
                                    ncol = focal.K),
                             matrix(data = NA,
                                    nrow = 1,
                                    ncol = max(site.lookup$K) - focal.K))
      
      # else, the state is 2 and hare *likely* has a CH
      # this will not be the case for indivs with imputed states
      # which I will ensure this handles appropriately
    } else {
      
      # find correct MR data
      mr.data.focal <- mr.data %>%
        
        filter(SiteName == x$Site &
               year == i) %>%
        
        # keep only the MR data
        dplyr::select(D1:D8)
      
      # extract vector of MRIDs
      MRID.vec <- unique(as.vector(as.matrix(mr.data.focal)))
      
      # is focal indiv in this list?
      # if so, extract its CH
      if (x$MRID %in% MRID.vec) {
        
        # extract its CH, going by occasion k
        indiv.ch.vec <- vector()
        
        for (k in 1:focal.K) {
          
          indiv.ch.vec[k] <- ifelse(length(which(mr.data.focal[, k] == x$MRID)) > 0,
                                    which(mr.data.focal[, k] == x$MRID),
                                    37)
          
        }
        
        # and bind in extra NAs, if needed
        indiv.ch.vec <- c(indiv.ch.vec, rep(NA, times = max(site.lookup$K) - focal.K))
        
        # and add into list as a matrix (row vector)
        indiv.ch[i, ] <- matrix(indiv.ch.vec, nrow = 1)
        
        # else (not trapped) give a J + 1 filled matrix as before
      } else {
        
        indiv.ch[i, ] <- cbind(matrix(data = 37, 
                                      nrow = 1, 
                                      ncol = focal.K),
                               matrix(data = NA,
                                      nrow = 1,
                                      ncol = max(site.lookup$K) - focal.K))
        
      }
      
    }
    
  }
  
  # return matrix
  return(indiv.ch)
  
}

# ______________________________________________________________________________
# 6b. Apply function and bind ----
# ______________________________________________________________________________

all.indiv.ch <- lapply(open.ch.split, extract_closed_CH)

# bind correctly
all.indiv.ch.arr <- array(NA, dim = c(length(all.indiv.ch), 
                                      max(site.lookup$K),
                                      4))

for (i in 1:length(all.indiv.ch)) {
  
  all.indiv.ch.arr[i , , 1] <- all.indiv.ch[[i]][1, ]
  all.indiv.ch.arr[i , , 2] <- all.indiv.ch[[i]][2, ]
  all.indiv.ch.arr[i , , 3] <- all.indiv.ch[[i]][3, ]
  all.indiv.ch.arr[i , , 4] <- all.indiv.ch[[i]][4, ]
  
}

# ______________________________________________________________________________
# 7. Binary trap response covariate ----
# ______________________________________________________________________________

# blank array
prev.cap.arr <- array(NA, dim = c(nrow(all.indiv.ch.arr), 
                                  max(site.lookup$K),
                                  4))

# loop through indivs i
for (i in 1:nrow(all.indiv.ch.arr)) {
  
  # loop through primary occasions y
  for (y in 1:4) {
    
    # if indiv was caught at all
    if (any(all.indiv.ch.arr[i , , y] %in% 1:36)) {
      
      first.cap.index <- which(all.indiv.ch.arr[i , , y] != 37)[1]
      
      # create a binary vector
      prev.cap.vec <- c(rep(0, times = first.cap.index - 1),
                        rep(1, times = max(site.lookup$K) - (first.cap.index - 1)))
      
      # add NAs
      prev.cap.vec[which(is.na(all.indiv.ch.arr[i , , y]))] <- NA
      
      # bind in
      prev.cap.arr[i, , y] <- prev.cap.vec
      
      # else indiv not caught, all NAs
    } else {
      
      prev.cap.arr[i, , y] <- rep(NA, times = max(site.lookup$K))
      
    }
    
  }
  
}

# ______________________________________________________________________________
# 8. Individual indices and covariates ----

# these will be sessionID, siteID, clusterID, treatment, year (for the session)
# and sex (for the individual)

# we need these to be perfectly aligned to the individuals in the previous 2 arrays

# ______________________________________________________________________________

# sort open.ch.1 based upon the list order
identical(names(open.ch.split), sort(open.ch.1$MRID.Site))

# we just need to arrange open.ch.1 in ascending order
open.ch.2 <- open.ch.1[order(open.ch.1$MRID.Site), ]

identical(names(open.ch.split), open.ch.2$MRID.Site)

# now we can just tack on covariate values in that order
# this will end up being a list of each covariate
# each matrix will be [i] or [i, YR]

# ______________________________________________________________________________
# 8a. siteID ----
# ______________________________________________________________________________

# as integer
site.to.int <- left_join(data.frame("SiteName" = open.ch.2$Site),
                         site.lookup[6:17 , c("SiteName", "siteID")])

# should be the same site per year
cov.site <- site.to.int$siteID

# ______________________________________________________________________________
# 8b. clusterID ----
# ______________________________________________________________________________

# as integer
cov.cluster <- case_when(site.to.int$siteID %in% c(1:3) ~ 1,
                         site.to.int$siteID %in% c(4:6) ~ 2,
                         site.to.int$siteID %in% c(7:9) ~ 3,
                         site.to.int$siteID %in% c(10:12) ~ 4)

# ______________________________________________________________________________
# 8c. Treatment ----

# varies by year

# ______________________________________________________________________________

cov.ret <- matrix(NA, nrow = nrow(open.ch.2), ncol = 4)
cov.pil <- matrix(NA, nrow = nrow(open.ch.2), ncol = 4)

# add zeroes for pre-treatment
cov.ret[ , c(1:2)] <- 0
cov.pil[ , c(1:2)] <- 0

# post
cov.ret[which(cov.site[ , 1] %in% c(1, 5, 8, 10)), c(3:4)] <- 1
cov.ret[which(cov.site[ , 1] %notin% c(1, 5, 8, 10)), c(3:4)] <- 0

cov.pil[which(cov.site[ , 1] %in% c(2, 4, 7, 11)), c(3:4)] <- 1
cov.pil[which(cov.site[ , 1] %notin% c(2, 4, 7, 11)), c(3:4)] <- 0

# ______________________________________________________________________________
# 8d. Sex ----

# this must be binary (F == 0, M == 1) with NAs for anything latent

# ______________________________________________________________________________

cov.sex <- case_when(open.ch.2$Sex == "F" ~ 0,
                     open.ch.2$Sex == "M" ~ 1,
                     open.ch.2$Sex == "U" ~ NA)

# ______________________________________________________________________________
# 8e. indivID ----

# this will be helpful to have as just an index

# ______________________________________________________________________________

cov.indivID <- 1:nrow(open.ch.2)

# ______________________________________________________________________________
# 8f. Bind into a list ----
# ______________________________________________________________________________

indiv.covs <- list(cov.site, cov.cluster, cov.ret, cov.pil, cov.sex, cov.indivID)

# ______________________________________________________________________________
# 9. Trap operation matrices ----

# n.site matrices
# J x K x YR array
# The values in here will be multiplied within the model to allow for or "zero out"
# trap-specific capture probabilities
# Cheekily, I will include 0.5 as a possible value to encode our uncertainty about
# trap availability given the trap was:
# (1) found closed or 
# (2) caught a bycatch species

# we won't know how long the trap was available for a hare capture, so we consider
# its availability to be a Bernoulli trial. Maybe it was open all night and closed
# or caught another critter in the morning after no hare was going to go in. Maybe
# it closed or caught another critter shortly after we set it.

# In effect, the multiplication by 0.5 allows the full capture probability for
# half of all MCMC draws, and zeroes it out for the other half

# any "escaped" hares that did not receive an AnimalID will also induce a 0.5 here

# ______________________________________________________________________________

trap.op.list <- list()

# loop through site
for (i in 1:12) {
  
  trap.op.array <- array(data = NA, dim = c(36, max(site.lookup$K), 4))
  
  # loop through year
  for (y in 1:4) {
    
    # if there aren't any, give all NAs
    if (nrow(mr.data %>% filter(siteID == i & year == y)) == 0) {
      
      trap.op.array[ , , y] <- matrix(NA, nrow = 36, ncol = max(site.lookup$K))
      
      # else, fill according to the rules
    } else {
      
      # subset MR session
      focal.mr <- mr.data %>% 
        
        filter(siteID == i & 
               year == y) %>%
        
        # keep MR data
        dplyr::select(D1:D8)
        
      # unlist and assign
      focal.mr.vect <- unlist(focal.mr)
      
      focal.mr.vect[which(focal.mr.vect %notin% c("X", "C", "B", "E"))] <- 1
      focal.mr.vect[which(focal.mr.vect %in% c("X"))] <- 0
      focal.mr.vect[which(focal.mr.vect %in% c("C", "B", "E"))] <- 0.5
      
      trap.op.array[ , , y] <- matrix(as.numeric(focal.mr.vect),
                                      nrow = 36,
                                      ncol = max(site.lookup$K))
      
    }
    
  }
  
  # bind into list
  trap.op.list[[i]] <- trap.op.array
  
}

# ______________________________________________________________________________
# 10. Occasions by session ----

# U x YR

# ______________________________________________________________________________

occ.sess <- as.matrix(
  
  site.lookup %>% 
  
  dplyr::select(K, siteID, year) %>%
  
  pivot_wider(names_from = year,
              values_from = K) %>%
  
  arrange(siteID) %>%
  
  dplyr::select(-siteID)
  
  )

# ______________________________________________________________________________
# 11. Write to file ----
# ______________________________________________________________________________

# open CH
open.ch.3 <- as.matrix(open.ch.2[ , c(4:7)])

saveRDS(open.ch.3, "for_model/open_ch.rds")

# closed CH
saveRDS(all.indiv.ch.arr, "for_model/closed_ch.rds")

# previous capture
saveRDS(prev.cap.arr, "for_model/prev_cap.rds")

# individual covariates
saveRDS(indiv.covs, "for_model/indiv_covs.rds")

# trap operation
saveRDS(trap.op.list, "for_model/trap_op.rds")

# occasions by session
saveRDS(occ.sess, "for_model/occ_sess.rds")
