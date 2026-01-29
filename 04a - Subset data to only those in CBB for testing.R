# PROJECT: Closed multi-session SCR
# SCRIPT: 04a - Subset data to only those in CBB for testing
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Jan 2026
# LAST MODIFIED: 29 Jan 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("for_model/constants.rds")
data.list <- readRDS("for_model/data.rds")

# ______________________________________________________________________________
# 3. CBB indices ----
# ______________________________________________________________________________

cbb.indices.M <- which(constant.list$cluster == 2)
cbb.indices.n <- cbb.indices.M[cbb.indices.M <= 876]

# ______________________________________________________________________________
# 3. Subset constants ----
# ______________________________________________________________________________

constants.cbb <- list(
  
  n = cbb.indices.n,
  M = length(cbb.indices.M),
  J = 36,
  U = 3,
  YR = 4,
  S.areas = constant.list$S.areas[4:6],
  K = constant.list$K[4:6, ],
  site = constant.list$site[cbb.indices.M]
  
  # no need for cluster or FT here
  
)

# replace site with 1:3
constants.cbb$site[constants.cbb$site == 4] <- 1
constants.cbb$site[constants.cbb$site == 5] <- 2
constants.cbb$site[constants.cbb$site == 6] <- 3

# ______________________________________________________________________________
# 4. Subset data ----
# ______________________________________________________________________________

data.cbb <- list(
  
  z = data.list$z[cbb.indices.M],
  ch = data.list$ch[cbb.indices.n],
  prev.cap = data.list$prev.cap[cbb.indices.M],
  trap.op = data.list$trap.op[ , , , 4:6],
  ret = data.list$ret[cbb.indices.M],
  pil = data.list$pil[cbb.indices.M],
  sex = data.list$sex[cbb.indices.M],
  zeroes = data.list$zeroes[cbb.indices.M, ]
  
)

# ______________________________________________________________________________
# 5. Write to file ----
# ______________________________________________________________________________

saveRDS(constants.cbb, "for_model/CBB_test/constants.rds")
saveRDS(data.cbb, "for_model/CBB_test/data.rds")

