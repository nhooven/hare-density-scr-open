# PROJECT: Open multi-session SCR
# SCRIPT: 08 - CMR summaries
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 
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
# 3. Total captures ----
# ______________________________________________________________________________

sum(data.list$ch != 37, na.rm = T)
