# PROJECT: Open multi-session SCR
# SCRIPT: 06 - Check convergence
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 Apr 2026
# COMPLETED: 27 Apr 2026
# LAST MODIFIED: 27 Apr 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(MCMCvis)
library(coda)
library(bayestestR)

#_______________________________________________________________________________
# 2. Read in posterior samples ----
#_______________________________________________________________________________

# density and sex ratio
params.1 <- readRDS("final_samples/samples_clean_1.rds")

# open parameters
params.2 <- readRDS("final_samples/samples_clean_2.rds")

# closed parameters
params.3 <- readRDS("final_samples/samples_clean_3.rds")

#_______________________________________________________________________________
# 3. Summarize ----
#_______________________________________________________________________________
# 3a. Density ----

colnames(params.1[[1]])
ncol(params.1[[1]])

#_______________________________________________________________________________

MCMCsummary(params.1, params = colnames(params.1[[1]])[1:48], ISB = F, HPD = T, hpd_prob = 0.90)

#_______________________________________________________________________________
# 3b. Open parameters ----

colnames(params.2[[1]])

#_______________________________________________________________________________

# phi
MCMCsummary(params.2, params = colnames(params.2[[1]])[1:15], ISB = F, HPD = T, hpd_prob = 0.90)

# rho
MCMCsummary(params.2, params = colnames(params.2[[1]])[16:29], ISB = F, HPD = T, hpd_prob = 0.90)

#_______________________________________________________________________________
# 3c. Closed parameters ----

colnames(params.3[[1]])

#_______________________________________________________________________________

# lam0
MCMCsummary(params.3, params = colnames(params.3[[1]])[1:4], ISB = F, HPD = T, hpd_prob = 0.90)

# alpha2
MCMCsummary(params.3, params = colnames(params.3[[1]])[5:8], ISB = F, HPD = T, hpd_prob = 0.90)

# alpha1
MCMCsummary(params.3, params = colnames(params.3[[1]])[9:12], ISB = F, HPD = T, hpd_prob = 0.90)

#_______________________________________________________________________________
# 4. Trace ----
#_______________________________________________________________________________
# 4a. Density ----
#_______________________________________________________________________________

MCMCtrace(params.1, params = colnames(params.1[[1]])[1:48], ISB = F, pdf = F)

#_______________________________________________________________________________
# 4b. Open parameters ----
#_______________________________________________________________________________

# phi
MCMCtrace(params.2, params = colnames(params.2[[1]])[1:15], ISB = F, pdf = F)

# rho
MCMCtrace(params.2, params = colnames(params.2[[1]])[16:29], ISB = F, pdf = F)

#_______________________________________________________________________________
# 4c. Closed parameters ----
#_______________________________________________________________________________

# lam0
MCMCtrace(params.3, params = colnames(params.3[[1]])[1:4], ISB = F, pdf = F)

# alpha2
MCMCtrace(params.3, params = colnames(params.3[[1]])[5:8], ISB = F, pdf = F)

# alpha1
MCMCtrace(params.3, params = colnames(params.3[[1]])[9:12], ISB = F, pdf = F)
