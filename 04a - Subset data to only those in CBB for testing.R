# PROJECT: Closed multi-session SCR
# SCRIPT: 04a - Subset data to only those in CBB for testing
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Jan 2026
# LAST MODIFIED: 16 Feb 2026
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
state.inits <- readRDS("for_model/state_inits.rds")

# ______________________________________________________________________________
# 3. CBB indices ----
# ______________________________________________________________________________

cbb.indices.M <- which(constant.list$cluster == 2)
cbb.indices.n <- cbb.indices.M[cbb.indices.M <= 876]

# ______________________________________________________________________________
# 3. Subset constants ----
# ______________________________________________________________________________

constants.cbb <- list(
  
  n = length(cbb.indices.n),
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
  
  z = data.list$z[cbb.indices.M, ],
  ch = data.list$ch[cbb.indices.n, , ],
  prev.cap = data.list$prev.cap[cbb.indices.M, , ],
  trap.deaths = data.list$trap.deaths[cbb.indices.M, , ],
  trap.op = data.list$trap.op[ , , , 4:6],
  S.lim = data.list$S.lim[4:6, , ],
  trap.coords = data.list$trap.coords[ , , 4:6],
  ret = data.list$ret[cbb.indices.M, ],
  pil = data.list$pil[cbb.indices.M, ], 
  sex = data.list$sex[cbb.indices.M],
  zeroes = data.list$zeroes[cbb.indices.M, ]
  
)

# ______________________________________________________________________________
# 5. Subset state inits ----
# ______________________________________________________________________________

state.inits.cbb <- list(
  
  state.inits[[1]][cbb.indices.M, ],
  state.inits[[2]][cbb.indices.M, ],
  state.inits[[3]][cbb.indices.M, ]
  
)

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
# 6. Write to file ----
# ______________________________________________________________________________

saveRDS(constants.cbb, "for_model/CBB_test/constants.rds")
saveRDS(data.cbb, "for_model/CBB_test/data.rds")
saveRDS(state.inits.cbb, "for_model/CBB_test/state_inits.rds")
