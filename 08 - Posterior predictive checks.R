# PROJECT: Open multi-session SCR
# SCRIPT: 08 - Posterior predictive checks
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 04 May 2026
# COMPLETED: 06 May 2026
# LAST MODIFIED: 06 May 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(nimble)
library(coda)
library(tidybayes)
library(abind)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

constant.list <- readRDS("for_model/constants.rds")
data.list <- readRDS("for_model/data.rds")

# samples for PPC
samples <- readRDS("final_samples/samples_clean_4.rds")

# ______________________________________________________________________________
# 3. Clean samples ----
# ______________________________________________________________________________

# create tibble
samples.1 <- samples |>
  
  # bind
  do.call(rbind, args = _) |>
  
  # data frame |> tibble
  as.data.frame() |> tibble() |>
  
  # keep every fifth row
  slice(
    
    sample(seq(1, n(), by = 5))
    
  )

# remove huge list
rm(samples)

# prepare for tidying - remove spaces from colnames
colnames(samples.1) <- gsub(" ", "", colnames(samples.1))

# tidy
# s
samples.s <- samples.1 |>
  
  posterior::as_draws_df() |>

  spread_draws(s[i, j, k]) |>
  
  # drop sampling indices
  dplyr::select(-c(.chain, .iteration)) |>
  
  # rename
  # i is the individual, j is the coordinate, k is the year
  rename(
    
    indiv = i,
    coord = j,
    year = k
    
  )

# z
samples.z <- samples.1 |>
  
  posterior::as_draws_df() |>
  
  spread_draws(z[i, j]) |>
  
  # drop sampling indices (keep draw)
  dplyr::select(-c(.chain, .iteration)) |>
  
  # rename
  # i is the individual, j is the year
  rename(
    
    indiv = i,
    year = j,
    
  )

# detection
samples.detec <- samples.1 |>
  
  dplyr::select(
    
    lam0_b0, `lam0_b[1]`, `lam0_b[2]`, `lam0_b[3]`,
    alpha2_b0, `alpha2_b[1]`, `alpha2_b[2]`, `alpha2_b[3]`,
    alpha1_b0, `alpha1_b[1]`, `alpha1_b[2]`, `alpha1_b[3]`
    
  ) |>
  
  # rename
  rename(
    
    lam0_b1 = `lam0_b[1]`, 
    lam0_b2 = `lam0_b[2]`, 
    lam0_b3 = `lam0_b[3]`,
    alpha2_b1 = `alpha2_b[1]`, 
    alpha2_b2 = `alpha2_b[2]`, 
    alpha2_b3 = `alpha2_b[3]`,
    alpha1_b1 = `alpha1_b[1]`, 
    alpha1_b2 = `alpha1_b[2]`, 
    alpha1_b3 = `alpha1_b[3]`
    
  ) |>
  
  posterior::as_draws_df() |>
  
  # drop sampling indices (keep draw)
  dplyr::select(-c(.chain, .iteration))

# detection
samples.sex <- samples.1 |>
  
  posterior::as_draws_df() |>
  
  spread_draws(sex[i]) |>
  
  # drop sampling indices (keep draw)
  dplyr::select(-c(.chain, .iteration)) %>%
  
  # rename
  rename(indiv = i)

# ______________________________________________________________________________
# 4. Set up ch.array ----
# ______________________________________________________________________________

# set up ch.array
# expand true capture histories [M, K, J, YR]
ch.array <- array(data = NA, dim = c(constant.list$M,
                                     8,
                                     36,
                                     4))


for (i in 1:nrow(data.list$ch)) {
  
  for (t in 1:4) {
    
    # only proceed if we trapped there
    if (constant.list$first.year[constant.list$site[i]] <= t) {
      
      for (k in 1:constant.list$K[constant.list$site[i], t]) {
        
        # which cap
        binary.vect <- rep(0, times = 36)
        
        if (data.list$ch[i, k, t] == 37) {
          
          ch.array[i, k, , t] <- binary.vect
          
        } else {
          
          binary.vect[data.list$ch[i, k, t]] <- 1
          
          ch.array[i, k, , t] <- binary.vect
          
        }
        
      } # K
      
    }
    
  } # YR
  
} # N

# ______________________________________________________________________________
# 5. Helper functions ----
# ______________________________________________________________________________
# 5a. Calculate probabilities (nimbleFunction) ----
# ______________________________________________________________________________

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

# ______________________________________________________________________________
# 5b. PPC function ----
# ______________________________________________________________________________

# we'll split by site-year so trap frequencies make sense
site.year.df <- data.frame(
  
  site = c(which(constant.list$first.year == 1),
           1:12,
           1:12,
           1:12),
  
  year = c(rep(1, 5), rep(2, 12), rep(3, 12), rep(4, 12))
  
) |>
  
  mutate(site.year = 1:n())

# split
site.year.split <- split(site.year.df, ~ site.year)

# function
# this returns the discrepancies for each draw
calc_discrep <- function (x) {
  
  # x input is a site.year
  focal.site <- x$site
  focal.year <- x$year
  
  focal.K <- unname(constant.list$K[focal.site, focal.year])
  
  # which individuals?
  site.year.indivs <- focal.samples.z %>%
    
    filter(
      
      indiv <= 842 &                                        # must be real
      indiv %in% which(constant.list$site == focal.site) &  # must be at the site
      year == focal.year &                                  # must be during the year
      z == 2                                                # must be present
      
    )
  
  # loop over individuals
  # initialize arrays
  p <- array(data = 0, dim = c(length(unique(site.year.indivs$indiv)),
                               8,
                               37))
  
  sim.ch <- array(data = 0, dim = c(length(unique(site.year.indivs$indiv)),
                                    8,
                                    36))
  
  for (i in 1:length(unique(site.year.indivs$indiv))) {
    
    focal.indiv = unique(site.year.indivs$indiv)[i]
    
    # activity center
    s.indiv <- c(
      
      focal.samples.s$s[focal.samples.s$indiv == focal.indiv &
                          focal.samples.s$coord == 1 &
                          focal.samples.s$year == focal.year],
      
      focal.samples.s$s[focal.samples.s$indiv == focal.indiv &
                          focal.samples.s$coord == 2 &
                          focal.samples.s$year == focal.year]
      
    )
    
    # loop over occasions
    prev.c <- 0 # initialize at zero
    
    for (k in 1:focal.K) {
      
      # switch prev.c to 1 if the animal was previously captured
      if (k > 1 & sim.ch.k != 37) { prev.c <- 1 }
      
      # calculate probabilities
      p[i, k, 1:37] <- cat_p(
        
        # baseline hazard of detection
        lam0_b0 = focal.samples.detec$lam0_b0,
        lam0_b = unlist(focal.samples.detec[c("lam0_b1", "lam0_b2", "lam0_b3")]),
        
        # previous capture effect
        alpha2_b0 = focal.samples.detec$alpha2_b0,
        alpha2_b = unlist(focal.samples.detec[c("alpha2_b1", "alpha2_b2", "alpha2_b3")]),
        
        # spatial scale of movement
        alpha1_b0 = focal.samples.detec$alpha1_b0,
        alpha1_b = unlist(focal.samples.detec[c("alpha1_b1", "alpha1_b2", "alpha1_b3")]),
        
        # constants
        sex = focal.samples.sex$sex[focal.samples.sex$indiv == focal.indiv],   # 0 == F
        ret = constant.list$det.ret[focal.indiv, focal.year],
        pil = constant.list$det.pil[focal.indiv, focal.year],
        z = 1,  # must be 1, not 2!
        prev.cap = prev.c,
        trap.deaths = constant.list$trap.deaths[focal.indiv, k, focal.year],  
        trap.op = constant.list$trap.op[ , k, focal.year, focal.site],
        J = 36,
        trap.coords = constant.list$trap.coords[ , 1:2, focal.site],
        
        # activity center coords
        s = s.indiv
        
      )
      
      # simulate capture histories
      # these will be categorical - we need to expand into a binary array
      sim.ch.k <- rcat(prob = p[i, k, ])
      
      # IF capped, flip to 1
      if (sim.ch.k != 37) { sim.ch[i, k, sim.ch.k] <- 1 }
      
    } # K 
    
  } # i
  
  # aggregate capture histories for test statistics
  # subset true ch [M, K, J, T]
  true.ch <- ch.array[unique(site.year.indivs$indiv), 1:8, , focal.year]
  
  # subset p for E.count
  p.ch <- p[ , , 1:36]
  
  # calculate discrepancy by site
  # FT1 - individual frequencies (along rows, individual counts)
  FT.indiv <- data.frame(
    
    # simulated counts
    count.sim = apply(sim.ch, 1, sum, na.rm = T),
    
    # real counts
    count.real = apply(true.ch, 1, sum, na.rm = T),
    
    # expected counts (must remove the "no capture" bin)
    E.count = apply(p.ch, 1, sum, na.rm = T)
    
  ) %>%
    
    mutate(
      
      discrep.sim = (sqrt(count.sim) - sqrt(E.count))^2,
      discrep.real = (sqrt(count.real) - sqrt(E.count))^2
      
    ) %>%
    
    summarize(
      
      discrep.sim.sum = sum(discrep.sim),
      discrep.real.sum = sum(discrep.real)
      
    )
  
  # FT2 - trap frequencies (along slices, trap counts)
  FT.trap <- data.frame(
    
    # simulated counts
    count.sim = apply(sim.ch, 3, sum, na.rm = T),
    
    # real counts
    count.real = apply(true.ch, 3, sum, na.rm = T),
    
    # expected counts (must remove the "no capture" bin)
    E.count = apply(p.ch, 3, sum, na.rm = T)
    
  ) %>%
    
    mutate(
      
      discrep.sim = (sqrt(count.sim) - sqrt(E.count))^2,
      discrep.real = (sqrt(count.real) - sqrt(E.count))^2
      
    ) %>%
    
    summarize(
      
      discrep.sim.sum = sum(discrep.sim),
      discrep.real.sum = sum(discrep.real)
      
    )
  
  # FT3 - individual x trap frequencies (along individual and trap)
  FT.indiv.trap <- list(
    
    # simulated counts
    count.sim = apply(sim.ch, c(1, 3), sum, na.rm = T),
    
    # real counts
    count.real = apply(true.ch, c(1, 3), sum, na.rm = T),
    
    # expected counts (must remove the "no capture" bin)
    E.count = apply(p.ch, c(1, 3), sum, na.rm = T)
    
  )
  
  # calculate by indiv-trap
  discrep.sim <- matrix(data = NA, nrow = nrow(sim.ch), ncol = 36)
  discrep.real <- matrix(data = NA, nrow = nrow(sim.ch), ncol = 36)
  
  for (i in 1:nrow(sim.ch)) {
    
    for (j in 1:36) {
      
      discrep.sim[i, j] <- (sqrt(FT.indiv.trap$count.sim[i, j]) - 
                              sqrt(FT.indiv.trap$E.count[i, j]))^2
      
      discrep.real[i, j] <- (sqrt(FT.indiv.trap$count.real[i, j]) - 
                               sqrt(FT.indiv.trap$E.count[i, j]))^2
      
    } # j
    
  } # i
  
  FT.indiv.trap.df <- data.frame(discrep.sim.sum = sum(discrep.sim),
                                 discrep.real.sum = sum(discrep.real))
  
  # bind all together
  FT.all <- data.frame(
    
    rbind(
      
      FT.indiv,
      FT.trap,
      FT.indiv.trap.df
      
    ),
    
    cbind("check" = c("FT_indiv", "FT_trap", "FT_indiv_trap"))
    
  )  |> 
    
    # add site
    mutate(site = focal.site,
           year = focal.year)
  
  return(FT.all)
  
} # f() 
  
# ______________________________________________________________________________
# 6. Goodness of fit tests ----

total.draws <- max(samples.s$.draw)

# ______________________________________________________________________________

# loop through iterations
FT.all.n <- data.frame()

# start time
start.time <- Sys.time()

for (n in 1:total.draws) {
  
  # subset
  focal.samples.s <- samples.s |> filter(.draw == n)
  focal.samples.z <- samples.z |> filter(.draw == n)
  focal.samples.detec <- samples.detec |> filter(.draw == n)
  focal.samples.sex <- samples.sex |> filter(.draw == n)
  
  # apply function (~15-20 s)
  FT.all.sites <- lapply(site.year.split, calc_discrep)
  
  # bind and add draw
  FT.all.sites.df <- do.call(rbind, FT.all.sites) |>
    
    mutate(.draw = n)
  
  # bind in
  FT.all.n <- rbind(FT.all.n, FT.all.sites.df)
  
  # print status
  if (n %% 100 == 0) {
    
    elapsed.time <- round(as.numeric(difftime(Sys.time(), 
                                              start.time, 
                                              units = "mins")), 
                          digits = 1)
    
    print(paste0("Completed draw ", 
                 n, 
                 " of ", 
                 total.draws, 
                 " - ", 
                 elapsed.time, 
                 " mins"))
    
  } # IF
  
} # n

# save
saveRDS(FT.all.n, "ppc/FT_all_n.rds")

# ______________________________________________________________________________
# 7. Posterior predictive p-values ----

# proportion sim discrep > obs discrep

# ______________________________________________________________________________

pppval.session <- FT.all.n |>
  
  mutate(Dsim.Dobs = discrep.sim.sum > discrep.real.sum) |>
  
  group_by(check, site, year) |>
  
  summarize(
    
    p = sum(Dsim.Dobs) / n()
    
  ) |>
  
  # factor labels
  mutate(check = factor(check,
                        levels = c("FT_indiv",
                                   "FT_trap",
                                   "FT_indiv_trap"),
                        labels = c("Individual encounters",
                                   "Trap encounters",
                                   "Individual by trap encounters")))

pppval.session

# histograms of pppvals
ggplot(pppval.session) +
  
  theme_bw() +
  
  facet_wrap(~ check) +
  
  geom_histogram(aes(x = p),
                 color = "white",
                 fill = "aquamarine4",
                 bins = 20) +
  
  # vlines
  geom_vline(xintercept = c(0.1, 0.9),
             linetype = "dashed",
             color = "gray25") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA)) +
  
  xlab("Proportion simulated > observed discrepancy") +
  
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
  coord_cartesian(ylim = c(0.39, 8))

# 589 x 251

# summarize
pppval.session |>
  
  ungroup() |>
  
  group_by(check) |>
  
  summarize(
    
    med = median(p),
    min = min(p),
    max = max(p),
    n.poor.fit = sum(p > 0.9) + sum(p < 0.1)
    
  )

  