# PROJECT: Open multi-session SCR
# SCRIPT: 07c - Visualization (closed parameters)
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
library(MCMCvis)
library(coda)
library(bayestestR)
library(ggridges)        # ridgeline plot

# ______________________________________________________________________________
# 2. Read in samples ----
# ______________________________________________________________________________

params.3 <- readRDS("final_samples/samples_clean_3.rds")

# ______________________________________________________________________________
# 3. Parameter estimates ----

# these will be very similar to the hazard model coefficients

# HELPER FUNCTIONS
# long draws (for densities)
long_draws <- function (x) {
  
  x.1 <- x |>
    
    as.data.frame() %>%
    
    pivot_longer(cols = colnames(x))
  
  return(x.1)
  
}

# function to extract medians and CIs
extract_med_ci <- function (x, ci.1 = 0.50, ci.2 = 0.90) {
  
  x.1 <- x %>%
    
    group_by(name) %>%
    
    summarize(
      
      med = median(value),
      sd = sd(value),
      lo.1 = as.numeric(hdi(value, ci = ci.1)[2]),
      up.1 = as.numeric(hdi(value, ci = ci.1)[3]),
      lo.2 = as.numeric(hdi(value, ci = ci.2)[2]),
      up.2 = as.numeric(hdi(value, ci = ci.2)[3])
      
    ) %>%
    
    ungroup()
  
  return (x.1)
  
}

# ______________________________________________________________________________
# 3a. Process ----
# ______________________________________________________________________________

params.long <- params.3 |>
  
  do.call(rbind, args = _) |>
  
  long_draws()

# extract summaries
params.summary <- extract_med_ci(params.long)

# prepare for plotting
params.long.plot <- params.long %>%
  
  mutate(
    
    param = factor(
      
      case_when(
      
      name %in% unique(params.long$name)[1:4] ~ "lam0",
      name %in% unique(params.long$name)[5:8] ~ "alpha2",
      name %in% unique(params.long$name)[9:12] ~ "alpha1"
      
    ),
    
    levels = c("lam0", "alpha2", "alpha1")
    
    ),
    
    coef = factor(
      
      case_when(
        
        name %in% unique(params.long$name)[c(1, 5, 9)] ~ "intercept",
        name %in% unique(params.long$name)[c(2, 6, 10)] ~ "male",
        name %in% unique(params.long$name)[c(3, 7, 11)] ~ "retention",
        name %in% unique(params.long$name)[c(4, 8, 12)] ~ "piling"
        
      ),
      
      levels = c("intercept", "male", "retention", "piling")
    
  )
  
  )

params.summary.plot <- params.summary %>%
  
  mutate(
    
    param = factor(
      
      case_when(
        
        name %in% unique(params.summary$name)[1:4] ~ "alpha1",
        name %in% unique(params.summary$name)[5:8] ~ "alpha2",
        name %in% unique(params.summary$name)[9:12] ~ "lam0"
        
      ),
      
      levels = c("lam0", "alpha2", "alpha1")
      
    ),
    
    coef = factor(
      
      case_when(
        
        name %in% unique(params.summary$name)[c(1, 5, 9)] ~ "intercept",
        name %in% unique(params.summary$name)[c(2, 6, 10)] ~ "male",
        name %in% unique(params.summary$name)[c(3, 7, 11)] ~ "retention",
        name %in% unique(params.summary$name)[c(4, 8, 12)] ~ "piling"
        
      ),
      
      levels = c("intercept", "male", "retention", "piling")
      
    )
  
  )

# ______________________________________________________________________________
# 3b. Half-eye plots ----

# labeller
greek.labels <- c(
  
  "lam0" = "lambda[0]",
  "alpha2" = "alpha[2]",
  "alpha1" = "alpha[1]"
  
)

# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ param,
             scales = "free_y",
             nrow = 3,
             labeller = labeller(param = as_labeller(greek.labels,
                                                     label_parsed))) +
  
  # vertical line at 0
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  
  # KDE
  geom_density_ridges(data = params.long.plot,
                      aes(x = value,
                          y = coef),
                      color = NA,
                      alpha = 0.5,
                      scale = 0.5,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  geom_errorbar(data = params.summary.plot,
                aes(xmin = lo.2,
                    xmax = up.2,
                    y = coef),
                alpha = 0.40,
                height = 0,
                linewidth = 1.35,
                position = position_nudge(y = -0.1)) +
  
  geom_errorbar(data = params.summary.plot,
                aes(xmin = lo.1,
                    xmax = up.1,
                    y = coef),
                height = 0,
                linewidth = 1.35,
                position = position_nudge(y = -0.1)) +
  
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "gray90",
                                        linetype = "blank")) +
  
  # x axis title
  xlab("Parameter estimate")

# 334 x 491

# ______________________________________________________________________________
# 4. Write table ----
# ______________________________________________________________________________

write.table(params.summary, "clipboard", sep = "\t")

# ______________________________________________________________________________
# 5. Detection function predictions ----
# ______________________________________________________________________________
# 5a. Function ----

# this accepts a matrix of posterior samples x
# and a data.frame to predict on y
  # includes indicators for covariates
  # and distance values

x.lam <- params.3 |> do.call(rbind, args = _) |> as.data.frame()

# a test with all zeroes
test.y <- data.frame(
  
  d = seq(0, 500, length.out = 500),
  sex = 0,
  ret = 0,
  pil = 0,
  prev.cap = 1
  
)

# ______________________________________________________________________________

predict_dfn <- function (x, y) {
  
  # change names
  names(x) <- c("lam0_b0", "lam0_b1", "lam0_b2", "lam0_b3",
                "alpha2_b0", "alpha2_b1", "alpha2_b2", "alpha2_b3",
                "alpha1_b0", "alpha1_b1", "alpha1_b2", "alpha1_b3")
  
  # calculate predictions 
  # alpha2 - linear scale
  if (y$prev.cap[1] == 1) {
    
    alpha2 <- x$alpha2_b0 +
      
      x$alpha1_b1 * y$sex +
      x$alpha1_b2 * y$ret +
      x$alpha1_b3 * y$pil
    
  } else {
    
    alpha2 = 0
    
  }
  
  # lam0 - hazard scale
  lam0 <- exp(
    
    x$lam0_b0 +
      
      x$lam0_b1 * y$sex +
      x$lam0_b2 * y$ret +
      x$lam0_b3 * y$pil +
      
      alpha2
    
  )
  
  # alpha1 - must inverse the logit
  alpha1 <- 1 / 
    
    (1 + exp(
      
      -1 * (
        
        x$alpha1_b0 +
          
          x$alpha1_b1 * y$sex +
          x$alpha1_b2 * y$ret +
          x$alpha1_b3 * y$pil
        
      )
      
      )
     
     ) 
  
  # calculate eta (un-normalized probability)
  # posterior for each distance value - d.alph
  eta <- matrix(data = NA, nrow = length(alpha1), ncol = length(y$d))
  
  for (i in 1:nrow(eta)) {
    
    for (j in 1:ncol(eta)) {
      
      eta[i, j] <- lam0[i] * exp(-alpha1[i] * y$d[j]) 
      
    }
    
  }
  
  # summarize for each distance
  y$pred.med <- apply(eta, 2, median)
  y$pred.l90 <- apply(eta, 2, quantile, prob = 0.05)
  y$pred.u90 <- apply(eta, 2, quantile, prob = 0.95)
  
  # return
  return(y)
  
}

# test function
predict_dfn(x.lam, test.y)

# ______________________________________________________________________________
# 5b. Predict ----
# ______________________________________________________________________________

# data.frames to predict on
all.df 


all.preds <- rbind(
  
  test.y <- data.frame(
    
    d = seq(0, 500, length.out = 100),
    sex = 0,
    ret = 0,
    pil = 0,
    prev.cap = 1
    
  )
  
)