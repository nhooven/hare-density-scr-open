# PROJECT: Open multi-session SCR
# SCRIPT: 07b - Visualization (open parameters)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Apr 2026
# COMPLETED: 28 Apr 2026
# LAST MODIFIED: 11 Nov 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(MCMCvis)
library(coda)
library(bayestestR)
library(ggridges)        # ridgeline plot
library(mefa4)

# ______________________________________________________________________________
# 2. Read in samples for phi ----
# ______________________________________________________________________________

params.2 <- readRDS("final_samples/samples_clean_2.rds")

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
# 3a. Process for phi ----
# ______________________________________________________________________________

# subset data
params.phi <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = 1:28) |>
  
  long_draws()

# extract summaries
phi.summary <- extract_med_ci(params.phi)

# prepare for plotting
# factors
params.phi$name <- factor(params.phi$name)
phi.summary$name <- factor(phi.summary$name)

# ______________________________________________________________________________
# 3b. Half-eye plots - phi ----

# which coefficients/groups?
unique(params.phi$name)

phi.coefs.full <- params.phi %>%
  
  # subset to only coefficients
  filter(name %in% unique(params.phi$name)[20:28]) %>%
  
  # grouping for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.phi$name)[23:26] ~ "trt",
        
        # other
        TRUE ~ "default"
        
      )
      
    )
    
  )

phi.coefs.summ <- phi.summary %>% 
  
  filter(name %in% unique(params.phi$name)[20:28]) %>%
  
  # grouping for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.phi$name)[23:26] ~ "trt",
        
        # other
        TRUE ~ "default"
        
      )
      
    )
    
  )

# labels (must be reversed)
phi.labels <- rev(c("Male", expression(tau[2]), expression(tau[3]), 
                    expression("RET x " * tau[2]),
                    expression("RET x " * tau[3]),
                    expression("PIL x " * tau[2]),
                    expression("PIL x " * tau[3]),
                    "DM", "OPEN"))
  
# ______________________________________________________________________________

ggplot() +
  
  theme_classic() +
  
  # vertical line at 0
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  
  # KDE
  geom_density_ridges(data = phi.coefs.full,
                      aes(x = value,
                          y = name,
                          fill = coef.group),
                      color = NA,
                      alpha = 0.5,
                      scale = 0.5,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  geom_errorbar(data = phi.coefs.summ,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = coef.group),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.35,
                 position = position_nudge(y = -0.1)) +
  
  geom_errorbar(data = phi.coefs.summ,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = coef.group),
                 height = 0,
                 linewidth = 1.35,
                 position = position_nudge(y = -0.1)) +
  
  # x axis title
  xlab("Parameter estimate") +
  
  # y-axis scale
  scale_y_discrete(labels = phi.labels,
                   limits = rev) +
  
  # colors
  scale_color_manual(values = c("gray45", "#CC0033")) +
  scale_fill_manual(values = c("gray45", "#CC0033")) +
  
  # coordinates
  coord_cartesian(xlim = c(-2.5, 2.5),
                  ylim = c(1.3, 9)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        legend.position = "none")

# 311 x 322

# move this to the supplement- only include predictions in the main text

# ______________________________________________________________________________
# 4. Predictions ----
# ______________________________________________________________________________
# 4a. Phi ----

# inverse logit function
inv_logit <- function (x) {
  
  y <- 1 / (1 + exp(-x))
  
  return(y)
  
}

# ______________________________________________________________________________

# sex-year-trt predictions
x.phi <- params.2 |> do.call(rbind, args = _) |> as.data.frame()

# CONTROLS
phi.ctrl.1.f <- inv_logit(x.phi$phi_b0)
phi.ctrl.2.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[2]`)
phi.ctrl.3.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[3]`)

phi.ctrl.1.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]`)
phi.ctrl.2.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[2]`)
phi.ctrl.3.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[3]`)

# RET
phi.ret.2.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[2]` + x.phi$`phi_b[4]`)
phi.ret.3.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[3]` + x.phi$`phi_b[5]`)

phi.ret.2.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[2]` + x.phi$`phi_b[4]`)
phi.ret.3.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[3]` + x.phi$`phi_b[5]`)

# PIL
phi.pil.2.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[2]` + x.phi$`phi_b[6]`)
phi.pil.3.f <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[3]` + x.phi$`phi_b[7]`)

phi.pil.2.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[2]` + x.phi$`phi_b[6]`)
phi.pil.3.m <- inv_logit(x.phi$phi_b0 + x.phi$`phi_b[1]` + x.phi$`phi_b[3]` + x.phi$`phi_b[7]`)

# summaries
phi.summ <- data.frame(
  
  # medians
  med = c(
    
    # CONTROL
    median(phi.ctrl.1.f), 
    median(phi.ctrl.2.f),
    median(phi.ctrl.3.f),
    median(phi.ctrl.1.m), 
    median(phi.ctrl.2.m),
    median(phi.ctrl.3.m),
    
    # RET
    median(phi.ret.2.f),
    median(phi.ret.3.f), 
    median(phi.ret.2.m),
    median(phi.ret.3.m),
    
    # PIL 
    median(phi.pil.2.f),
    median(phi.pil.3.f), 
    median(phi.pil.2.m),
    median(phi.pil.3.m)
    
  ),
  
  # l90
  l90 = c(
    
    # CONTROL
    as.numeric(hdi(phi.ctrl.1.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.ctrl.2.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.ctrl.3.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.ctrl.1.m, ci = 0.90))[2], 
    as.numeric(hdi(phi.ctrl.2.m, ci = 0.90))[2], 
    as.numeric(hdi(phi.ctrl.3.m, ci = 0.90))[2], 
    
    # RET
    as.numeric(hdi(phi.ret.2.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.ret.3.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.ret.2.m, ci = 0.90))[2], 
    as.numeric(hdi(phi.ret.3.m, ci = 0.90))[2],  
    
    # PIL
    as.numeric(hdi(phi.pil.2.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.pil.3.f, ci = 0.90))[2], 
    as.numeric(hdi(phi.pil.2.m, ci = 0.90))[2], 
    as.numeric(hdi(phi.pil.3.m, ci = 0.90))[2]  
    
  ),
  
  # u90
  u90 = c(
    
    # CONTROL
    as.numeric(hdi(phi.ctrl.1.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.ctrl.2.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.ctrl.3.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.ctrl.1.m, ci = 0.90))[3], 
    as.numeric(hdi(phi.ctrl.2.m, ci = 0.90))[3], 
    as.numeric(hdi(phi.ctrl.3.m, ci = 0.90))[3], 
    
    # RET
    as.numeric(hdi(phi.ret.2.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.ret.3.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.ret.2.m, ci = 0.90))[3], 
    as.numeric(hdi(phi.ret.3.m, ci = 0.90))[3],  
    
    # PIL
    as.numeric(hdi(phi.pil.2.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.pil.3.f, ci = 0.90))[3], 
    as.numeric(hdi(phi.pil.2.m, ci = 0.90))[3], 
    as.numeric(hdi(phi.pil.3.m, ci = 0.90))[3]
    
  ),
  
  # identifiers
  period = c(1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3),
  trt = factor(c(rep("unthinned", 6), rep("retention", 4), rep("piling", 4)),
               levels = c("unthinned", "retention", "piling")),
  sex = factor(c("F", "F", "F", 
                 "M", "M", "M", 
                 "F", "F", "M", "M", 
                 "F", "F", "M", "M")
               
  )
  
)

# prediction plot
ggplot(data = phi.summ) +
  
  theme_bw() +
  
  facet_grid(~ trt) +
  
  # connecting lines
  geom_line(aes(x = period,
                y = med,
                color = sex,
                linetype = sex),
            linewidth = 0.8) +
  
  # CIs
  geom_errorbar(aes(x = period,
                    y = med,
                    ymin = l90,
                    ymax = u90,
                    color = sex),
                alpha = 0.45,
                width = 0,
                linewidth = 1.5) +
  
  # points
  geom_point(aes(x = period,
                 y = med,
                 shape = sex,
                 color = sex,
                 size = sex),
             fill = "white",
             stroke = 0.8) +
  
  # scales
  scale_shape_manual(values = c(21, 23)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = c(0.85, 0.85),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "gray90",
                                        linetype = "blank"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        plot.margin = margin(0.1, 0.1, 0.02, 0.1, unit = "cm")) +
  
  # greek letters
  scale_x_continuous(
    
    breaks = 1:3,
    
    labels = c("1" = expression(tau[1]),
               "2" = expression(tau[2]),
               "3" = expression(tau[3])),
    
    limits = c(0.75, 3.25)
    
    ) +
  
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25, 0.35, 0.45)) +
  
  # color and size
  scale_color_manual(values = c("#FF3300", "gray25")) +
  scale_size_manual(values = c(2, 1.8)) +
  
  # labels
  ylab("Persistence probability") +
  xlab("Transition period")
  
# 408 x 251
# we'll probably not use this one standalone

# ______________________________________________________________________________
# 4b. Rho ----
# ______________________________________________________________________________

# which ones?
real.rho <- c(30, 31, 33, 34, 36:49, 51:55, 57, 58, 60, 61, 63, 64)

# subset data
params.rho <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = real.rho) |>
  
  long_draws()

# extract summaries
rho.summary <- extract_med_ci(params.rho)

# add site and tau
rho.summary <- rho.summary |>
  
  mutate(
    
    site = c(1, 1, 10, 10, 11, 11, 12, 12, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9, 9),
    tau = c(2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3)
    
  ) |>
  
  # trt and cluster
  mutate(
    
    trt = case_when(
      
      site %in% c(3, 6, 9, 12) ~ "control",
      site %in% c(1, 5, 8, 10) ~ "retention",
      site %in% c(2, 4, 7, 11) ~ "piling"
      
    ),
    
    cluster = case_when(
      
      site %in% 1:3 ~ "cluster 1",
      site %in% 4:6 ~ "cluster 2",
      site %in% 7:9 ~ "cluster 3",
      site %in% 10:12 ~ "cluster 4"
      
    )
    
  ) |>
  
  mutate(trt = factor(trt, levels = c("control", "retention", "piling")))

# trt x tau summaries
trt.tau.sum <- params.rho |>
  
  # join in trt and tau
  left_join(
    
    rho.summary |> dplyr::select(name, trt, tau)
    
  ) |>
  
  # summarize
  group_by(trt, tau) |>
  
  summarize(med = median(value),
            lo.2 = as.numeric(hdi(value, ci = 0.90))[2],
            up.2 = as.numeric(hdi(value, ci = 0.90))[3]) |>
  
  ungroup() |>
  
  mutate(cluster = "all",
         name = paste0(trt, tau))

# join in
rho.summary.1 <- bind_rows(rho.summary, trt.tau.sum) |>
  
  # factor levels
  mutate(cluster = factor(cluster, levels = c("all", 
                                              "cluster 1", 
                                              "cluster 2", 
                                              "cluster 3", 
                                              "cluster 4", 
                                              "cluster 5"))) |>
  
  # for colors
  mutate(group = ifelse(cluster == "all",
                        "all",
                        "clusters"))

# plot
ggplot(data = rho.summary.1) +
  
  theme_bw() +
  
  facet_grid(trt ~ cluster) +
  
  # individual sites
  geom_errorbar(data = rho.summary.1,
                aes(x = tau,
                    y = med, 
                    ymin = lo.2,
                    ymax = up.2,
                    color = group),
                width = 0,
                alpha = 0.35,
                linewidth = 1.5) +
  
  # connecting line
  geom_line(aes(x = tau,
                y = med,
                color = group),
            linewidth = 0.8) +
  
  geom_point(data = rho.summary.1,
             aes(x = tau,
                 y = med,
                 color = group),
             shape = 21,
             size = 1.25,
             stroke = 0.8,
             fill = "white") +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "gray90",
                                        linetype = "blank"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        plot.margin = margin(0.1, 0.1, 0.02, 0.1, unit = "cm")) +
  
  # greek letters
  scale_x_continuous(
    
    breaks = 1:3,
    
    labels = c("1" = expression(tau[1]),
               "2" = expression(tau[2]),
               "3" = expression(tau[3])),
    
    limits = c(0.75, 3.25)
    
  ) +
  
  scale_color_manual(values = c("gray45", "black")) +
  
  scale_y_continuous(breaks = c(0.5, 1.5, 2.5)) +
  
  # labels
  ylab("Realized per capita recruitment") +
  xlab("Transition period")

# 598 x 385
# supplement? or not?

# ______________________________________________________________________________
# 5. Random intercepts ----
# ______________________________________________________________________________
# 5a. Phi ----
# ______________________________________________________________________________

params.phi.clust <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = 1:19) |>
  
  as.data.frame() %>%
  
  mutate(
    
    c1 = inv_logit(`phi_c[1]`),
    c2 = inv_logit(`phi_c[2]`),
    c3 = inv_logit(`phi_c[3]`),
    c4 = inv_logit(`phi_c[4]`),
    u1 = inv_logit(`phi_c[1]` + `phi_u_sd` * `phi_u_z[1]`),
    u2 = inv_logit(`phi_c[1]` + `phi_u_sd` * `phi_u_z[2]`),
    u3 = inv_logit(`phi_c[1]` + `phi_u_sd` * `phi_u_z[3]`),
    u4 = inv_logit(`phi_c[2]` + `phi_u_sd` * `phi_u_z[4]`),
    u5 = inv_logit(`phi_c[2]` + `phi_u_sd` * `phi_u_z[5]`),
    u6 = inv_logit(`phi_c[2]` + `phi_u_sd` * `phi_u_z[6]`),
    u7 = inv_logit(`phi_c[3]` + `phi_u_sd` * `phi_u_z[7]`),
    u8 = inv_logit(`phi_c[3]` + `phi_u_sd` * `phi_u_z[8]`),
    u9 = inv_logit(`phi_c[3]` + `phi_u_sd` * `phi_u_z[9]`),
    u10 = inv_logit(`phi_c[4]` + `phi_u_sd` * `phi_u_z[10]`),
    u11 = inv_logit(`phi_c[4]` + `phi_u_sd` * `phi_u_z[11]`),
    u12 = inv_logit(`phi_c[4]` + `phi_u_sd` * `phi_u_z[12]`)
    
  ) %>%
  
  # keep intercept columns
  dplyr::select(c1:u12) %>%
  
  # summarize
  pivot_longer(cols = c1:u12) %>%
  
  group_by(name) %>%
  
  summarize(
    
    med = median(value),
    sd = sd(value),
    lo.1 = as.numeric(hdi(value, ci = 0.5)[2]),
    up.1 = as.numeric(hdi(value, ci = 0.5)[3]),
    lo.2 = as.numeric(hdi(value, ci = 0.9)[2]),
    up.2 = as.numeric(hdi(value, ci = 0.9)[3])
    
  ) |>
  
  # add ID columns for plotting
  mutate(cluster = case_when(
    
    name %in% c("c1", "u1", "u2", "u3") ~ "cluster 1",
    name %in% c("c2", "u4", "u5", "u6") ~ "cluster 2",
    name %in% c("c3", "u7", "u8", "u9") ~ "cluster 3",
    name %in% c("c4", "u10", "u11", "u12") ~ "cluster 4"),
    
    level = case_when(
      
      name %in% c("c1", "c2", "c3", "c4") ~ "cluster",
      TRUE ~ "site"
      
    )
  )

# random intercept plot
ggplot(data = params.phi.clust) +
  
  theme_bw() +
  
  facet_wrap(~ cluster, 
             scales = "free_x",
             nrow = 1) +
  
  # 90%
  geom_errorbar(aes(x = name,
                    ymin = lo.2,
                    ymax = up.2,
                    color = level),
                alpha = 0.4,
                width = 0,
                linewidth = 1.75) +
  
  # 50% CI
  geom_errorbar(aes(x = name,
                    ymin = lo.1,
                    ymax = up.1,
                    color = level),
                width = 0,
                linewidth = 1.75) +
  
  # median
  geom_point(aes(x = name,
                 y = med,
                 color = level),
             shape = 21,
             fill = "white",
             size = 2,
             stroke = 0.8) +
  
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  
  scale_color_manual(values = c("black", "gray55")) +

  ylab(expression(phi))

# 591 x 192

# ______________________________________________________________________________
# 6. Write tables ----
# ______________________________________________________________________________

write.table(phi.summary, "clipboard", sep = "\t")
write.table(rho.summary, "clipboard", sep = "\t")
write.table(rho.summary.1, "clipboard", sep = "\t")

# ______________________________________________________________________________
# 7. Derived delta D ----
# ______________________________________________________________________________
# 7a. Process ----
# ______________________________________________________________________________

# density
params.1 <- readRDS("final_samples/samples_clean_1.rds")

# just the N values
d.df <- params.1 |> do.call(rbind, args = _) |>
  
  as.data.frame() |>
  
  dplyr::select(1:48)

# S areas
S.areas <- readRDS("for_model/S_area.rds")

# constants
constant.list <- readRDS("for_model/constants.rds")

# site year lookup (these are in order)
site.year <- expand.grid(site = 1:12, year = 1:4)

# loop through each site-year and calculate D
all.D <- matrix(data = NA,
                nrow = nrow(d.df),
                ncol = ncol(d.df))

for (j in 1:48) {
  
  # focal site.year
  focal.col = site.year[j, ]
  
  # summarize (only if we trapped that year)
  if (focal.col$year == 1 & constant.list$first.year[focal.col$site] == 2) {
    
    focal.D <- vector(length = nrow(d.df))
    
  } else {
    
    # posterior density
    all.D[ , j] <- d.df[ , j] / S.areas$area.rect[focal.col$site]
    
  }
  
} # j

# ______________________________________________________________________________
# 7c. Calculate deltaD and percent deltaD - site-tau ----
# ______________________________________________________________________________

# site lookup
site.lookup <- read.csv("for_model/site_lookup.csv") |>
  
  group_by(SiteName) |>
  slice(1)  # keep only 1

# function
calc_deltaD <- function (.siteName,
                         .tau) {
  
  # translate .site to integer
  .siteID <- site.lookup$siteID[site.lookup$SiteName == .siteName]
  
  # if .tau == 1 and first.year == 2, throw an error
  if (.tau == 1 & constant.list$first.year[.siteID] == 2) { 
    
    print("Site not trapped during t1!") 
    
    # if good, proceed
  } else {
      
    # correct years
    if (.tau == 1) { t1 <- 1 ; t2 <- 2 }
    if (.tau == 2) { t1 <- 2 ; t2 <- 3 }
    if (.tau == 3) { t1 <- 3 ; t2 <- 4 }
    
    # pick correct columns
    col.t1 <- all.D[ , which(site.year$site == .siteID & site.year$year == t1)]
    col.t2 <- all.D[ , which(site.year$site == .siteID & site.year$year == t2)]
    
    # change metrics
    deltaD <- col.t2 - col.t1
    perc.deltaD <- (deltaD / col.t1) * 100
    
    # summarize and bind in
    suppressMessages(
    focal.deltaD <- data.frame(
      
      # IDs
      site = .siteName,
      tau = .tau,
      
      deltaD.med = median(deltaD),
      deltaD.l90 = as.numeric(hdi(deltaD, ci = 0.9)[2]),
      deltaD.u90 = as.numeric(hdi(deltaD, ci = 0.9)[3]),
      
      perc.deltaD.med = median(perc.deltaD),
      perc.deltaD.l90 = as.numeric(hdi(perc.deltaD, ci = 0.9)[2]),
      perc.deltaD.u90 = as.numeric(hdi(perc.deltaD, ci = 0.9)[3])
      
    )
    )
    
    return(focal.deltaD)
    
    }
  
} # fn

# bind together
all.deltaD <- bind_rows(
  
  calc_deltaD("2A", 1),
  calc_deltaD("2B", 1),
  calc_deltaD("2C", 1),
  calc_deltaD("3A", 1),
  calc_deltaD("3C", 1),
  
  calc_deltaD("1A", 2),
  calc_deltaD("1B", 2),
  calc_deltaD("1C", 2),
  calc_deltaD("2A", 2),
  calc_deltaD("2B", 2),
  calc_deltaD("2C", 2),
  calc_deltaD("3A", 2),
  calc_deltaD("3B", 2),
  calc_deltaD("3C", 2),
  calc_deltaD("4A", 2),
  calc_deltaD("4B", 2),
  calc_deltaD("4C", 2),
  
  calc_deltaD("1A", 3),
  calc_deltaD("1B", 3),
  calc_deltaD("1C", 3),
  calc_deltaD("2A", 3),
  calc_deltaD("2B", 3),
  calc_deltaD("2C", 3),
  calc_deltaD("3A", 3),
  calc_deltaD("3B", 3),
  calc_deltaD("3C", 3),
  calc_deltaD("4A", 3),
  calc_deltaD("4B", 3),
  calc_deltaD("4C", 3)
  
)

# ______________________________________________________________________________
# 7d. Calculate deltaD and percent deltaD - treatment-tau ----
# ______________________________________________________________________________

# function
calc_deltaD_trt <- function (.trt,
                             .tau) {
  
  # which sites?
  if (.trt == "untreated" & .tau == 1) { focal.sites <- which(constant.list$first.year == 1) }
  if (.trt == "untreated" & .tau != 1) { focal.sites <- c(3, 6, 9, 12) }
  if (.trt == "retention" & .tau != 1) { focal.sites <- c(1, 5, 8, 10) }
  if (.trt == "piling" & .tau != 1) { focal.sites <- c(2, 4, 7, 11) }
  if (.trt != "untreated" & .tau == 1) { stop("Not valid!") }
  
  # correct years
  if (.tau == 1) { t1 <- 1 ; t2 <- 2 }
  if (.tau == 2) { t1 <- 2 ; t2 <- 3 }
  if (.tau == 3) { t1 <- 3 ; t2 <- 4 }
  
  # loop through focal sites, calculating deltaD and perc.delta D
  focal.deltaD <- data.frame()
  
  for (i in focal.sites) {
    
    # pick correct columns
    col.t1 <- all.D[ , which(site.year$site == i & site.year$year == t1)]
    col.t2 <- all.D[ , which(site.year$site == i & site.year$year == t2)]
    
    # change metrics
    deltaD <- col.t2 - col.t1
    perc.deltaD <- (deltaD / col.t1) * 100
    
    # bind into "hyper"posterior (?)
    focal.deltaD <- rbind(focal.deltaD, cbind(deltaD, perc.deltaD))
    
  }
  
   # summarize
    suppressMessages(
      focal.deltaD.summary <- data.frame(
        
        # IDs
        trt = .trt,
        tau = .tau,
        
        deltaD.med = median(focal.deltaD$deltaD),
        deltaD.l90 = as.numeric(hdi(focal.deltaD$deltaD, ci = 0.9)[2]),
        deltaD.u90 = as.numeric(hdi(focal.deltaD$deltaD, ci = 0.9)[3]),
        
        perc.deltaD.med = median(focal.deltaD$perc.deltaD),
        perc.deltaD.l90 = as.numeric(hdi(focal.deltaD$perc.deltaD, ci = 0.9)[2]),
        perc.deltaD.u90 = as.numeric(hdi(focal.deltaD$perc.deltaD, ci = 0.9)[3])
        
      )
    )
    
    return(focal.deltaD.summary)
  
} # fn

# bind together
all.deltaD.trt <- bind_rows(
  
  calc_deltaD_trt("untreated", 1),
  
  calc_deltaD_trt("untreated", 2),
  calc_deltaD_trt("retention", 2),
  calc_deltaD_trt("piling", 2),
  
  calc_deltaD_trt("untreated", 3),
  calc_deltaD_trt("retention", 3),
  calc_deltaD_trt("piling", 3)
  
) |>
  
  # factor levels
  mutate(trt = factor(trt, levels = c("untreated", "retention", "piling")))

# ______________________________________________________________________________
# 8. Final 3-panel plot ----

# theme
theme_3panel <- function () {
  
  theme_bw() +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"))
  
}

# ______________________________________________________________________________
# 8a. ROW 1 - Percent delta D by treatment ----
# ______________________________________________________________________________

plot.deltaD <- ggplot() +
  
  theme_3panel() +
  
  facet_wrap(~ trt) +
  
  # hline at 0
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  # connecting lines
  geom_line(data = all.deltaD.trt,
            aes(x = tau,
                y = perc.deltaD.med),
            color = "gray45",
            linewidth = 0.8) +
  
  # error bars
  geom_errorbar(data = all.deltaD.trt,
                aes(x = tau,
                    y = perc.deltaD.med,
                    ymin = perc.deltaD.l90,
                    ymax = perc.deltaD.u90),
                color = "gray45",
                alpha = 0.45,
                width = 0,
                linewidth = 1.5) +
  
  geom_point(data = all.deltaD.trt,
             aes(x = tau,
                 y = perc.deltaD.med),
             color = "gray45",
             size = 2,
             shape = 21,
             fill = "white",
             stroke = 0.8) +
    
  # additional theme arguments
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = margin(t = 2, b = 0, r = 2, l = 5)) +
  
  scale_x_continuous(breaks = 1:3, limits = c(0.75, 3.25)) +
  scale_y_continuous(breaks = c(-50, 0, 50, 100, 150)) +
  
  # titles
  ylab(expression("% "*Delta*"D"))

# ______________________________________________________________________________
# 8b. ROW 2 - Predicted phi by treatment ----
# ______________________________________________________________________________

plot.phi <- ggplot(data = phi.summ) +
  
  theme_3panel() +
  
  facet_wrap(~ trt) +
  
  # connecting lines
  geom_line(aes(x = period,
                y = med,
                color = sex,
                linetype = sex),
            linewidth = 0.8) +
  
  # error bars
  geom_errorbar(aes(x = period,
                    y = med,
                    ymin = l90,
                    ymax = u90,
                    color = sex),
                alpha = 0.45,
                width = 0,
                linewidth = 1.5) +
  
  geom_point(aes(x = period,
                 y = med,
                 shape = sex,
                 color = sex,
                 size = sex),
             fill = "white",
             stroke = 0.8) +
  
  # additional theme arguments
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 11),
        
        legend.position = c(0.75, 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.key.spacing.x = unit(0.1, "cm"),
        legend.key.spacing.y = unit(0.001, "cm"),
        legend.key.size = unit(0.5, "cm"),
        
        plot.margin = margin(t = 0, b = 0, r = 2, l = 2)) +
  
  # scales
  scale_shape_manual(values = c(21, 23)) +
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25, 0.35, 0.45)) +
  scale_x_continuous(breaks = 1:3, limits = c(0.75, 3.25)) +
  scale_color_manual(values = c("#FF3300", "gray25")) +
  scale_size_manual(values = c(2, 1.6)) +
  
  # titles
  ylab(expression(phi))

plot.phi

# ______________________________________________________________________________
# 8c. ROW 3 - Derived rho by treatment ----

rho.summ <- params.rho |>
  
  # join in trt and tau
  left_join(
    
    rho.summary |> dplyr::select(name, trt, tau)
    
  ) |>
      
    # correct trt to mirror the other 2 panels
    mutate(trt.1 = case_when(
        
      trt == "control" ~ "untreated",
      tau == 1 ~ "untreated",
      tau != 1 & trt == "retention" ~ "retention",
      tau != 1 & trt == "piling" ~ "piling"
        
      )
      
    ) |>
  
  # summarize
  group_by(trt.1, tau) |>
  
  summarize(med = median(value),
            l90 = as.numeric(hdi(value, ci = 0.90))[2],
            u90 = as.numeric(hdi(value, ci = 0.90))[3]) |>
  
  ungroup() |>
  
  # factor levels
  mutate(trt = factor(trt.1, levels = c("untreated", "retention", "piling")))

# ______________________________________________________________________________

plot.rho <- ggplot(data = rho.summ) +
  
  theme_3panel() +
  
  facet_wrap(~ trt) +
  
  # connecting lines
  geom_line(aes(x = tau,
                y = med),
            color = "gray45",
            linewidth = 0.8) +
  
  # error bars
  geom_errorbar(aes(x = tau,
                    y = med,
                    ymin = l90,
                    ymax = u90),
                color = "gray45",
                alpha = 0.45,
                width = 0,
                linewidth = 1.5) +
  
  geom_point(aes(x = tau,
                 y = med),
             color = "gray45",
             fill = "white",
             size = 2,
             shape = 21,
             stroke = 0.8) +
  
  # additional theme arguments
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0, b = 2, r = 2, l = 8.5)) +
  
  # scales
  scale_x_continuous(
    
    breaks = 1:3,
    
    labels = c("1" = expression(tau[1]),
               "2" = expression(tau[2]),
               "3" = expression(tau[3])),
    
    limits = c(0.75, 3.25)
    
  ) +
  
  # titles
  xlab("Transition period") +
  ylab(expression(rho))

# ______________________________________________________________________________
# 8d. Plot all ----
# ______________________________________________________________________________

cowplot::plot_grid(plot.deltaD, plot.phi, plot.rho,
          nrow = 3,
          rel_heights = c(1.0, 0.8, 1.05))

# 401 x 396
ggsave("final_plots/Figure_tau.pdf",
       device = "pdf",
       width = 401,
       height = 396,
       scale = 3.5,
       units = "px")

# ______________________________________________________________________________
# 9. Write to tables ----
# ______________________________________________________________________________

write.table(all.deltaD.trt, "clipboard", sep = "\t")
write.table(phi.summ, "clipboard", sep = "\t")
write.table(rho.summ, "clipboard", sep = "\t")

# ______________________________________________________________________________
# 11. Site-by-site scaled open parameters ----

# instead of making separate graphs, let's do it all in one facetted plot
# similar to the site-by-site density

# ______________________________________________________________________________
# 11a. Process delta D and rho ----
# ______________________________________________________________________________

# deltaD and rho are easy - we have them already, just need to scale
scaled.deltaD <- all.deltaD |>
  
  dplyr::select(site, tau, perc.deltaD.med, perc.deltaD.l90, perc.deltaD.u90) |>
  
  mutate(med = (perc.deltaD.med - mean(all.deltaD$perc.deltaD.med)) / sd(all.deltaD$perc.deltaD.med),
         l90 = (perc.deltaD.l90 - mean(all.deltaD$perc.deltaD.med)) / sd(all.deltaD$perc.deltaD.med),
         u90 = (perc.deltaD.u90 - mean(all.deltaD$perc.deltaD.med)) / sd(all.deltaD$perc.deltaD.med)) |>
  
  dplyr::select(site, tau, med, l90, u90) |>
  
  # variable
  mutate(var = "perc.deltaD") |>
  
  # overwrite site
  left_join(
    
    site.lookup |> dplyr::select(SiteName, siteID) |>
      
      rename(site = SiteName)
    
  ) |>
  
  # columns in correct order, drop site, rename
  dplyr::select(siteID, tau:var) |>
  rename(site = siteID)

scaled.rho <- rho.summary |>
  
  dplyr::select(site, tau, med, lo.2, up.2) |>
  
  mutate(med = (med - mean(rho.summary$med)) / sd(rho.summary$med),
         l90 = (lo.2 - mean(rho.summary$med)) / sd(rho.summary$med),
         u90 = (up.2 - mean(rho.summary$med)) / sd(rho.summary$med)) |>
  
  dplyr::select(-c(lo.2, up.2)) |>
  
  # variable
  mutate(var = "rho")

# bind together
scaled.D.rho <- bind_rows(scaled.deltaD, scaled.rho)

# ______________________________________________________________________________
# 11b. Process phi (conditional predictions) ----
# ______________________________________________________________________________

# predict df
phi.pred.df <- data.frame(
  
  site = c(4, 5, 6, 7, 9,
           1:12,
           1:12),
  cluster = c(2, 2, 2, 3, 3,
              rep(1:4, each = 3),
              rep(1:4, each = 3)),
  tau = c(rep(1, 5), rep(2, 12), rep(3, 12)),
  ret = c(0, 1, 0, 0, 0,
          1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
          1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0),
  pil = c(1, 0, 0, 1, 0,
          0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
          0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0),
  post1 = c(0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  post2 = c(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  
) |>
  
  # join in DM and OPEN
  left_join(
    
    data.frame(
      
      site = constant.list$site,
      dm = constant.list$op.dm,
      open = constant.list$op.o
      
    ) |>
      
      group_by(site) |>
      slice(1)

  )

# split by row
phi.pred.df$rowname <- rownames(phi.pred.df)
phi.pred.split <- split(phi.pred.df, phi.pred.df$rowname)

# function to calculate site-tau predictions and summaries
# split x.phi so I don't have to mess with names
x.phi.means.sds <- x.phi[ , 1:3]
x.phi.c <- x.phi[, 4:7]
x.phi.u <- x.phi[ , 8:19]
x.phi.b <- x.phi[ , 20:28]

pred_phi_site <- function (.row) {
  
  row.pred <- inv_logit(
    
    x.phi.c[ , .row$cluster] + x.phi.means.sds$phi_u_sd * x.phi.u[ , .row$s] +
      
      x.phi.b[, 2] * .row$post1 +
      x.phi.b[, 3] * .row$post2 +
      x.phi.b[, 4] * .row$ret * .row$post1 +
      x.phi.b[, 5] * .row$ret * .row$post2 +
      x.phi.b[, 6] * .row$pil * .row$post1 +
      x.phi.b[, 7] * .row$pil * .row$post2 +
      x.phi.b[, 8] * .row$dm +
      x.phi.b[, 9] * .row$open
    
  )
  
  # summary df
  row.pred.summ <- data.frame(
    
    site = .row$site,
    tau = .row$tau,
    med = median(row.pred),
    l90 = as.numeric(hdi(row.pred, ci = 0.90))[2],
    u90 = as.numeric(hdi(row.pred, ci = 0.90))[3]
    
  )
  
}

# apply function
phi.pred.site <- do.call(rbind, lapply(phi.pred.split, pred_phi_site))

# scale and add var name
scaled.phi <- phi.pred.site |>
  
  mutate(med = (med - mean(phi.pred.site$med)) / sd(phi.pred.site$med),
         l90 = (l90 - mean(phi.pred.site$med)) / sd(phi.pred.site$med),
         u90 = (u90 - mean(phi.pred.site$med)) / sd(phi.pred.site$med)) |>
  
  # variable
  mutate(var = "phi")

# ______________________________________________________________________________
# 11c. Bind everything and add in cluster and treatment ----
# ______________________________________________________________________________

scaled.all <- bind_rows(
  
  scaled.D.rho, scaled.phi
  
) |>
  
  mutate(
    
    cluster = case_when(
      
      site %in% 1:3 ~ "cluster 1",
      site %in% 4:6 ~ "cluster 2",
      site %in% 7:9 ~ "cluster 3",
      site %in% 10:12 ~ "cluster 4"
      
    ),
    
    trt = case_when(
      
      site %in% c(3, 6, 9, 12) ~ "control",
      site %in% c(1, 5, 8, 10) ~ "retention",
      site %in% c(2, 4, 7, 11) ~ "piling"
      
    )
    
  ) |>
  
  mutate(trt = factor(trt, levels = c("control", "retention", "piling")))
  
# ______________________________________________________________________________
# 11d. Plot ----
# ______________________________________________________________________________

ggplot(data = scaled.all) +
  
  theme_bw() +
  
  facet_grid(trt ~ cluster) +
  
  geom_errorbar(aes(x = tau,
                    y = med,
                    ymin = l90,
                    ymax = u90,
                    color = var),
                width = 0,
                linewidth = 1,
                alpha = 0.45) +
  
  geom_line(aes(x = tau,
                 y = med,
                 color = var)) +
  
  geom_point(aes(x = tau,
                 y = med,
                 color = var,
                 shape = var),
             fill = "white",
             stroke = 0.8) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0.65, "cm")) +
  
  # scales
  scale_x_continuous(
    
    breaks = 1:3,
    
    labels = c("1" = expression(tau[1]),
               "2" = expression(tau[2]),
               "3" = expression(tau[3])),
    
    limits = c(0.75, 3.25)
    
  ) +
  
  scale_color_manual(labels = c(expression("% "*Delta*"D"),
                                expression(phi),
                                expression(rho)),
                     values = c("gray45", "purple", "orange")) +
  
  scale_shape_manual(labels = c(expression("% "*Delta*"D"),
                                expression(phi),
                                expression(rho)),
                     values = c(21, 23, 24)) +
  
  # titles
  xlab("Transition period") +
  ylab("Standardized estimate")

# 506 x 456

# ______________________________________________________________________________
# 12. Scatter plots of tau parameters ----

# which one is more closely related to density change? obviously rho

# ______________________________________________________________________________
# 12a. Process ----
# ______________________________________________________________________________

df.for.scat <- rho.summary |>
  
  dplyr::select(site, tau, med, lo.2, up.2) |>
  rename(med.rho = med,
         l90.rho = lo.2,
         u90.rho = up.2) |>
  
  left_join(
    
    all.deltaD |>
      
      left_join(site.lookup |> dplyr::select(SiteName, siteID) |> rename(site = SiteName)) |>
      dplyr::select(siteID, tau, perc.deltaD.med, perc.deltaD.l90, perc.deltaD.u90) |>
      rename(site = siteID,
             med.D = perc.deltaD.med,
             l90.D = perc.deltaD.l90,
             u90.D = perc.deltaD.u90)
    
  ) |>
  
  left_join(
    
    phi.pred.site |>
      
      rename(med.phi = med,
             l90.phi = l90,
             u90.phi = u90)
    
  )

# ______________________________________________________________________________
# 12b. Plot ----
# ______________________________________________________________________________

scat.plot.phi <- ggplot(data = df.for.scat) +
  
  theme_classic() +
  
  # error bars
  # x
  geom_errorbar(aes(x = med.phi,
                    xmin = l90.phi,
                    xmax = u90.phi,
                    y = med.D),
                linewidth = 0.75,
                width = 0,
                alpha = 0.25) +
  
  # y
  geom_errorbar(aes(x = med.phi,
                    y = med.D,
                    ymin = l90.D,
                    ymax = u90.D),
                linewidth = 0.75,
                width = 0,
                alpha = 0.25) +
  
  # median points
  geom_point(aes(x = med.phi,
                 y = med.D),
             shape = 21,
             stroke = 0.8,
             fill = "white",
             size = 1.1) +
  
  # labels
  xlab(expression(phi)) +
  ylab(expression("% "*Delta*"D"))

scat.plot.phi

scat.plot.rho <- ggplot(data = df.for.scat) +
  
  theme_classic() +
  
  # error bars
  # x
  geom_errorbar(aes(x = med.rho,
                    xmin = l90.rho,
                    xmax = u90.rho,
                    y = med.D),
                linewidth = 0.75,
                width = 0,
                alpha = 0.25) +
  
  # y
  geom_errorbar(aes(x = med.rho,
                    y = med.D,
                    ymin = l90.D,
                    ymax = u90.D),
                linewidth = 0.75,
                width = 0,
                alpha = 0.25) +
  
  # median points
  geom_point(aes(x = med.rho,
                 y = med.D),
             shape = 21,
             stroke = 0.8,
             fill = "white",
             size = 1.1) +
  
  # labels
  xlab(expression(rho)) +
  
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

scat.plot.rho

# plot together
cowplot::plot_grid(scat.plot.phi, scat.plot.rho,
                   rel_widths = c(1.0, 0.86))

# 516 x 275

# ______________________________________________________________________________
# 13. Table of rho and deltaD ----
# ______________________________________________________________________________

write.table(df.for.scat, "clipboard", sep = "\t")
