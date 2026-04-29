# PROJECT: Open multi-session SCR
# SCRIPT: 07b - Visualization (open parameters)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Apr 2026
# COMPLETED: 28 Apr 2026
# LAST MODIFIED: 29 Apr 2026
# R VERSION: 4.4.3

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
# 2. Read in samples ----
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
  subset(select = 1:15) |>
  
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
  filter(name %in% unique(params.phi$name)[7:15]) %>%
  
  # grouping for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.phi$name)[10:13] ~ "trt",
        
        # landscape
        name %in% unique(params.phi$name)[14:15] ~ "ls",
        
        # other
        TRUE ~ "default"
        
      )
      
    )
    
  )

phi.coefs.summ <- phi.summary %>% 
  
  filter(name %in% unique(params.phi$name)[7:15]) %>%
  
  # grouping for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.phi$name)[10:13] ~ "trt",
        
        # landscape
        name %in% unique(params.phi$name)[14:15] ~ "ls",
        
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
                    "DM", "O"))
  
# ______________________________________________________________________________

phi.est <- ggplot() +
  
  theme_classic() +
  
  # white background
  theme_bw() +
  
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
  scale_color_manual(values = c("gray45", "olivedrab", "firebrick")) +
  scale_fill_manual(values = c("gray45", "olivedrab", "firebrick")) +
  
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
        strip.background = element_blank(),
        plot.margin = margin(t = 0.2, b = 0.2, l = 0.1, r = 0, unit = "cm"),
        legend.position = "none") +
  
  # annotation
  annotate("text", x = 2, y = 9, label = "a")

# 311 x 322

# ______________________________________________________________________________
# 3c. Process for rho ----
# ______________________________________________________________________________

# subset data
params.rho <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = 16:29) |>
  
  long_draws()

# extract summaries
rho.summary <- extract_med_ci(params.rho)

# prepare for plotting
# factors
params.rho$name <- factor(params.rho$name)

rho.summary$name <- factor(rho.summary$name)

# ______________________________________________________________________________
# 3d. Half-eye plots - rho ----

# prepare for plotting
rho.coefs.full <- params.rho %>% 
  
  # subset to only coefficients
  filter(name %in% unique(params.rho$name)[7:14]) %>%
  
  # group for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.rho$name)[9:12] ~ "trt",
        
        # landscape
        name %in% unique(params.rho$name)[13:14] ~ "ls",
        
        # other
        TRUE ~ "default"
        
      )
      
    )
    
  )

rho.coefs.summ <- rho.summary %>% 
  
  filter(name %in% unique(params.rho$name)[7:14]) %>%
  
  # group for colors
  mutate(
    
    coef.group = factor(
      
      case_when(
        
        # treatment
        name %in% unique(params.rho$name)[9:12] ~ "trt",
        
        # landscape
        name %in% unique(params.rho$name)[13:14] ~ "ls",
        
        # other
        TRUE ~ "default"
        
      )
      
    )
    
  )

# labels (must be reversed)
rho.labels <- rev(c(expression(tau[2]), expression(tau[3]), 
                    expression("RET x " * tau[2]),
                    expression("RET x " * tau[3]),
                    expression("PIL x " * tau[2]),
                    expression("PIL x " * tau[3]),
                    "DM", "O"))

# ______________________________________________________________________________

rho.est <- ggplot() +
  
  theme_classic() +
  
  # white background
  theme_bw() +
  
  # vertical line at 0
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  
  # KDE
  geom_density_ridges(data = rho.coefs.full,
                      aes(x = value,
                          y = name,
                          fill = coef.group),
                      color = NA,
                      alpha = 0.5,
                      scale = 0.5,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  geom_errorbar(data = rho.coefs.summ,
                aes(xmin = lo.2,
                    xmax = up.2,
                    y = name,
                    color = coef.group),
                alpha = 0.40,
                height = 0,
                linewidth = 1.35,
                position = position_nudge(y = -0.1)) +
  
  geom_errorbar(data = rho.coefs.summ,
                aes(xmin = lo.1,
                    xmax = up.1,
                    y = name,
                    color = coef.group),
                height = 0,
                linewidth = 1.35,
                position = position_nudge(y = -0.1)) +
  
  # colors
  scale_color_manual(values = c("gray45", "olivedrab", "firebrick")) +
  scale_fill_manual(values = c("gray45", "olivedrab", "firebrick")) +
  
  # x axis title
  xlab("Parameter estimate") +
  
  # y-axis scale
  scale_y_discrete(limits = rev,
                   labels = rho.labels) +
  
  # coordinates
  coord_cartesian(xlim = c(-2.5, 2.5),
                  ylim = c(1.3, 8)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(t = 1.17, b = 0.21, l = 0, r = 0.1, unit = "cm"),
        legend.position = "none") +
  
  # annotation
  annotate("text", x = 2, y = 8, label = "b")

# 311 x 286

# ______________________________________________________________________________
# 3e. Plot together ----
# ______________________________________________________________________________

cowplot::plot_grid(phi.est, 
                   rho.est, 
                   nrow = 1)

# 510 x 312

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
phi.pred.plot <- ggplot(data = phi.summ) +
  
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
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
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
  
# 511 x 345

# ______________________________________________________________________________
# 4b. Rho ----
# ______________________________________________________________________________

# sex-year-trt predictions
x.rho <- params.2 |> do.call(rbind, args = _) |> as.data.frame()

# CONTROLS
rho.ctrl.1 <- exp(x.rho$rho_b0)
rho.ctrl.2 <- exp(x.rho$rho_b0 + x.rho$`rho_b[1]`)
rho.ctrl.3 <- exp(x.rho$rho_b0 + x.rho$`rho_b[2]`)

# RET
rho.ret.2 <- exp(x.rho$rho_b0 + x.rho$`rho_b[1]` + x.rho$`rho_b[3]`)
rho.ret.3 <- exp(x.rho$rho_b0 + x.rho$`rho_b[2]` + x.rho$`rho_b[4]`)

# PIL
rho.pil.2 <- exp(x.rho$rho_b0 + x.rho$`rho_b[1]` + x.rho$`rho_b[5]`)
rho.pil.3 <- exp(x.rho$rho_b0 + x.rho$`rho_b[2]` + x.rho$`rho_b[6]`)

# summaries
rho.summ <- data.frame(
  
  # medians
  med = c(
    
    # CONTROL
    median(rho.ctrl.1), 
    median(rho.ctrl.2),
    median(rho.ctrl.3),
    
    # RET
    median(rho.ret.2),
    median(rho.ret.3), 
    
    # PIL 
    median(rho.pil.2),
    median(rho.pil.3)
    
  ),
  
  # l90
  l90 = c(
    
    # CONTROL
    as.numeric(hdi(rho.ctrl.1, ci = 0.90))[2], 
    as.numeric(hdi(rho.ctrl.2, ci = 0.90))[2], 
    as.numeric(hdi(rho.ctrl.3, ci = 0.90))[2], 
    
    # RET
    as.numeric(hdi(rho.ret.2, ci = 0.90))[2], 
    as.numeric(hdi(rho.ret.3, ci = 0.90))[2],  
    
    # PIL
    as.numeric(hdi(rho.pil.2, ci = 0.90))[2], 
    as.numeric(hdi(rho.pil.3, ci = 0.90))[2]
    
  ),
  
  # u90
  u90 = c(
    
    # CONTROL
    as.numeric(hdi(rho.ctrl.1, ci = 0.90))[3], 
    as.numeric(hdi(rho.ctrl.2, ci = 0.90))[3], 
    as.numeric(hdi(rho.ctrl.3, ci = 0.90))[3], 
    
    # RET
    as.numeric(hdi(rho.ret.2, ci = 0.90))[3], 
    as.numeric(hdi(rho.ret.3, ci = 0.90))[3], 
    
    # PIL
    as.numeric(hdi(rho.pil.2, ci = 0.90))[3], 
    as.numeric(hdi(rho.pil.3, ci = 0.90))[3]
    
  ),
  
  # identifiers
  period = c(1, 2, 3, 2, 3, 2, 3),
  trt = factor(c(rep("unthinned", 3), rep("retention", 2), rep("piling", 2)),
               levels = c("unthinned", "retention", "piling"))
  
)

# prediction plot
rho.pred.plot <- ggplot(data = rho.summ) +
  
  theme_bw() +
  
  facet_grid(~ trt) +
  
  # connecting lines
  geom_line(aes(x = period,
                y = med),
            linewidth = 0.8) +
  
  # CIs
  geom_errorbar(aes(x = period,
                    y = med,
                    ymin = l90,
                    ymax = u90),
                alpha = 0.35,
                width = 0,
                linewidth = 1.5,
                position = position_dodge(width = 0.9)) +
  
  # points
  geom_point(aes(x = period,
                 y = med),
             shape = 21,
             size = 2,
             fill = "white",
             stroke = 0.8,
             position = position_dodge(width = 0.9)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.27, unit = "cm")) +
  
  # greek letters
  scale_x_continuous(
    
    breaks = 1:3,
    
    labels = c("1" = expression(tau[1]),
               "2" = expression(tau[2]),
               "3" = expression(tau[3])),
    
    limits = c(0.75, 3.25)
    
    ) +
  
  # labels
  ylab("Per capita recruitment") +
  xlab("Transition period")
  
# 523 x 230

# ______________________________________________________________________________
# 4c. Plot together ----
# ______________________________________________________________________________

cowplot::plot_grid(phi.pred.plot, 
                   rho.pred.plot, 
                   nrow = 2)

# 397 x 447

# ______________________________________________________________________________
# 5. Cluster-specific random intercepts ----
# ______________________________________________________________________________
# 5a. Phi ----
# ______________________________________________________________________________

params.phi.clust <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = 1:6) |>
  
  as.data.frame() %>%
  
  mutate(
    
    int1 = inv_logit(phi_b0 + phi_sd * `phi_c[1]`),
    int2 = inv_logit(phi_b0 + phi_sd * `phi_c[2]`),
    int3 = inv_logit(phi_b0 + phi_sd * `phi_c[3]`),
    int4 = inv_logit(phi_b0 + phi_sd * `phi_c[4]`)
    
  ) %>%
  
  # keep intercept columns
  dplyr::select(int1:int4) %>%
  
  # summarize
  pivot_longer(cols = int1:int4) %>%
  
  group_by(name) %>%
  
  summarize(
    
    med = median(value),
    sd = sd(value),
    lo.1 = as.numeric(hdi(value, ci = 0.5)[2]),
    up.1 = as.numeric(hdi(value, ci = 0.5)[3]),
    lo.2 = as.numeric(hdi(value, ci = 0.9)[2]),
    up.2 = as.numeric(hdi(value, ci = 0.9)[3])
    
  )

# random intercept plot
phi.RI <- ggplot() +
  
  theme_classic() +
  
  # 95%
  geom_errorbar(data = params.phi.clust,
                aes(x = name,
                    ymin = lo.2,
                    ymax = up.2),
                alpha = 0.4,
                width = 0,
                linewidth = 1.75) +
  
  # 50% CI
  geom_errorbar(data = params.phi.clust,
                aes(x = name,
                    ymin = lo.1,
                    ymax = up.1),
                width = 0,
                linewidth = 1.75) +
  
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   vjust = 0.5)) +

  ylab("Persistence probability") +
  
  scale_y_continuous(breaks = c(0.25, 0.35, 0.45, 0.55)) +
  
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))

# ______________________________________________________________________________
# 5b. Rho ----
# ______________________________________________________________________________

params.rho.clust <- params.2 |>
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset
  subset(select = 16:21) |>
  
  as.data.frame() %>%
  
  mutate(
    
    int1 = exp(rho_b0 + rho_sd * `rho_c[1]`),
    int2 = exp(rho_b0 + rho_sd * `rho_c[2]`),
    int3 = exp(rho_b0 + rho_sd * `rho_c[3]`),
    int4 = exp(rho_b0 + rho_sd * `rho_c[4]`)
    
  ) %>%
  
  # keep intercept columns
  dplyr::select(int1:int4) %>%
  
  # summarize
  pivot_longer(cols = int1:int4) %>%
  
  group_by(name) %>%
  
  summarize(
    
    med = median(value),
    sd = sd(value),
    lo.1 = as.numeric(hdi(value, ci = 0.5)[2]),
    up.1 = as.numeric(hdi(value, ci = 0.5)[3]),
    lo.2 = as.numeric(hdi(value, ci = 0.9)[2]),
    up.2 = as.numeric(hdi(value, ci = 0.9)[3])
    
  )

# plot 
rho.RI <- ggplot() +
  
  theme_classic() +
  
  # 95%
  geom_errorbar(data = params.rho.clust,
                aes(x = name,
                    ymin = lo.2,
                    ymax = up.2),
                alpha = 0.4,
                width = 0,
                linewidth = 1.75) +
  
  # 50% CI
  geom_errorbar(data = params.rho.clust,
                aes(x = name,
                    ymin = lo.1,
                    ymax = up.1),
                width = 0,
                linewidth = 1.75) +
  
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   vjust = 0.5)) +
  
  ylab("Per capita recruitment") +
  
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))

# plot together
cowplot::plot_grid(phi.RI, rho.RI)

# 351 x 263

# ______________________________________________________________________________
# 6. Write tables ----
# ______________________________________________________________________________

write.table(phi.summary, "clipboard", sep = "\t")
write.table(rho.summary, "clipboard", sep = "\t")
