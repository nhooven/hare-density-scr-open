# PROJECT: Open multi-session SCR
# SCRIPT: 07a - Visualization (density and sex ratio)
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
# 2. Read in data ----
#_______________________________________________________________________________

# density
params.1 <- readRDS("final_samples/samples_clean_1.rds")

# sex ratio
params.SR <- readRDS("final_samples/samples_clean_SR.rds")

# constants
constant.list <- readRDS("for_model/constants.rds")

# S areas
S.areas <- readRDS("for_model/S_area.rds")

#_______________________________________________________________________________
# 3. Calculations ----

site.year <- expand.grid(site = 1:12, year = 1:4)

#_______________________________________________________________________________
# 3a. Density ----
#_______________________________________________________________________________

N.only <- params.1 |> 
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset to just N
  subset(select = 1:48)

# loop through columns
D.mat <- matrix(data = NA, nrow = nrow(N.only), ncol = ncol(N.only))

for (j in 1:ncol(N.only)) {
  
  # colname
  focal.col = site.year[j, ]
  
  # calculate D (only if we trapped that year)
  if (focal.col$year == 1 & constant.list$first.year[focal.col$site] == 2) {
    
    focal.D <- matrix(data = NA, nrow = nrow(N.only), ncol = 1)
    
  } else {
    
    focal.D <- N.only[, j] / S.areas$area.rect[focal.col$site]
    
  }
  
  # bind in
  D.mat[, j] <- focal.D
  
}

# 90% HPD
all.D.summary <- data.frame()

for (j in 1:ncol(N.only)) {
  
  # colname
  focal.col = site.year[j, ]
  
  # summarize (only if we trapped that year)
  if (focal.col$year == 1 & constant.list$first.year[focal.col$site] == 2) {
    
    focal.D.summary <- NULL
    
  } else {
    
    focal.D.summary <- data.frame(
      
      site = focal.col$site,
      year = focal.col$year,
      med = median(D.mat[ ,j]),
      l50 = as.numeric(hdi(D.mat[ ,j], ci = 0.50)[2]),
      u50 = as.numeric(hdi(D.mat[ ,j], ci = 0.50)[3]),
      l90 = as.numeric(hdi(D.mat[ ,j], ci = 0.90)[2]),
      u90 = as.numeric(hdi(D.mat[ ,j], ci = 0.90)[3])
      
    )
    
  }
  
  # bind in
  all.D.summary <- rbind(all.D.summary, focal.D.summary)
  
}

# write to clipboard
write.table(all.D.summary, "clipboard", sep = "\t")

#_______________________________________________________________________________
# 3b. Sex ratio ----
#_______________________________________________________________________________

sex.only <- params.SR |> 
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset to just N
  subset(select = 49:3150)  # last sex

z.only <- params.SR |> 
  
  # bind together
  do.call(rbind, args = _) |>
  
  # subset to just N
  subset(select = 3151:ncol(params.SR[[1]]))

# ensure we have the right amount in each
ncol(sex.only)    # 3102 - the total M
ncol(z.only)      # 3102 * 4 = 12408 - total states

# loop through each site-year, subsetting individuals and calculating SR
all.SR.summary <- data.frame()

for (j in 1:48) {
  
  # focal site.year
  focal.col = site.year[j, ]
  
  # summarize (only if we trapped that year)
  if (focal.col$year == 1 & constant.list$first.year[focal.col$site] == 2) {
    
    focal.SR.summary <- NULL
    
  } else {
    
    # extract individuals for the focal site
    which.inds <- which(constant.list$site == focal.col$site)
    
    # subset sex.only
    focal.sex.only <- sex.only[ , which.inds]
    
    # subset z.only
    focal.z.only <- z.only[ , which.inds * focal.col$year]
    
    # posterior sex ratio index (proportion male) - rowsums
    focal.SR <- rowSums(focal.sex.only) / rowSums(focal.z.only)
    
    # calculate summary
    focal.SR.summary <- data.frame(
      
      site = focal.col$site,
      year = focal.col$year,
      med = median(focal.SR),
      l50 = as.numeric(hdi(focal.SR, ci = 0.50)[2]),
      u50 = as.numeric(hdi(focal.SR, ci = 0.50)[3]),
      l90 = as.numeric(hdi(focal.SR, ci = 0.90)[2]),
      u90 = as.numeric(hdi(focal.SR, ci = 0.90)[3])
      
    )
    
  }
  
  # bind in
  all.SR.summary <- rbind(all.SR.summary, focal.SR.summary)
  
}

# write to clipboard
write.table(all.SR.summary, "clipboard", sep = "\t")

#_______________________________________________________________________________
# 4. Prepare for plotting ----
#_______________________________________________________________________________
# 4a. Density ----
#_______________________________________________________________________________

D.plot.df <- all.D.summary %>%
  
  mutate(trt = case_when(site %in% c(3, 6, 9, 12) ~ "control",
                         site %in% c(1, 5, 8, 10) ~ "retention",
                         site %in% c(2, 4, 7, 11) ~ "piling"),
         clusterID = case_when(site %in% c(1:3) ~ 1,
                               site %in% c(4:6) ~ 2,
                               site %in% c(7:9) ~ 3,
                               site %in% c(10:12) ~ 4)) %>%
  
  # reorder factors and add life zone
  mutate(trt = factor(trt, levels = c("control", "retention", "piling")),
         lz = ifelse(clusterID == 4, "XMC", "SFL"),
         clust = factor(clusterID, labels = c("cluster 1",
                                              "cluster 2",
                                              "cluster 3",
                                              "cluster 4")))

#_______________________________________________________________________________
# 4b. Sex ratio ----
#_______________________________________________________________________________

SR.plot.df <- all.SR.summary %>%
  
  mutate(trt = case_when(site %in% c(3, 6, 9, 12) ~ "control",
                         site %in% c(1, 5, 8, 10) ~ "retention",
                         site %in% c(2, 4, 7, 11) ~ "piling"),
         clusterID = case_when(site %in% c(1:3) ~ 1,
                               site %in% c(4:6) ~ 2,
                               site %in% c(7:9) ~ 3,
                               site %in% c(10:12) ~ 4)) %>%
  
  # reorder factors and add life zone
  mutate(trt = factor(trt, levels = c("control", "retention", "piling")),
         lz = ifelse(clusterID == 4, "XMC", "SFL"),
         clust = factor(clusterID, labels = c("cluster 1",
                                              "cluster 2",
                                              "cluster 3",
                                              "cluster 4")))

#_______________________________________________________________________________
# 4c. Both ----
#_______________________________________________________________________________

both.plot.df <- D.plot.df %>%
  
  # drop 50%
  dplyr::select(-c(l50, u50)) %>%
  
  # rename
  rename(
    
    D.med = med,
    D.l90 = l90,
    D.u90 = u90
    
  ) %>%
  
  # join
  
  left_join(
    
    SR.plot.df %>%
      
      # drop 50%
      dplyr::select(-c(l50, u50)) %>%
      
      # rename
      rename(
        
        SR.med = med,
        SR.l90 = l90,
        SR.u90 = u90
        
      )
    
  )

#_______________________________________________________________________________
# 5. Plot ----
#_______________________________________________________________________________
# 5a. Density ----
#_______________________________________________________________________________

# 3 x 4 - gridded facets, points and CIs, connecting lines
ggplot(data = D.plot.df) +
  
  theme_bw() +
  
  facet_grid(trt ~ clust) +
  
  # line
  geom_line(aes(x = year,
                y = med,
                color = lz),
            linewidth = 0.8) +
  
  # 90% CI
  geom_errorbar(aes(x = year,
                    ymin = l90,
                    ymax = u90,
                    color = lz),
                width = 0,
                linewidth = 1.25,
                alpha = 0.45) +
  
  # median
  geom_point(aes(x = year,
                 y = med,
                 color = lz),
             shape = 21,
             size = 1.25,
             stroke = 0.8,
             fill = "white") +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = 0.5,
                                   size = 8),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "gray90",
                                        linetype = "blank"),
        legend.position = "none") +
  
  # axis title
  ylab("Density (hares/ha)") +
  
  # coordinates
  # let's add some padding on either side
  coord_cartesian(xlim = c(0.75, 4.25)) +
  
  # y axis
  scale_y_continuous(breaks = c(0.5, 1.5, 2.5)) +
  
  # x labels
  scale_x_continuous(breaks = seq(1, 4, 1),
                     labels = c("PRE 1", "PRE 2", "POST 1", "POST 2")) +
  
  # color
  scale_color_manual(values = c("#003300", "#669900"))

# 493 x 403

#_______________________________________________________________________________
# 5b. Sex ratio (proportion male) ----
#_______________________________________________________________________________

# 3 x 4 - gridded facets, points and CIs, connecting lines
ggplot(data = SR.plot.df) +
  
  theme_bw() +
  
  facet_grid(trt ~ clust) +
  
  # horizontal line at parity
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "gray50") +
  
  # line
  geom_line(aes(x = year,
                y = med,
                color = lz),
            linewidth = 0.8) +
  
  # 90% CI
  geom_errorbar(aes(x = year,
                    ymin = l90,
                    ymax = u90,
                    color = lz),
                width = 0,
                linewidth = 1.25,
                alpha = 0.45) +
  
  # median
  geom_point(aes(x = year,
                 y = med,
                 color = lz),
             shape = 21,
             size = 1.25,
             stroke = 0.8,
             fill = "white") +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = 0.5,
                                   size = 8),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "gray90",
                                        linetype = "blank"),
        legend.position = "none") +
  
  # axis title
  ylab("Proportion male") +
  
  # coordinates
  # let's add some padding on either side
  coord_cartesian(xlim = c(0.75, 4.25)) +
  
  # y axis
  #scale_y_continuous(breaks = c(0.25, 0.35, 0.45, 0.55)) +
  
  # x labels
  scale_x_continuous(breaks = seq(1, 4, 1),
                     labels = c("PRE 1", "PRE 2", "POST 1", "POST 2")) +
  
  # color
  scale_color_manual(values = c("#003300", "#669900"))

# 493 x 403

#_______________________________________________________________________________
# 5c. Both correlation ----
#_______________________________________________________________________________

ggplot(data = both.plot.df) +
  
  theme_classic() +
  
  # errors
  geom_errorbar(aes(x = D.med,
                    ymin = SR.l90,
                    ymax = SR.u90,
                    color = lz),
                width = 0,
                linewidth = 1,
                alpha = 0.25) +
  
  geom_errorbar(aes(y = SR.med,
                    xmin = D.l90,
                    xmax = D.u90,
                    color = lz),
                width = 0,
                linewidth = 1,
                alpha = 0.25) +
  
  # point estimates
  geom_point(aes(x = D.med,
                 y = SR.med),
             shape = 21,
             size = 1.25,
             stroke = 0.8,
             fill = "white") +
  
  # theme
  theme(legend.position = "none") +
  
  # axis titles
  xlab("Density (hares/ha)") +
  ylab("Proportion male") +
  
  scale_color_manual(values = c("#003300", "#669900"))

# 300 x 300
