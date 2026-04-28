# PROJECT: Open multi-session SCR
# SCRIPT: 07a - Visualization (density and sex ratio)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 Apr 2026
# COMPLETED: 27 Apr 2026
# LAST MODIFIED: 28 Apr 2026
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

# loop through each site-year, subsetting individuals and calculating D/SR
all.D.summary <- data.frame()

for (j in 1:48) {
  
  # focal site.year
  focal.col = site.year[j, ]
  
  # summarize (only if we trapped that year)
  if (focal.col$year == 1 & constant.list$first.year[focal.col$site] == 2) {
    
    focal.D.summary <- NULL
    
  } else {
    
    # extract individuals for the focal site
    which.inds <- which(constant.list$site == focal.col$site)
    
    # subset z.only (change to binary)
    focal.z.only <- z.only[ , which.inds + (3102 * (focal.col$year - 1))]
    focal.z.only <- ifelse(focal.z.only == 2, 1, 0)
    
    # subset sex.only (MUST BE multiplied by the presence indicator)
    focal.sex.only <- sex.only[ , which.inds] * focal.z.only

    # posterior density
    focal.D <- rowSums(focal.z.only) / S.areas$area.rect[focal.col$site]
    focal.D.m <- rowSums(focal.sex.only) / S.areas$area.rect[focal.col$site]
    focal.D.f <- (rowSums(focal.z.only) - rowSums(focal.sex.only)) / 
                 S.areas$area.rect[focal.col$site]
    
    # posterior sex ratio index (proportion male) - rowsums
    focal.SR <- rowSums(focal.sex.only) / rowSums(focal.z.only)
    
    # calculate summary
    focal.D.summary <- data.frame(
      
      # indices
      site = focal.col$site,
      year = focal.col$year,
      
      # which parameter?
      param = c("D", "D.m", "D.f", "SR"),
      
      # median
      med = c(median(focal.D), 
              median(focal.D.m),
              median(focal.D.f),
              median(focal.SR)),
      
      # l50
      l50 = c(as.numeric(hdi(focal.D, ci = 0.50)[2]),
              as.numeric(hdi(focal.D.m, ci = 0.50)[2]),
              as.numeric(hdi(focal.D.f, ci = 0.50)[2]),
              as.numeric(hdi(focal.SR, ci = 0.50)[2])),
      
      # u50
      u50 = c(as.numeric(hdi(focal.D, ci = 0.50)[3]),
              as.numeric(hdi(focal.D.m, ci = 0.50)[3]),
              as.numeric(hdi(focal.D.f, ci = 0.50)[3]),
              as.numeric(hdi(focal.SR, ci = 0.50)[3])),
      
      # l90
      l90 = c(as.numeric(hdi(focal.D, ci = 0.90)[2]),
              as.numeric(hdi(focal.D.m, ci = 0.90)[2]),
              as.numeric(hdi(focal.D.f, ci = 0.90)[2]),
              as.numeric(hdi(focal.SR, ci = 0.90)[2])),
      
      # u90
      u90 = c(as.numeric(hdi(focal.D, ci = 0.90)[3]),
              as.numeric(hdi(focal.D.m, ci = 0.90)[3]),
              as.numeric(hdi(focal.D.f, ci = 0.90)[3]),
              as.numeric(hdi(focal.SR, ci = 0.90)[3]))
      
    )
    
  }
  
  # bind in
  all.D.summary <- rbind(all.D.summary, focal.D.summary)
  
}

#_______________________________________________________________________________
# 3b. Write to table ----
#_______________________________________________________________________________

all.D.table <- all.D.summary %>%
  
  pivot_wider(names_from = param,
              values_from = med:u90) %>%
  
  mutate(
    
    Cluster = case_when(
      
      site %in% c(1, 2, 3) ~ 1,
      site %in% c(4, 5, 6) ~ 2,
      site %in% c(7, 8, 9) ~ 3,
      site %in% c(10, 11, 12) ~ 4
      
    ),
    
    Treatment = case_when(
      
      site %in% c(1, 5, 8, 10) ~ "R",
      site %in% c(2, 4, 7, 11) ~ "P",
      site %in% c(3, 6, 9, 12) ~ "C"
      
    ),
    
    Year = case_when(
      
      year == 1 ~ "PRE 1",
      year == 2 ~ "PRE 2",
      year == 3 ~ "POST 1",
      year ==4 ~ "POST 2"
      
    )
    
  ) %>%
  
  # column order
  dplyr::select(
    
    Cluster,
    Treatment,
    Year,
    med_D, l50_D, u50_D, l90_D, u90_D,
    med_D.m, l50_D.m, u50_D.m, l90_D.m, u90_D.m,
    med_D.f, l50_D.f, u50_D.f, l90_D.f, u90_D.f,
    med_SR, l50_SR, u50_SR, l90_SR, u90_SR
    
  ) %>% 
  
  # arrange
  arrange(Cluster, Treatment, Year)

# write to clipboard
write.table(all.D.table, "clipboard", sep = "\t")

#_______________________________________________________________________________
# 4. Prepare for plotting ----
#_______________________________________________________________________________

plot.df <- all.D.summary %>%
  
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
# 5. Plot ----
#_______________________________________________________________________________
# 5a. Density ----

plot.df.D <- plot.df %>% filter(param == "D") 

#_______________________________________________________________________________

# 3 x 4 - gridded facets, points and CIs, connecting lines
ggplot(data = plot.df.D) +
  
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
# 5b. Sex-specific density ----

plot.df.Dsex <- plot.df %>% 
  
  filter(param %in% c("D.m", "D.f")) %>%
  
  mutate(Sex = ifelse(param == "D.m", "M", "F"))

#_______________________________________________________________________________

# 3 x 4 - gridded facets, points and CIs, connecting lines
ggplot(data = plot.df.Dsex) +
  
  theme_bw() +
  
  facet_grid(trt ~ clust) +
  
  # line
  geom_line(aes(x = year,
                y = med,
                color = Sex),
            linewidth = 0.8) +
  
  # 90% CI
  geom_errorbar(aes(x = year,
                    ymin = l90,
                    ymax = u90,
                    color = Sex),
                width = 0,
                linewidth = 1.25,
                alpha = 0.45) +
  
  # median
  geom_point(aes(x = year,
                 y = med,
                 color = Sex),
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
        legend.position = "top") +
  
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
  scale_color_manual(values = c("#FF3300", "gray25"))

# 493 x 442

#_______________________________________________________________________________
# 5c. Sex ratio (proportion male) ----

plot.df.SR <- plot.df %>% filter(param == "SR")

#_______________________________________________________________________________

# 3 x 4 - gridded facets, points and CIs, connecting lines
ggplot(data = plot.df.SR) +
  
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
