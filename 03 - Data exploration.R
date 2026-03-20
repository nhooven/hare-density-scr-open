# PROJECT: Open multi-session SCR
# SCRIPT: 03 - Data exploration
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 20 Mar 2026
# LAST MODIFIED: 20 Mar 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# It will be helpful to examine spatial capture histories, as well as total
# spatial captures, as well as summarizing various information on recaptures, etc.

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(mefa4)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# trap locations
trap.kml <- st_read(dsn = paste0(getwd(), "/spatial_data/", lyr = "final_live.kml"))

# closed CH
ch <- readRDS("for_model/closed_ch.rds")

# covariates
cov <- readRDS("for_model/indiv_covs.rds")

# site lookup table
site.lookup <- read.csv("for_model/site_lookup.csv")

# ______________________________________________________________________________
# 4. Cleaning ----
# ______________________________________________________________________________
# 4a. Trap location data ----
# ______________________________________________________________________________

traps.sf <- trap.kml %>%
  
  # add site name
  mutate(site = rep(c("1A", "1C", "1B",      # 1B and 1C were flipped
                      "2A", "2B", "2C",
                      "3A", "3B", "3C",
                      "4A", "4B", "4C"),
                    each = 36)) %>%
  
  # keep only site and geometry
  dplyr::select(site, geometry) %>%
  
  # arrange by site
  arrange(site) %>%
  
  # add trap number
  mutate(trap = rep(1:36, times = 12)) %>%
  
  # attribute numeric siteID
  left_join(
    
    site.lookup %>% dplyr::select(SiteName, siteID) %>%
      
      rename(site = SiteName) %>%
      
      group_by(site) %>%
      
      slice(1) %>%
      
      ungroup()
    
    ) %>%
  
  # transform to UTM
  st_transform(crs = "epsg:32611")

# ______________________________________________________________________________
# 4b. Add year name to lookup ----
# ______________________________________________________________________________

site.lookup <- site.lookup %>%
  
  mutate(yearName = case_when(
    
    year == 1 ~ "PRE 1",
    year == 2 ~ "PRE 2",
    year == 3 ~ "POST 1",
    year == 4 ~ "POST 2"
    
  )
  
)

# ______________________________________________________________________________
# 5. Plot spatial capture histories by site-year session ----

# Here we want to plot individual centroids as points, and lines connecting them
# to traps if they were caught > 1 time

# this should be extremely general

# ______________________________________________________________________________
# 5a. Function ----
# ______________________________________________________________________________

plot_spatial_ch <- function (
    
  .site = 1,
  .year = 2
  
) {
  
  # _____________________________
  # a. input error catching
  # _____________________________
  
  # did not trap site during year
  if (.site %notin% site.lookup$siteID[site.lookup$year == .year]) {
    
    return(paste0("Error: Did not trap at site ", .site, " during year ", .year))
    
  } 

  # _____________________________
  # b. subset trap data
  # _____________________________
  
  traps.sf.1 <- traps.sf %>% filter(siteID == .site) %>%
    
    # remove unnecessary columns
    dplyr::select(trap, geometry)
  
  # un-spatialize coordinates for plotting
  traps.df <- as.data.frame(cbind(st_coordinates(traps.sf.1)))
  
  # _____________________________
  # c. subset individual data
  # _____________________________
  
  # 03-19-2026
  # in the future it would be nice to attribute MRID / AnimalID here
  # I'll need to output something from script 01
  
  # which individuals occupy the site?
  indivs.site <- which(cov[[1]] == .site)
  
  # subset CH by individual and year
  ch.year <- ch[indivs.site , , .year]
  
  # which individuals were caught > 0 times? 
  # this is nifty
  indivs.site.year <- indivs.site[apply(ch.year, 1, function (x) !all(x == 37, na.rm = T))]
  
  # subset CH to only those caught
  ch.site.year <- ch[indivs.site.year, , .year]
  
  # sex covariate
  sex.site.year = cov[[9]][indivs.site.year]
  
  # _____________________________
  # d. tidy spatial capture histories
  # _____________________________ 
  
  # loop through individuals
  spatial.ch.traps <- data.frame()
  spatial.ch.cent <- data.frame()
  
  for (i in 1:length(indivs.site.year)) {
    
    # unique traps (no NAs)
    indiv.traps.unq <- unique(ch.site.year[i, ])
    indiv.traps.unq <- indiv.traps.unq[!is.na(indiv.traps.unq) & indiv.traps.unq != 37]
    
    # loop through unique traps
    for (j in 1:length(indiv.traps.unq)) {
      
      trap.j <- traps.sf.1 %>% 
        
        filter(trap == indiv.traps.unq[j]) %>% 
        
        # drop trap
        dplyr::select(-trap) %>%
        
        # add individual information
        mutate(indiv = i,
               indiv.trap = paste0(i, "_", j),
               sex = sex.site.year[i],
               type = "trap")
      
      # bind in
      spatial.ch.traps <- rbind(spatial.ch.traps, trap.j)
      
    } # j
    
    # compute centroids
    # if just one trap, use that coordinate
    if (length(spatial.ch.traps$indiv == i) == 1) {
      
      cent.i.sf <- traps.sf.1 %>% filter(trap == indiv.traps.unq[1])
        
        # df
        cent.i <- data.frame(
          
          indiv = i,
          indiv.trap = NA,
          sex = spatial.ch.traps$sex[spatial.ch.traps$indiv == 1][1],
          type = "cent",
          X = st_coordinates(cent.i.sf)[ ,1],
          Y = st_coordinates(cent.i.sf)[ ,2]
          
        )
      
      # bind in
      spatial.ch.cent <- rbind(spatial.ch.cent, cent.i)
      
      # else, compute centroid and add
    } else {
      
      # just take the mean of coordinates
      focal.coords <- st_coordinates(spatial.ch.traps %>% filter(indiv == i))
      
      mean.coords <- apply(focal.coords, 2, mean)
      
      cent.i <- data.frame(
        
        indiv = i,
        indiv.trap = NA,
        sex = spatial.ch.traps$sex[spatial.ch.traps$indiv == 1][1],
        type = "cent",
        X = mean.coords[1],
        Y = mean.coords[2]
        
      )
        
      # bind in
      spatial.ch.cent <- rbind(spatial.ch.cent, cent.i)
      
    }
    
  } # i
  
  # change sex to categorical
  spatial.ch <- spatial.ch %>%
    
    mutate(sex = case_when(sex == 0 ~ "F",
                           sex == 1 ~ "M",
                           is.na(sex) == T ~ "U"))
  
  # un-spatialize for plotting
  spatial.ch.traps.df <- as.data.frame(cbind(spatial.ch.traps[ , 1:4], 
                                             st_coordinates(spatial.ch.traps))) %>%
    
    dplyr::select(-geometry)
  
  # segments for plotting
  # traps
  spatial.ch.seg.trap <- spatial.ch.traps.df %>%
    
    # drop "type"
    dplyr::select(-type) %>%
    
    # keep distinct rows
    distinct() %>%
    
    # rename
    rename(X.end = X,
           Y.end = Y)
    
  # centroids
  suppressMessages(
  
  spatial.ch.seg <- spatial.ch.cent %>%
    
    # keep only indiv and coords
    dplyr::select(c(indiv, X, Y)) %>%
    
    # join in 
    right_join(spatial.ch.seg.trap)
  
  )
  
  # _____________________________
  # e. plot
  # _____________________________ 
  
  ggplot() +
    
    theme_classic() +
    
    # traps (pluses)
    geom_sf(data = traps.sf.1,
            shape = 3,
            size = 3) +
    
    # segments connecting to traps
    geom_segment(data = spatial.ch.seg,
                 aes(x = X,
                     y = Y,
                     xend = X.end,
                     yend = Y.end,
                     color = as.factor(indiv)),
                 linewidth = 0.6) +
    
    # individual centroids
    geom_point(data = spatial.ch.cent,
               aes(x = X,
                   y = Y,
                   fill = as.factor(indiv)),
               color = "black",
               shape = 21,
               size = 2.5) +
    
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    
    # theme
    theme(
      
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
      
    ) +
  
    # UTM coordinates
    coord_sf(datum = st_crs(32611)) +
    
    # add title
    ggtitle(paste0(site.lookup$SiteName[site.lookup$siteID == .site][1],
                   " - ",
                   site.lookup$yearName[site.lookup$year == .year][1]))
  

  }
  
# ______________________________________________________________________________
# 5b. Use function ----
# ______________________________________________________________________________

# PRE 1
plot_spatial_ch(.site = 4, .year = 1)
plot_spatial_ch(.site = 5, .year = 1)
plot_spatial_ch(.site = 6, .year = 1)
plot_spatial_ch(.site = 7, .year = 1)
plot_spatial_ch(.site = 9, .year = 1)

# PRE 2
plot_spatial_ch(.site = 1, .year = 2)
plot_spatial_ch(.site = 2, .year = 2)
plot_spatial_ch(.site = 3, .year = 2)
plot_spatial_ch(.site = 4, .year = 2)
plot_spatial_ch(.site = 5, .year = 2)
plot_spatial_ch(.site = 6, .year = 2)
plot_spatial_ch(.site = 7, .year = 2)
plot_spatial_ch(.site = 8, .year = 2)
plot_spatial_ch(.site = 9, .year = 2)
plot_spatial_ch(.site = 10, .year = 2)
plot_spatial_ch(.site = 11, .year = 2)
plot_spatial_ch(.site = 12, .year = 2)

# POST 1
plot_spatial_ch(.site = 1, .year = 3)
plot_spatial_ch(.site = 2, .year = 3)
plot_spatial_ch(.site = 3, .year = 3)
plot_spatial_ch(.site = 4, .year = 3)
plot_spatial_ch(.site = 5, .year = 3)
plot_spatial_ch(.site = 6, .year = 3)
plot_spatial_ch(.site = 7, .year = 3)
plot_spatial_ch(.site = 8, .year = 3)
plot_spatial_ch(.site = 9, .year = 3)
plot_spatial_ch(.site = 10, .year = 3)
plot_spatial_ch(.site = 11, .year = 3)
plot_spatial_ch(.site = 12, .year = 3)

# POST 2
plot_spatial_ch(.site = 1, .year = 4)
plot_spatial_ch(.site = 2, .year = 4)
plot_spatial_ch(.site = 3, .year = 4)
plot_spatial_ch(.site = 4, .year = 4)
plot_spatial_ch(.site = 5, .year = 4)
plot_spatial_ch(.site = 6, .year = 4)
plot_spatial_ch(.site = 7, .year = 4)
plot_spatial_ch(.site = 8, .year = 4)
plot_spatial_ch(.site = 9, .year = 4)
plot_spatial_ch(.site = 10, .year = 4)
plot_spatial_ch(.site = 11, .year = 4)
plot_spatial_ch(.site = 12, .year = 4)
