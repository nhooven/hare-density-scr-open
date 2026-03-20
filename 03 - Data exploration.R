# PROJECT: Open multi-session SCR
# SCRIPT: 03 - Data exploration
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 
# LAST MODIFIED: 19 Mar 2026
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

plot_spatial_ch <- function (
    
  .site = 1,
  .year = 2
  
) {
  
  # _____________________________
  # a. input error catching
  # _____________________________
  
  # did not trap site during year
  if (.site %notin% site.lookup$siteID[site.lookup$year == .year]) {
    
    return(paste0("Error: Did not trap at site ", site, " during year ", year))
    
  } 

  # _____________________________
  # b. subset trap data
  # _____________________________
  
  traps.sf.1 <- traps.sf %>% filter(siteID == .site) %>%
    
    # remove unnecessary columns
    dplyr::select(trap, geometry)
  
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
  spatial.ch <- data.frame()
  
  for (i in 1:length(indivs.site.year)) {
    
    # unique traps (no NAs)
    indiv.traps.unq <- unique(ch.site.year[i, ])
    indiv.traps.unq <- indiv.traps.unq[!is.na(indiv.traps.unq)]
    
    # loop through unique traps
    for (j in indiv.traps.unq) {
      
      trap.j <- traps.sf.1 %>% 
        
        filter(trap == j) %>% 
        
        # drop trap
        dplyr::select(-trap) %>%
        
        # add individual information
        mutate(indiv = i,
               sex = sex.site.year[i],
               type = "trap")
      
      # bind in
      spatial.ch <- rbind(spatial.ch, trap.j)
      
    } # j
    
    # compute centroids
    # if just one trap, use that coordinate
    if (length(spatial.ch$indiv == i) == 1) {
      
      cent.i <- traps.sf.1 %>% 
        
        filter(trap == indiv.traps.unq[1]) %>% 
        
        # drop trap
        dplyr::select(-trap) %>%
        
        # add individual information
        mutate(indiv = i,
               sex = sex.site.year[i],
               type = "cent")
      
      # bind in
      spatial.ch <- rbind(spatial.ch, cent.i)
      
      # else, compute centroid and add
    } else {
      
      cent.i <- spatial.ch %>% 
        
        filter(indiv == i) %>% 
        
        # centroid, we just need one
        st_centroid() %>%
        slice(1) %>%
        
        # add individual information
        mutate(indiv = i,
               sex = sex.site.year[i],
               type = "cent")
      
      # bind in
      spatial.ch <- rbind(spatial.ch, cent.i)
      
    }
    
  } # i
  
  # change sex to categorical
  spatial.ch <- spatial.ch %>%
    
    mutate(sex = case_when(sex == 0 ~ "F",
                           sex == 1 ~ "M",
                           is.na(sex) == T ~ "U"))
  
  
  

  }
  

plot_spatial_ch(1, 1)







