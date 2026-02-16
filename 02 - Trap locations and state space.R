# PROJECT: Closed multi-session SCR
# SCRIPT: 02 - Trap locations and state space
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 14 Jan 2026
# COMPLETED: 14 Jan 2026
# LAST MODIFIED: 16 Feb 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Purpose and conventions ----
# ______________________________________________________________________________

# Here I'll read in the trap locations (xj),
# keeping a shapefile and a .csv of x-y coordinates.
# Then, I'll prescribe a rectangular state space
# so that I can use the unif(x, y) prior on ACs si
# I'll save the xlim and ylim for each unit's S
# assume constant across years (reasonable)

# I'll also create tighter S for each unit, and a binary "mask"
# to better limit possible ACs on irregular-shaped units
# ultimately, it would be nice if this:
  # (1) improved sampling speed
  # (2) reduced data augmentation needs

# ______________________________________________________________________________
# 2. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(terra)

# ______________________________________________________________________________
# 3. Read in data ----
# ______________________________________________________________________________

# trap locations
trap.kml <- st_read(dsn = paste0(getwd(), "/spatial_data/", lyr = "final_live.kml"))

# unit boundaries (for visualization)
units <- st_read(dsn = paste0(getwd(), "/spatial_data/units_fixed_utm/", lyr = "units_fixed_utm.shp"))

# ______________________________________________________________________________
# 4. Clean trap location data ----
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
  
  # transform to UTM
  st_transform(crs = "epsg:32611")

# ______________________________________________________________________________
# 5. Clean unit boundaries ----
# ______________________________________________________________________________

units.sf <- units %>%
  
  # keep only necessary columns
  dplyr::select(name, geometry) %>%
  
  # rename name to site
  rename(site = name) %>%
  
  # arrange by site
  arrange(site)

# ______________________________________________________________________________
# 6. Plot both ----

library(cowplot)

# ______________________________________________________________________________

all.units <- list()

for (i in 1:12) {
  
  focal.unit <- units.sf %>% slice(i)
  focal.traps <- traps.sf %>% filter(site == focal.unit$site)
  
  all.units[[i]] <- ggplot() +
    
    theme_bw() +
    
    geom_sf(data = focal.unit) +
    
    geom_sf(data = focal.traps,
            shape = 21) +
    
    coord_sf(datum = st_crs(32611)) +
    
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
}

plot_grid(plotlist = all.units, nrow = 4)

# ______________________________________________________________________________
# 7. State space ----

# buffering the bounding box of the trap grid seems the most reasonable here
# we'll pack these into a list with lapply
traps.split <- split(traps.sf, traps.sf$site)

# from preliminary results using the secr package, it looks like we get
# diminishing returns beyond about 125 m
# let's do 150 to be safe

# unsurprisingly, this is ~5 sigma (~32 m)

# ______________________________________________________________________________
# 7a. Rectangular state space ----

# units are irregularly shaped, so many of them have annoying voids

# ______________________________________________________________________________

make_state_space_rect <- function (
  
  x,  
  buffer = 150
  
) {
  
  # bounding box of traps
  focal.bbox <- spatialEco::bbox_poly(x)
  
  # buffer and extract THAT bbox to make sure it's actually a rectangle
  focal.S <- spatialEco::bbox_poly(st_buffer(focal.bbox, dist = buffer))
  
  return(focal.S)
  
}

# use function
all.S.rect <- lapply(traps.split, make_state_space_rect)

# sanity check
all.S.rect.plots <- list()

for (i in 1:12) {
  
  focal.unit <- units.sf %>% slice(i)
  focal.traps <- traps.sf %>% filter(site == focal.unit$site)
  
  all.S.rect.plots[[i]] <- ggplot() +
    
    theme_bw() +
    
    geom_sf(data = all.S.rect[[i]]) +
    
    geom_sf(data = focal.unit) +
    
    geom_sf(data = focal.traps,
            shape = 21) +
    
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
}

plot_grid(plotlist = all.S.rect.plots, nrow = 4)

# ______________________________________________________________________________
# 7b. Irregular state space ----

# we also need to turn these polygons into binary "masks"
# we'll start with terra rasters, then convert to matrices

# COORDINATES
# these need to correspond to the indices of the matrix [row, col]
# so this needs to be (y, x), descending
# meaning that (0, 0) has to be TOP LEFT


# pack traps and rectangular S into a list
traps.rect.S <- list() 

for (i in 1:12) {
  
  traps.rect.S[[i]] <- list(traps.split[[i]], all.S.rect[[i]])
  
}

# ______________________________________________________________________________

make_state_space_irr <- function (
    
  x,    # traps.rect.S
  buffer = 150
  
) {
  
  # buffer traps and dissolve
  focal.poly <- spatialEco::sf_dissolve(st_buffer(x[[1]], dist = buffer))
  
  # create raster
  
  # in a nimbleFunction, we'll need to evaluate the mask using a 1s vector
  # Turek and Eacker do this:
    # hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)]
  # implying that that indices of the mask + 1 correspond to the possible integers
  # of S
  # to do it, we need the TOP LEFT to be (0, 0) instead of the centroid
  # the x-axis is the SAME (the columns)
  # the y-coords need to be TRANSLATED over the x-axis for the spatial layers
  
  # bounding box
  focal.bbox <- st_bbox(x[[2]])
  
  # subtract coordinates
    # x = j
    # y = i
  # rectangular S
  new.rect <- matrix(c(st_coordinates(x[[2]])[ , 1] - focal.bbox[1],
                       abs(st_coordinates(x[[2]])[ , 2] - focal.bbox[4])),
                     ncol = 2)
  
  # traps
  new.traps <- matrix(c(st_coordinates(x[[1]])[ , 1] - focal.bbox[1],
                        abs(st_coordinates(x[[1]])[ , 2] - focal.bbox[4])),
                      ncol = 2)
  
  # mask
  new.mask <- matrix(c(st_coordinates(focal.poly)[ , 1] - focal.bbox[1],
                       abs(st_coordinates(focal.poly)[ , 2] - focal.bbox[4])),
                     ncol = 2)
  
  # new sf objects
  # rect S
  new.rect.sf <- st_as_sf(as.data.frame(new.rect),
                          coords = c("V1", "V2")) %>%
    
    summarize(geometry = st_combine(geometry)) %>%
    
    st_cast("POLYGON")
  
  # traps
  new.traps.sf <- st_as_sf(as.data.frame(new.traps),
                           coords = c("V1", "V2"))
  
  # mask
  new.mask.sf <- st_as_sf(as.data.frame(new.mask),
                          coords = c("V1", "V2")) %>%
    
    summarize(geometry = st_combine(geometry)) %>%
    
    st_cast("POLYGON")
  
  #plot(st_geometry(new.rect.sf))
  #plot(st_geometry(new.mask.sf), add = T)
  #plot(st_geometry(new.traps.sf), add = T)
  
  # create matrix
  # let's see what a 10-m resolution looks like
  x.width = focal.bbox[3] - focal.bbox[1]
  y.width = focal.bbox[4] - focal.bbox[2]
  
  # round up to make sure these are integers
  mat.ncol = round(x.width / 10)
  mat.nrow = round(y.width / 10)
  
  mask.mat <- matrix(0, nrow = mat.nrow, mat.ncol)
  
  # now we need to fill with 1s as per the mask
  mask.rast <- rast(resolution = 10,
                    xmin = st_bbox(new.rect.sf)[1],
                    xmax = st_bbox(new.rect.sf)[3],
                    ymin = st_bbox(new.rect.sf)[2],
                    ymax = st_bbox(new.rect.sf)[4],
                    nrow = mat.nrow,
                    ncol = mat.ncol,
                    vals = 1)
  
  # mask it
  mask.rast.1 <- mask(mask.rast, vect(new.mask.sf))
  
  # and merge with a zero raster
  mask.rast.0 <- rast(mask.rast, vals = 0)
  
  mask.rast.merge <- merge(mask.rast.1, mask.rast.0)
  
  # and convert to matrix
  # importantly, this matrix is "upside down" so the [row, col] are properly indexed
  mask.mat.final <- matrix(mask.rast.merge,
                           nrow = mat.nrow,
                           ncol = mat.ncol,
                           byrow = T)
  
  # pack into list and return
  all.out <- list(
    
    # rectangular S poly
    new.rect.sf,
    
    # traps
    st_coordinates(new.traps.sf),
    
    # irregular S matrix
    mask.mat.final,
    
    # irregular S poly
    new.mask.sf
    
  )
  
  return(all.out)
  
}

# use function
all.S.irr <- lapply(traps.rect.S, make_state_space_irr)

# ______________________________________________________________________________
# 8. Center S and traps on (0, 0) ----

# UNUSED
# even if I don't end up using the mask, we can keep that transform of the coordinates

# ______________________________________________________________________________

all.S.center <- list()
traps.center <- list()

for (i in 1:12) {
  
  focal.S <- all.rect.S[[i]]
  focal.traps <- traps.split[[i]]
  
  # center coordinates (what will end up being (0, 0))
  cent.x = (max(st_coordinates(focal.S)[ , 1]) + min(st_coordinates(focal.S)[ , 1])) / 2
  cent.y = (max(st_coordinates(focal.S)[ , 2]) + min(st_coordinates(focal.S)[ , 2])) / 2
  
  # S
  # center S coords
  S.coords.uncenter <- st_coordinates(focal.S)[ , c(1:2)]
    
  S.coords.center <- matrix(c(S.coords.uncenter[ , 1] - cent.x,
                              S.coords.uncenter[ , 2] - cent.y),
                            ncol = 2)
  
  # promote to sf
  all.S.center[[i]] <- st_as_sf(as.data.frame(S.coords.center),
                                coords = c("V1", "V2")) %>%
    
    summarize(geometry = st_combine(geometry)) %>%
    
    st_cast("POLYGON")
  
  # TRAPS
  # center trap coords
  traps.coords.uncenter <- st_coordinates(focal.traps)[ , c(1:2)]
    
  traps.coords.center <- matrix(c(traps.coords.uncenter[ , 1] - cent.x,
                                  traps.coords.uncenter[ , 2] - cent.y),
                                ncol = 2)
    
  # promote to sf
  traps.center[[i]] <- st_as_sf(as.data.frame(traps.coords.center),
                                coords = c("V1", "V2"))
  
}

# ______________________________________________________________________________
# 9. x and y limits for S ----

# each will be a matrix [site, 2]
# we can then append to an array with depth 2
# we'll then save as an .rds

# ______________________________________________________________________________

# xlim
S.xlim <- matrix(
  
  data = c(
    
    st_bbox(all.S.irr[[1]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[2]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[3]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[4]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[5]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[6]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[7]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[8]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[9]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[10]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[11]][[1]])[c(1, 3)],
    st_bbox(all.S.irr[[12]][[1]])[c(1, 3)]
    
  ),
  
  nrow = 12,
  ncol = 2,
  byrow = T
  
)

# ylim
S.ylim <- matrix(
  
  data = c(
    
    st_bbox(all.S.irr[[1]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[2]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[3]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[4]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[5]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[6]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[7]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[8]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[9]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[10]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[11]][[1]])[c(2, 4)],
    st_bbox(all.S.irr[[12]][[1]])[c(2, 4)]
    
  ),
  
  nrow = 12,
  ncol = 2,
  byrow = T
  
)

# make array
S.lim <- array(data = c(S.xlim, S.ylim),
               dim = c(12, 2, 2))

# write file
saveRDS(S.lim, file = "for_model/S_lim.rds")

# ______________________________________________________________________________
# 12. S areas ----

# each will be in ha

# ______________________________________________________________________________

S.areas.df <- data.frame()

for (i in 1:12) {
  
  # rectangular S
  area.rect <- st_area(all.S.irr[[i]][[1]]) * 0.0001
  
  # irregular S
  area.irr <- st_area(all.S.irr[[i]][[4]]) * 0.0001
  
  # unit
  area.unit <- as.numeric(st_area(units.sf[i, ]) * 0.0001)
  
  # df
  focal.S.areas.df <- data.frame(
    
    unit = i,
    area.rect = area.rect,
    area.irr = area.irr,
    area.unit = area.unit,
    aug.rect = round(area.rect / area.unit),
    aug.irr = round(area.irr / area.unit)
    
  )
  
  # bind in
  S.areas.df <- rbind(S.areas.df, focal.S.areas.df)
  
}

# write file
saveRDS(S.areas.df, "for_model/S_area.rds")

# ______________________________________________________________________________
# 11. Save trap coordinates ----

# each will be a matrix [36, 2]
# we can then append to an array with depth 12
# we'll then save as an .rds

# ______________________________________________________________________________

trap.coords <- array(data = NA, dim = c(36, 2, 12))

for (i in 1:12) {
  
  trap.coords[ , , i] <- all.S.irr[[i]][[2]]
  
}

# write file
saveRDS(trap.coords, file = "for_model/trap_coords.rds")

# ______________________________________________________________________________
# 12. Final sanity check ----
# ______________________________________________________________________________


