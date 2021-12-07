
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)
library(patchwork)
library(fasterize)

## TODO: Run for 2.5, 5 and 10 arc minutes

res <- 2.5
projection <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs " # ETRS89-extended / LAEA Europe (epsg:3035)
years_ref <- 1986:2015

# Get data ----------------------------------------------------------------

### Disturbance data ###

## The following code uses disturbance maps downloaded from 10.5281/zenodo.3924381 (version1.0), which are downloaded using following code (if not existent)

# TODO: Automatically download disturbance maps

# Paths to disturbance and forest data
path_dist_maps <- "data/disturbances/version1.0/"
path_forest <- "data/disturbances/version1.0/"

# Get countries to be included in the study
countries <- list.files(path_dist_maps)

### Study area ###

## The following code uses global administrative boundareis (https://gadm.org), which are downloaded using following code (if not existent)

# Get names of countries included and matc to ISO3 names
countries <- getData("ISO3") %>%
  mutate(NAME = gsub(" ", "", tolower(NAME))) %>%
  mutate(NAME = case_when(
    NAME == "bosniaandherzegovina" ~ "bosniaherzegovina",
    NAME == "czechrepublic" ~ "czechia",
    TRUE ~ NAME
  )) %>%
  filter(NAME %in% countries)

# If not already done, download shapefiles of each country
dir.create("data/studyregion", showWarnings = FALSE)
if (!file.exists("data/studyregion/gadm36_ALB_0_sp.rds")) {
  countries_sf <- countries %>%
    split(.$NAME) %>%
    map(., ~ getData("GADM", country = .$ISO3, level = 0, path = "data/studyregion")) %>%
    map(st_as_sf) %>%
    do.call(what = rbind, args = .)
} else {
  countries_sf <- list.files("data/studyregion", full.names = TRUE) %>%
    map(readRDS) %>%
    map(st_as_sf) %>%
    do.call(what = rbind, args = .)
}
  
### Climate data ###

## The following code uses BIOCLIM variables (https://www.worldclim.org/data/bioclim.html), which are downloaded using following code (if not existent)

# Download data
if (!file.exists("data/bioclim/wc2-5/bio1.bil") & res == 2.5) {
  getData("worldclim", var = 'bio', res = res, path = "data/bioclim")
}

### TODO: Download for other resolutions

# Create reference grid for analysis based on BIOCLIM and save to disc (don't do it if file already exists) 
if (!file.exists(paste0("data/bioclim/bioclim_reference_grid_", gsub("\\.", "-", res), "_europe.gpkg"))) {
  grid <- raster("data/bioclim/wc2-5/bio1.bil") # Select any layer from CHELSA
  grid_eur <- crop(grid, countries_sf)
  grid_eur[!is.na(grid_eur)] <- 1:ncell(grid_eur[!is.na(grid_eur)])
  writeRaster(grid_eur, "data/bioclim/bioclim_reference_grid_2-5_europe.tif")
  grid_eur_sf <- rasterToPolygons(grid_eur) %>% st_as_sf()
  names(grid_eur_sf) <- c("gridindex", "geometry")
  grid_eur_sf <- grid_eur_sf %>% st_transform(., crs = projection)
  write_sf(grid_eur_sf, paste0("data/bioclim/bioclim_reference_grid_", gsub("\\.", "-", res), "_europe.gpkg"))
} else {
  grid_eur_sf <- read_sf(paste0("data/bioclim/bioclim_reference_grid_", gsub("\\.", "-", res), "_europe.gpkg"))
}

# Aggregate to climate grid -----------------------------------------------

for (i in 1:length(countries$NAME)) {
  
  cntr <- countries$NAME[i]
  
  print(cntr)
  
  disturbance <- raster(paste0(path_dist_maps, cntr, "/disturbance_year_filtered_", cntr, ".tif"))
  forest <- raster(paste0(path_forest, cntr, "/prediction_forestcover_", cntr, ".tif"))
  
  years <- unique(values(disturbance))
  years <- sort(na.omit(years))
  years <- years[years %in% years_ref]
  
  patches <- vector("list", length(years))
  edges <- vector("list", length(years))
  k <- 0
  
  for (y in years) {
    print(paste0("Calculating patches and edges: ", y))
    k <- k + 1  
    disturbance_tmp <- disturbance == y
    patches[[k]] <- clump(disturbance_tmp)
    edges[[k]] <- boundaries(disturbance_tmp)
  }
  
  ext <- as(extent(disturbance), 'SpatialPolygons')
  ext <- st_as_sf(ext)
  st_crs(ext) <- st_crs(grid_eur_sf)
  
  grid_sel <- st_intersection(grid_eur_sf, st_as_sf(ext))
  grid_sel_ras <- fasterize(grid_sel, disturbance, field = "gridindex")
  grid_values <- values(grid_sel_ras)
  
  # Calculate disturbance area per grid cell
  disturbance_area <- data.frame(gridindex = grid_values,
                                  disturbance = values(disturbance),
                                  country = cntr) %>%
    na.omit(.) %>%
    group_by(gridindex, year = disturbance) %>%
    summarize(disturbance_area_ha = n() * 0.09,
              country = unique(country)) %>%
    ungroup(.)
  
  # Calculate number of disturbance patches per grid cell
  disturbance_patches <- patches %>%
    map(., ~ data.frame(gridindex = grid_values,
                        patches = values(.),
                        country = cntr) %>%
          na.omit(.) %>%
          group_by(gridindex) %>%
          summarize(disturbances_patches_n = length(unique(patches)),
                    country = unique(country))) %>%
    set_names(years) %>%
    bind_rows(.id = "year") %>% 
    mutate(year = as.integer(year))
  
  # Calculate number of disturbance edge pixels per grid cell
  disturbance_edges <- edges %>%
    map(., ~ data.frame(gridindex = grid_values,
                        edges = values(.),
                        country = cntr) %>%
          na.omit(.) %>%
          group_by(gridindex) %>%
          summarize(disturbances_edges_n = sum(edges),
                    country = unique(country))) %>%
    set_names(years) %>%
    bind_rows(.id = "year") %>% 
    mutate(year = as.integer(year))
  
  # Caculate forest area per grid cell
  forest <- data.frame(gridindex = grid_values,
                       forest = values(forest),
                       country = cntr) %>%
    filter(!is.na(forest)) %>%
    group_by(gridindex) %>%
    summarize(forest_ha = sum(forest == 1, na.rm = TRUE) * 0.09,
              land_ha = n() * 0.09,
              country = unique(country)) %>%
    ungroup(.)
  
  dat <- disturbance_area %>%
    left_join(disturbance_patches, by = c("gridindex", "year", "country")) %>%
    left_join(disturbance_edges, by = c("gridindex", "year", "country")) %>%
    left_join(forest, by = c("gridindex", "country"))
  
  write_csv(dat, paste0("data/disturbances/aggregated_to_grid/disturbances_aggregated_to_grid_", cntr, ".csv"))
  
}

# Combine for all countires and calculate averages ------------------------

dat <- list.files("data/disturbances/aggregated_to_grid", pattern = ".csv$") %>%
  map(read_csv) %>%
  bind_rows()

dat <- dat %>%
  group_by(gridindex, year) %>%
  summarize(disturbance_area_ha = sum(disturbance_area_ha),
            disturbances_patches_n = sum(disturbances_patches_n),
            disturbances_edges_n = sum(disturbances_edges_n),
            forest_ha = sum(forest_ha),
            land_ha = sum(land_ha)) %>%
  ungroup() %>%
  mutate(disturbance_rate = disturbance_area_ha / forest_ha,
         disturbance_patch_density = disturbances_patches_n / forest_ha,
         disturbance_edge_density = disturbances_edges_n / forest_ha) %>%
  group_by(gridindex) %>%
  summarize_at(.vars = vars(disturbance_rate:disturbance_edge_density), 
               .funs = list(mean = mean, varcof = function(x) sd(x) / mean(x))) %>%
  ungroup()

# Translate into raster ---------------------------------------------------


dat <- read_csv(paste0("data/disturbances/aggregated_to_grid/disturbances_aggregated_to_grid_", cntr, ".csv"))

grid_eur_sf_dat <- grid_eur_sf %>%
  left_join(dat)

if (!exists("grid_eur")) {
  grid_eur <- raster("data/bioclim/bioclim_reference_grid_2-5_europe.tif")
}

fasterize(grid_eur_sf_dat, field = "disturbance_patch_density")
