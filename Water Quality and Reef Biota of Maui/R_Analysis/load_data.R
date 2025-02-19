# load_data.R

library(terra)
library(sf)

load_data <- function(tif_pattern = "20130827.*\\.TIF$") {
  
  # 1. Load and process Landsat raster
  tif_files <- list.files(pattern = tif_pattern, full.names = TRUE)
  landsat_stack <- rast(tif_files)
  names(landsat_stack) <- gsub("LC08_L2SP_063046_20130827_20200913_02_T1_SR_", 
                               "", names(landsat_stack))
  
  # 2. Load and process water chemistry data
  water_chem <- read.csv("Water Chemistry.csv", header = TRUE, 
                         check.names = FALSE, 
                         stringsAsFactors = FALSE)
  water_df <- setNames(water_chem[1:129, ], 
                       c("Sample.Name", "Location", "Sample.Type", "Date", "Time", "Latitude", "Longitude",
                         "Salinity", "P", "N", "PO43", "SiO44", "NO3", "NO2", "NH4", "DIN", "15N.NO3"))
  
  # 2.1 Extract unique locations
  water_long_lat <- unique(water_df[, c("Longitude", "Latitude")])
  selected_indices <- c(1, 24, 39, 69, 92, 109)
  x <- water_long_lat[selected_indices, ]
  
  # 2.2 Convert to spatial object and transform CRS
  x_sf <- st_as_sf(x, coords = c("Longitude", "Latitude"), crs = 4326)
  x_sf <- st_transform(x_sf, crs(landsat_stack))
  
  # 2.3 Create named list of locations
  locations <- setNames(split(st_coordinates(x_sf), seq(nrow(x_sf))),
                        c("Honolua", "Honomanu", "Kahului", "Maalaea", "Kuau", "Waiehu"))
  
  # 3 Load and process benthic data
  benthic = read.csv("benthic_assesment.csv", 
                     header = TRUE, 
                     check.names = TRUE, 
                     stringsAsFactors = FALSE,)
  # Return list of objects
  return(list(landsat = landsat_stack, water_df = water_df, locations = locations, benthic = benthic))
}
