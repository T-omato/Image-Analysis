## Libraries
library(terra)
library(sf)
library(gdistance)
library(factoextra)
################################################################################
setwd("C:\\Users\\pablo\\OneDrive\\Documents\\Personal\\Projects\\Carbon ating")
###############################################################################
source("load_data.R")
#Generic plot for plotting B4, B3, B2 bandwidths and showing natural color.
plot_natural_color <- function(landsat_SpatRaster){
  plotRGB(landsat_SpatRaster, r = "B4", g = "B3", b = "B2", stretch = "lin", axes = TRUE,
          main = "Landsat Colorized")
}
# Function for cropping a spatraster object. The AOI is modifiable, but is fixed to crop maoi. 
crop_landsat <- function(landsat_SpatRaster, xmin = 115000.0, xmax = 192350.0, ymin = 2274000, ymax = 2330500){
  ### Create area of interest (AOI) to select island of Maoi
  aoi = ext(xmin, xmax, ymin, ymax) 
  
  ### Crop the raster brick using the defined aoi
  maoi_2013 = crop(landsat_SpatRaster, aoi)
  maoi_2013 = na.omit(maoi_2013)
  
  ### Change coordinates to visualize final AOI with correct coordinates
  ### Change the CRS of the SpatRaster with coordinates for EPSG:4326 
  maoi_2013_wgs84 <- project(maoi_2013, "EPSG:4326")
  
  #Set the breaks to be visualized on the map. 
  lon_breaks <- seq(ext(maoi_2013_wgs84)[1], ext(maoi_2013_wgs84)[2], length.out = 3)
  lat_breaks <- seq(ext(maoi_2013_wgs84)[3], ext(maoi_2013_wgs84)[4], length.out = 3)
  
  ###Show the progression of data cropping
  par(mfrow=c(1,3))
  plot_natural_color(landsat_SpatRaster)
  plot(aoi,
       main = "Crop Extent",
       axes = TRUE,
       border = "green")
  plotRGB(maoi_2013_wgs84,
          r = "B4", g = "B3", b = "B2",
          stretch = "lin",
          axes = TRUE,
          smooth = TRUE,
          mar = c(3.1, 3.1, 2.1, 2.1),
          cex = 1.1,
          box = TRUE,
          clip = TRUE,
          main = "Cropped Landsat")
  return(maoi_2013)
}
# Calculates the frequency of pixels per band in the spatraster object to evaluate which band to use for kmeans clustering.
pixel_frequency_per_band = function(band){
  ### Set up a 2-row, 4-column layout for subplots
  par(mfrow = c(2,4))
  
  ### Loop through each band and plot the pixel frequency
  for (band in names(maoi_2013)) {
    hist(maoi_2013[[band]], breaks = 200, xlim = c(6500, 35000),#no pixels detected below 6500
         main = paste("Histogram of", band))
  }
}
# Retrieve a sample of the pixels and determine visually the optimal number of cluster with fviz_nbclust (within cluster sums of squares)
raster_clustering <- function(band){
  # Extract and clean band data
  landsat_band <- na.omit(band)
  # Scale values
  landsat_band_scaled <- scale(landsat_band) 
  # Sample data to reduce computational load
  landsat_band_sample <- spatSample(landsat_band_scaled, 
                                    size = 10000, 
                                    method = "regular", 
                                    cells = TRUE, 
                                    xy = TRUE)
  landsat_band_sample <- na.omit(landsat_band_sample)
  # Visualize optimal number of clusters using within-cluster sum of squares (WSS)
  fviz_nbclust(landsat_band_sample[4], kmeans, method = "wss")
}
# kmeans on the band of interest, using the previously determined amount of optimal number of clusters (center_points)
raster_classification <- function(band,center_points){
  landsat_band <- na.omit(band)
  # Extract values and apply k-means clustering
  band_means <- kmeans(na.omit(values(landsat_band)), 
                       centers = center_points, 
                       nstart = 25)
  # Sort cluster centers in ascending order
  band_km_centers <- sort(band_means$centers[, 1])
  # Output cluster centers
  return(band_km_centers = band_km_centers)
}
# Creates a reclass matrix using the center points as ranges to classify the groups numerically into 1 ~ length(center_points)
create_reclass_matrix <- function(band_km_centers) {
  # Sort centers to ensure correct binning
  centers <- sort(band_km_centers)
  # Initialize reclassification matrix with first range (0 to first center)
  reclass_list <- list(0, round(centers[1]), 1)
  
  for (i in 1:(length(centers) - 1)) {
    reclass_list <- c(reclass_list, 
                      round(centers[i]), round(centers[i + 1]), i + 1)
  }
  
  # Correct final classification step to Inf
  reclass_list[length(reclass_list) - 1] <- Inf  # Replace last upper bound with Inf
  
  # Convert to vector
  reclass_matrix <- c(unlist(reclass_list), ncol = 3, byrow = TRUE)
  
  return(reclass_matrix[-2,])
}
# Reclassifies the spatraster with the raster_classification
raster_reclassification <- function(reclass_df, band){
  reclass_m <- matrix(reclass_df,# Creates a reclassification matrix
                      ncol = 3,
                      byrow = TRUE)
  
  ### reclassify the raster using the reclass object - reclass_m
  maoi <- classify(band,
                   reclass_m)
  return(maoi)
}
# Visualizes the pixel frequency for each cluster and visualizes how the new classified spatraster looks.
visualize_reclass <- function(reclass_raster){
  ### view number of pixels in each class and visualize the classified raster
  par(mfrow=c(1,2))
  barplot(maoi,
          col = terrain.colors(length(band_km_centers)),
          main = "Number of pixels in each class")
  plot(maoi,
       col = terrain.colors(length(band_km_centers)),   #c("red", "blue", "green","yellow", "orange"),
       main = "Reclassified Maoi")
}
# uses the "raster()" function to rasterize the spatraster, since spatraster objects aren't usable for the following steps.
# Classifies land and water
rasterize_maoi <- function(SpatRaster){
  maoi_raster <- raster(maoi)
  ### Mask the land (set land to NA, water to 1)
  maoi_raster[maoi_raster > 1] <- NA  # Land as NA 
  maoi_raster[maoi_raster <= 1] <- 1  # Water as 1
  return(maoi_raster)
}
# Extracts the coastline from the map. Creates a buffered coastline, to increase the width of the coastline to ensure continuity of the boundary
coastline_extraction <- function(maoi_raster){
  
  coastline <- boundaries(maoi_raster, inner = TRUE)
  # Mask everything except coastline
  coastline[coastline == 0] <- NA
  par(mfrow=c(1,1))
  plot(maoi_raster, col=terrain.colors(10), main="Extracted Coastline")  # Background map
  plot(coastline, add=TRUE, col="red", lwd=3)  # Bold red overlay
  
  # Define a buffer distance (adjust as needed, in map units)
  buffer_distance <- 5  # Example: 5 meters
  
  # Convert the coastline raster to a vector 
  coastline_vec <- as.polygons(rast(coastline))
  
  # Apply a buffer around the extracted coastline
  buffered_coastline <- buffer(coastline_vec, width = buffer_distance)
  
  # Plot the original and buffered coastline for comparison
  plot(coastline_vec, main="Original (Black) vs Buffered (Blue) Coastline")
  
  plot(buffered_coastline, add=TRUE, col="blue", border="blue", lwd=2)
  
  return(buffered_coastline)
}
# Creates a transition matrix for the coastline to calculate distances on that boundary.
transition_matrix <- function(maoi_raster, buffered_coastline){
  # Define the resolution of the raster (match original raster)
  raster_template <- rast(maoi_raster)  # Use original raster as template
  res(raster_template) <- res(maoi_raster)  # Ensure resolution matches
  
  # Rasterize the buffered coastline
  buffered_raster <- rasterize(buffered_coastline, raster_template, field = 1, background = NA)
  
  # Assign conductance values (1 for coastline, NA for non-coastline)
  buffered_raster[!is.na(buffered_raster)] <- 1
  buffered_raster[is.na(buffered_raster)] <- 0
  par(mfrow=c(1,1))
  plot(buffered_raster,#Visualize the conductance values
       main = "Conductance values, 1 - Coastline, 0 - Non-Coastline")
  
  # Create a transition matrix (allowing movement along the coastline)
  tr <- transition(raster(buffered_raster), transitionFunction=mean, directions=8)
  tr <- geoCorrection(tr, type="c")  # Apply correction for true distances
  
  return(list(buffered_raster = buffered_raster, transition_layer = tr))
}
#Visualizes the shortest path between two points (p1, p2) using the transition matrix, on a colorized map of the AOI
visualize_distance <- function(buffered_raster, transition_layer, p1, p2, band) {
  # Ensure points are in the same CRS as the transition matrix
  points <- vect(rbind(p1, p2), crs=crs(buffered_raster))
  
  # Convert points to matrix format (gdistance requires this)
  coords <- as.matrix(geom(points)[, c("x", "y")])  # Extract XY coordinates
  
  # Function to find the nearest water cell to a given point
  find_nearest_water <- function(raster, point_coords) {
    # Ensure the coordinates are a matrix with one row
    point_coords <- matrix(point_coords, ncol=2)
    
    # Get the raster value at the point
    cell_value <- extract(raster, point_coords)
    
    # If the point is already in water, return it
    if (!is.na(cell_value) && cell_value == 1) {
      return(point_coords)
    }
    
    # Identify all water cells
    water_cells <- which(values(raster) == 1, arr.ind=TRUE)
    
    # Get coordinates of all water cells
    water_coords <- xyFromCell(raster, water_cells)
    
    # Compute distances from the point to all water cells
    dists <- spDistsN1(as.matrix(water_coords), point_coords, longlat = FALSE)
    
    # Find the nearest water cell
    nearest_index <- which.min(dists)
    
    # Return the coordinates of the nearest water cell
    return(water_coords[nearest_index, , drop=FALSE])
  }
  
  # Identify nearest water cells for origin and goal
  adjusted_origin <- find_nearest_water(raster(buffered_raster), coords[1,])
  adjusted_goal <- find_nearest_water(raster(buffered_raster), coords[2,])
  
  # Ensure these points are inside the valid transition area
  cost_map <- accCost(transition_layer, adjusted_origin)
  
  par(mfrow=c(1,1))
  if (is.finite(extract(cost_map, adjusted_goal))) {
    # Compute shortest path
    shortest_path <- shortestPath(transition_layer, origin = adjusted_origin, goal = adjusted_goal, output = "SpatialLines")
    
    # Plot results
    p1_name <- names(locations)[which(sapply(locations, identical, p1))]
    p2_name <- names(locations)[which(sapply(locations, identical, p2))]
    
    plotRGB(maoi_2013,
            r = "B5", g = "B4", b = "B3",
            main = paste0(p1_name, " as Blue, ", p2_name, " as Green"),
            stretch = "lin",
            legend=paste0("the distance is ",shortest_path))
    plot(shortest_path, add = TRUE, col = "yellow", lwd = 2)  # Overlay shortest path
    points(adjusted_origin[1], adjusted_origin[2], col = "blue", pch = 16)  # Mark adjusted origin
    points(adjusted_goal[1], adjusted_goal[2], col = "green", pch = 16)  # Mark adjusted goal
  } else {
    print("Adjusted goal is unreachable from the adjusted origin.")
  }
}


data_list <- load_data()

# Access loaded data
landsat_stack <- data_list$landsat
water_df <- data_list$water_df
locations <- data_list$locations

crop_landsat(landsat_stack)

pixel_frequency_per_band(landsat_stack)

find_cluster_numbers = raster_clustering(maoi_2013$B6)

band_km_centers = raster_classification(maoi_2013$B6, 5)

reclass_df = create_reclass_matrix(band_km_centers)

maoi = raster_reclassification(reclass_df, maoi_2013$B5)

visualize_reclass(maoi)

maoi_raster = rasterize_maoi(maoi)

buffered_coastline = coastline_extraction(maoi_raster)

tr_matrix = transition_matrix(maoi_raster, buffered_coastline)

visualize_distance(tr_matrix$buffered_raster, 
                   tr_matrix$transition_layer, 
                   locations$Honomanu, 
                   locations$Kahului, 
                   maoi_2013)






