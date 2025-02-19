## Libraries
library(terra)
library(tmap)
library(sp)
library(sf)
library(gdistance)
library(factoextra)
################################################################################
setwd( #set working directory#
  )
###############################################################################
source("load_data.R")
data_list <- load_data()

# Access loaded data
landsat_stack <- data_list$landsat
water_df <- data_list$water_df
locations <- data_list$locations
benthic <- data_list$benthic
# Example: Plot the raster
plotRGB(landsat_stack, r = "B4", g = "B3", b = "B2", stretch = "lin", axes = TRUE)

# Cropping

### Create area of interest (AOI) to select island of Maoi
aoi = ext(115000.0, 192350.0, 2274000,2330500) 

### Crop the raster brick using the defined aoi
maoi_2013 = crop(landsat_stack, aoi)
maoi_2013 = na.omit(maoi_2013)
par(mfrow=c(1,1))
maoi_2013_wgs84 <- project(maoi_2013, "EPSG:4326")
lon_breaks <- seq(ext(maoi_2013_wgs84)[1], ext(maoi_2013_wgs84)[2], length.out = 3)
lat_breaks <- seq(ext(maoi_2013_wgs84)[3], ext(maoi_2013_wgs84)[4], length.out = 3)
plotRGB(maoi_2013_wgs84,
        r = "B4", g = "B3", b = "B2",
        stretch = "lin",
        axes = TRUE,
        main = "Natural Color",
        smooth = TRUE,
        mar = c(3.1, 3.1, 2.1, 2.1),
        cex = 1.1,
        box = TRUE,
        clip = TRUE
        )
plot(maoi_2013$B6)
### Set up a 2-row, 3-column layout for subplots
par(mfrow = c(2,4))

### Loop through each band and plot the pixel frequency
for (band in names(maoi_2013)) {
  hist(maoi_2013[[band]], breaks = 200, xlim = c(6500, 35000),
       main = paste("Histogram of", band))
}

# Raster sampling, clustering, reclassification.

# Extract and clean B5 band data
landsat_b5 <- na.omit(maoi_2013$B6) # Remove NA values

# Scale values
landsat_b5_scaled <- scale(landsat_b5) 

# Sample data to reduce computational load
landsatb5_sample <- spatSample(landsat_b5_scaled, 
                               size = 10000, 
                               method = "regular", 
                               cells = TRUE, 
                               xy = TRUE)
landsatb5_sample <- na.omit(landsatb5_sample)
# Visualize optimal number of clusters using within-cluster sum of squares (WSS)
fviz_nbclust(landsatb5_sample[4], kmeans, method = "wss")

# Extract values and apply k-means clustering
b5_means <- kmeans(na.omit(values(landsat_b5)), 
                   centers = 5, 
                   nstart = 25)

# Sort cluster centers in ascending order
b5_km_centers <- sort(b5_means$centers[, 1])

# Output cluster centers
b5_km_centers

### Create classification matrix based on the kmeans clustering points.
reclass_df <- c(0, round(b5_km_centers[[1]]), 1,
                round(b5_km_centers[[1]]), round(b5_km_centers[[2]]), 2,
                round(b5_km_centers[[2]]), round(b5_km_centers[[3]]), 3,
                round(b5_km_centers[[3]]), round(b5_km_centers[[4]]), 4,
                round(b5_km_centers[[4]]), Inf, 5)

reclass_m <- matrix(reclass_df,# Creates a reclassification matrix
                    ncol = 3,
                    byrow = TRUE)

### reclassify the raster using the reclass object - reclass_m
maoi <- classify(maoi_2013$B5,
                 reclass_m)
### view number of pixels in each class and visualize the classified raster
par(mfrow=c(1,2))
barplot(maoi,
        main = "Number of pixels in each class")
plot(maoi,
     col = c("red", "blue", "green","yellow", "orange"),
     main = "Reclassified Maoi")
maoi_raster <- raster(maoi)
### Mask the land (set land to NA, water to 1)
maoi_raster[maoi_raster > 1] <- NA  # Land as NA 
maoi_raster[maoi_raster <= 1] <- 1  # Water as 1

# Coastline Extraction & Masking


coastline <- boundaries(maoi_raster, inner = TRUE)
# Mask everything except coastline
coastline[coastline == 0] <- NA
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

# Define the resolution of the raster (match original raster)
raster_template <- rast(maoi_raster)  # Use original raster as template
res(raster_template) <- res(maoi_raster)  # Ensure resolution matches

# Rasterize the buffered coastline
buffered_raster <- rasterize(buffered_coastline, raster_template, field = 1, background = NA)

# Assign conductance values (1 for coastline, NA for non-coastline)
buffered_raster[!is.na(buffered_raster)] <- 1
buffered_raster[is.na(buffered_raster)] <- 0
par(mfrow=c(1,1))
plot(buffered_raster)
# Create a transition matrix (allowing movement along the coastline)
tr <- transition(raster(buffered_raster), transitionFunction=mean, directions=8)
tr <- geoCorrection(tr, type="c")  # Apply correction for true distances


# Ensure points are in the same CRS as the transition matrix
points <- vect(rbind(locations$Honolua, locations$Waiehu), crs=crs(buffered_raster))

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
cost_map <- accCost(tr, adjusted_origin)
par(mfrow=c(1,1))
if (is.finite(extract(cost_map, adjusted_goal))) {
  # Compute shortest path
  shortest_path <- shortestPath(tr, origin = adjusted_origin, goal = adjusted_goal, output = "SpatialLines")
  
  # Plot results
  plot(maoi_2013$B6,
       main = "RGB composite image\n Landsat Bands B5, B4, B3\n\"Infrared ( Vegetation)\"")
  plot(shortest_path, add = TRUE, col = "yellow", lwd = 2)  # Overlay shortest path
  points(adjusted_origin[1], adjusted_origin[2], col = "blue", pch = 16)  # Mark adjusted origin
  points(adjusted_goal[1], adjusted_goal[2], col = "green", pch = 16)  # Mark adjusted goal
} else {
  print("Adjusted goal is unreachable from the adjusted origin.")
}
