# Analysis performed in R

## load_data.R

Opens and stores the files used for the subsequent analyses.
load_data is a function that only takes images with the following  = tif_pattern = "20130827.*\\.TIF$"

It performs 3 actions.

1) Loads and processes Landsat data by loading the .TIF files as SpatRaster objects, followed by removing the verbose language of the titles leaving only the names of the bands as "B1, [...], B7".

2) Loads and processes the "Water Chemistry.csv" file and renames the columns into identifiable indices for R. 
2.1) Uses the data set to extract the unique locations of Maoi used for the analysis.

3) Loads and processes the "Benthic.cvs" file, although nothing was used from this file, it was used for the data analysis in the paper.

4)Returns the list of objects as callable variables: "Landsat - SpatRaster object", "water_df - Water Chemistry.csv", "locations - unique locations where analyses were performed", "benthic - benthic data set".

## Geospatial_Analysis.R

Crops the Landsat image data from 2013 selecting only the island of Maui being studied.
Creates a transition map along the coastline of the island and uses gdistance to calculate the distance between points. To do this it uses 12 functions:

### Function 1 
#### Generic plot for plotting B4, B3, B2 bandwidths and showing natural color.
plot_natural_color <- function(landsat_SpatRaster)

### Function 2 
#### Function for cropping a spatraster object. The AOI is modifiable, but is fixed to crop maoi. 
crop_landsat <- function(landsat_SpatRaster, xmin = 115000.0, xmax = 192350.0, ymin = 2274000, ymax = 2330500)

### Function 3 
#### Calculates the frequency of pixels per band in the spatraster object to evaluate which band to use for kmeans clustering.
pixel_frequency_per_band = function(band)

### Function 4
#### Retrieve a sample of the pixels and determine visually the optimal number of cluster with fviz_nbclust (within cluster sums of squares)
raster_clustering <- function(band)

### Function 5
#### kmeans on the band of interest, using the previously determined amount of optimal number of clusters (center_points)
raster_classification <- function(band,center_points)

### Function 6
#### Creates a reclass matrix using the center points as ranges to classify the groups numerically into 1 ~ length(center_points)
create_reclass_matrix <- function(band_km_centers) 

### Function 7
#### Reclassifies the spatraster with the raster_classification
raster_reclassification <- function(reclass_df, band)

### Function 8
#### Visualizes the pixel frequency for each cluster and visualizes how the new classified spatraster looks.
visualize_reclass <- function(reclass_raster)

### Function 9
#### uses the "raster()" function to rasterize the spatraster, since spatraster objects aren't usable for the following steps.
#### Classifies land and water
rasterize_maoi <- function(SpatRaster)

### Function 10
#### Extracts the coastline from the map. Creates a buffered coastline, to increase the width of the coastline to ensure continuity of the boundary
coastline_extraction <- function(maoi_raster)

### Function 11
#### Creates a transition matrix for the coastline to calculate distances on that boundary.
transition_matrix <- function(maoi_raster, buffered_coastline)

### Function 12
#### Visualizes the shortest path between two points (p1, p2) using the transition matrix, on a colorized map of the AOI
visualize_distance <- function(buffered_raster, transition_layer, p1, p2, band) 
