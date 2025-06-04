source("Geospatial_Analysis.R")
source("tiling_functions.R")

### Select folder of the region to be studied,the root of the folder, and select
### the name of the band except for "B#" (B1, B2, B3 ~ B7)
prepare_landsat_stack(
  landsat_image_folder = "Nebraska1",
  root = "C:/Projects/Landsat/",
  landsat_root_name = "LC08_L2SP_030039_20240327_20240410_02_T1_SR_"
)
### Tiles the previously made landsat stack file into the specified tile size in the folder
tiling_landsat(
  landsat_stack = landsat_stack, 
  tile_size = 512, 
  root = "C:/Projects/Landsat/",
  output_folder = "Nebraska_Tiles")

# Create a folder to store images within the folder with tiled images
setup_outputFolder(
  root = "C:Nebraska_Tiles",
  output_dir = "tile_analysis_plots",
  file_dir = root
)

# Analyze Tiles

par(mfrow=c(1,1))
pixel_frequency_per_band(na.omit(landsat_stack), 0, 30000)
raster_clustering(landsat_stack$B5)
band_km_centers = raster_classification(landsat_stack$B5, 7)
#Stores the k means centers in a CSV for future reference as k-means is 
#performed on the original landsat image, not the tile. This is computationally expensive.
write.csv(band_km_centers, "Nebraska1_B5_centers.csv", row.names = TRUE)
reclass_df = create_reclass_matrix(band_km_centers)
analyze_landsat_tiles(
  tile_dir = "Nebraska_Tiles",
  band_name = "B5",
  reclass_centers = band_km_centers,
  reclass_df = reclass_df
)
