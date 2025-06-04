library(usmap)
library(terra)
library(sf)
### Function to rasterize the landsat images. Removes verbose from band names,
### resulting in B1 ~ B7. Returns image vegetation index of the landsat stacked, and
### the rasterized version: landsat_stack (bands 1 ~ 7).
prepare_landsat_stack = function(landsat_image_folder,root,landsat_root_name){
  landsat_folder = landsat_image_folder
  root_folder = root
  ### Set working directory to root_folder/landsat_folder 
  setwd(paste0(root_folder, landsat_folder, "\\"))
  ### Load, stack, remove verbose from Landsat raster, and visualize the raster
  tif_files <- list.files(pattern = ".TIF", full.names = TRUE)
  landsat_stack <- rast(tif_files)
  names(landsat_stack) <- gsub(landsat_root_name, 
                               "", names(landsat_stack))
  ### Original Extent
  par(mfrow=c(1,1))
  vis_agro(landsat_stack)
  return(landsat_stack)
}

# Creates an output folder where to store images of the polygonized tiles.
# Designates file directory where to retrieve tiles from for further analysis.
setup_outputFolder = function(root, output_dir_name){
  setwd = root
  output_dir = output_dir_name
  dir.create(output_dir, showWarnings = FALSE)
  file_dir = root
  tif_tile = list.files(path = file_dir, pattern = "\\.tif$", full.names = TRUE)
  return(file_dir, tif_tile)
}
### Vis Agro Function
### Visualizes the landsat image to emphasize vegetation index, r = band5,
### g = band 4, b = band 3
vis_agro  <- function(landsat_stack_image){
  plotRGB(landsat_stack_image,
          r = "B5", g = "B4", b = "B3",
          stretch = "lin",
          axes = TRUE,
          smooth = TRUE,
          mar = c(3.1, 3.1, 2.1, 2.1),
          cex = 1.1,
          box = TRUE,
          clip = TRUE)
}

### Crop agro will crop all the bands associated with the landsat image
### cropping is defined by xmin, ymin and xmax, ymax, as well as the landsat image to be cropped.
### x and y refer to the image coordinates. 
crop_agro <- function(landsat_stack, xmin, xmax, ymin, ymax){
  aoi = ext(xmin, xmax, ymin, ymax)
  landsat_analyzed = crop(landsat_stack, aoi)
  return(landsat_analyzed)
}

### tiling_landsat is a function that converts an entire landsate image composite,
### into 256, 512, or n*2 square tiles, resulting in "Z" amounts of tiles.
### Ultimately the user looks at the tiles and determines which tile will be easiest,
### to classify by kmeans and subsequent QGis classification.

tiling_landsat = function(landsat_stack, tile_size, root, output_folder){
  # Define the tile size (choose one)
  tile_size <- tile_size  # Change to 256 or 512 as needed
  
  # Get dimensions
  nrows <- nrow(landsat_stack)
  ncols <- ncol(landsat_stack)
  
  # Get resolution and image extent
  res_x <- res(landsat_stack)[1]  # Pixel width
  res_y <- res(landsat_stack)[2]  # Pixel height
  
  # Get bounds of the raster
  xmin <- xmin(landsat_stack)
  xmax <- xmax(landsat_stack)
  ymin <- ymin(landsat_stack)
  ymax <- ymax(landsat_stack)
  
  
  # Create a directory to save the tiles
  output_dir <- paste0(root,output_folder)
  dir.create(output_dir, showWarnings = TRUE)
  
  # Loop through and crop images into tiles To automate the fragmentation of the Landsat image 
  # while minimizing overlap, systematically loop through the raster using a grid-based approach.
  # The grid is defined as: tile_size * res_x to step through the X Direction and tile_size * res_y
  # to step through in the Y Direction. Seq() is used to ensure each tile starts immediately after the prior.
  #
  tile_id <- 1
  for (x_start in seq(xmin, xmax - (tile_size * res_x), by = tile_size * res_x)) {
    for (y_start in seq(ymax, ymin + (tile_size * res_y), by = -tile_size * res_y)) {
      
      # Compute the extent for the current tile
      x_end <- x_start + (tile_size * res_x)
      y_end <- y_start - (tile_size * res_y)  # Subtract because Y decreases downward
      
      crop_extent <- ext(x_start, x_end, y_end, y_start)
      
      # Crop the raster
      cropped_raster <- crop(landsat_stack, crop_extent)
      
      # Check percentage of NA values in the tile
      na_ratio <- sum(is.na(values(cropped_raster))) / length(values(cropped_raster))
      
      # Check if the tile contains NA values
      if (na_ratio > 0.01) {
        next  # Skip if more than 1% of pixels are NA
      }
      
      # Save the tile
      writeRaster(cropped_raster, 
                  filename = file.path(output_dir, paste0("tile_", tile_id, ".tif")),
                  filetype = "Gtiff",
                  overwrite=TRUE
      )
      
      tile_id <- tile_id + 1
    }
    print("Tiling complete. Images saved in the landsat_tiles directory.")
  }
}


combine_SpatRaster = function(file_dir, landID){
  
  file_dir = file_dir
  
  tif_list = list.files(path = file_dir, pattern = "\\.TIF$", full.names = TRUE)
  
  landsat_stack = rast(tif_list)
  
  names(landsat_stack) = gsub(landID, "", names(landsat_stack))
  
  return(landsat_stack)
}

visualize_tile_polygon_reclass = function(band, filename_no_ext, tile_file){
  tf_reclassified = raster_reclassification(reclass_df, band)
  reclass_df_matrix = matrix(reclass_df, ncol = 3, byrow = TRUE)
  
  par(mfrow = c(1, 3))
  visualize_reclass(tf_reclassified, band_km_centers, filename_no_ext)
  vis_agro(tile_file)
}

analyze_landsat_tiles <- function(tile_dir, output_dir = "tile_analysis_plots", 
                                  band_name = "B5", reclass_centers, reclass_df) {
  # Ensure output directory exists
  dir.create(file.path(tile_dir, output_dir), showWarnings = FALSE, recursive = TRUE)
  
  # Get list of .tif files
  tif_tile <- list.files(path = tile_dir, pattern = "\\.tif$", full.names = TRUE)
  
  for (i in seq_along(tif_tile)) {
    tf <- rast(tif_tile[i])
    
    # Check band existence
    if (!(band_name %in% names(tf))) {
      warning(paste("Band", band_name, "not found in", tif_tile[i]))
      next
    }
    
    band <- tf[[band_name]]
    filename <- basename(tif_tile[i])
    filename_no_ext <- tools::file_path_sans_ext(filename)
    
    output_filename <- file.path(tile_dir, output_dir, paste0(filename_no_ext, "_analyzed.png"))
    
    png(filename = output_filename, width = 4200, height = 2400)
    par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
    
    # This assumes the function `visualize_tile_polygon_reclass()` is in your environment
    visualize_tile_polygon_reclass(band, filename_no_ext, tf)
    
    dev.off()
  }
}
