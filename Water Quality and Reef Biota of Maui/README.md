# Impact of Submarine Groundwater Discharge on Marine Water Quality and Reef Biota of Maui

## Raster Manipulation & Image Processing

- Reading geospatial data (raster and vector formats)

### Visualization of landsat bands 1 - 6
![Landsat images](ImageAnalysis/Bands_1-6.png)

## Raster Manipulation & Image Processing

- Stacking raster layers to create multispectral images.
- Cropping and masking images (reducing image size, isolating AOIs).
- Extracting features from raster datasets

![Raster Classification](ImageAnalysis/cropping.png)

- Histogram-based analysis for feature identification.
- K-means clustering for land classification

![Raster Classification](ImageAnalysis/Raster_classification.png)

## Visualizing Maui study points

![Locations on Maui studied](ImageAnalysis/Maoi_sites_studied.png)
A - Honolua,
B - Honomanū,
C - Kahului,
D - Māʻalaea,
E - Kūʻau,
F - Waiehu

## Coastline Extraction, Masking & Distance Calculation

- Identifying land-water boundaries using classification.
- Creating masks for land and water.
- Extracting coastline features using spatial operations.
- Buffering coastlines to define zones.
- Rasterizing vector data for further analysis.
- Computing least-cost paths and transition matrices for spatial movement.

![Coastline Extraction, Masking & Distance Calculation](ImageAnalysis/Transition_Map.png)

