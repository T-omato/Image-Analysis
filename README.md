# Image-Analysis
This project analyzes Landsat 8 satellite imagery to classify land, extract coastlines, and compute least-cost paths using `terra`, `sf`, and `gdistance`. 

## Technologies Used  
- R (`terra`, `factoextra`, `gdistance`, `sp`)  
- Satellite Imagery (Landsat 8)  
- Spatial Analysis (K-means clustering, raster classification, transition matrix, least cost distance)
- QGIS for further layer classification
- Neural Networks

## Projects

### Impact of Submarine Groundwater Discharge on Marine Water Quality and Reef Biota of Maui

The authors aimed at studying the effects of agricultural leachoff on marine life at different points on the island of Maui. They used ArcGIS to add features and locate distances on the map. My personal objective was perform the geo-imaging on R instead of ArcGIS. 
The project digitizes landsat data from the island on the year the study was performed. It crops, cleans, and projects CRS onto the data. It creates a transition matrix on the island's coastline in order to evaluate quantitatively the effects distances have on benthic life with respect to large leachoff areas. 
The final product of this project shows and stores the distance between points on the coastline. 








