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
Creates a transition map along the coastline of the island and uses gdistance to calculate the distance between points.
