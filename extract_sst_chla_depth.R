#### Paper Sub-Antarctic Marine Isoscapes - Riccialdelli et al. 2024 - Progresss in Oceanography -under review ####

############ Code to extract data from satellite images #############
## Last modification June 2024 - Sami Dodino ####

rm(list=ls())
ls()

library(terra)
library(sf)
library(ggplot2)


# First step: .nc to .tiff
# data downloaded from https://oceandata.sci.gsfc.nasa.gov/api/file_search)

r1 <- rast("Imagenes/Chla/ENERO_2011.L3m_MO_CHL_chlor_a_4km.nc")
writeRaster(r1, "Imagenes/Chla/ENERO_2011.L3m_MO_CHL_chlor_a_4km.tiff", filetype = "GTiff")

r2 <- rast("Imagenes/SST/ENERO_2011.L3m_MO_SST4_sst4_4km.nc")
writeRaster(r2, "Imagenes/SST/ENERO_2011.L3m_MO_SST4_sst4_4km.tiff", filetype = "GTiff")


# Read .tiff files
r_chla <- rast("Imagenes/Chla/ENERO_2011.L3m_MO_CHL_chlor_a_4km.tiff")

r_sst <- rast("Imagenes/SST/ENERO_2011.L3m_MO_SST4_sst4_4km.tiff")

r_depth<-rast("batimetria/Bat_mosaico.tif")

# Stack the rasters
climStack <- c(r_sst, r_chla)

# Read the lat/long positions for each month/year

table<-read.csv("baseline_data_GIS_8Junio2023_isoscapes.csv",header=T)
names(table)
str(table)

#filter per month/year

filter_points <- table %>% 
  group_by(ID) %>% 
  filter(Month2=="ene"& Year=="2011")

filter_points <- filter_points[, c("Long", "Lat")]

# transfor points to sf object

points_sf <- st_as_sf(filter_points, coords = c("Long", "Lat"), crs = 4326)

# Extract values from raster using a buffer of 4km
buffer_4km_chla_sst <- extract(climStack, vect(points_sf), fun = mean, buffer = 4000, df = TRUE)
buffer_4km_depth <- extract(r_depth, vect(points_sf), fun = mean, buffer = 4000, df = TRUE)


# Merge with lat/long positions
bioclimat <- cbind(points, buffer_4km_chla_sst, buffer_4km_depth)

# Save
write.csv(bioclimat, "extract_chla_sst/ene2011_sst_chla_depth.csv", row.names = FALSE)

