#### Paper Sub-Antarctic Marine Isoscapes - Riccialdelli et al. 2024 - Progresss in Oceanography -under review ####

######################### Merge bathymetry rasters ############################################
##Import all to make one bathymetry raster
##Last modification July 2024 - Sami Dodino ####

rm(list = ls())
ls()

library(terra)
library(sf)
library(dplyr)

# Load bathymetry rasters
bat_IBCSO <- rast("batimetria/Batimetria_IBCSO/IBCO_ARG.tif")

crs(bat_IBCSO)
summary(bat_IBCSO)
ext(bat_IBCSO)

bat_beagle <- rast("batimetria/Canal_Beagle_interpol_reproyectado.tif")
crs(bat_beagle)

# Raster with different projections, reproject to the larger raster
bat_beagle_proj <- project(bat_beagle, crs(bat_IBCSO))

writeRaster(bat_beagle_proj, "Bat_Beagle_reproject.tif", filetype = "GTiff", overwrite = TRUE)

# Mosaic rasters
bat_beagle_proj_again <- rast("batimetria/Bat_Beagle_reproject.tif")

# Resample to match the resolution and extent of the larger raster to the smaller one
r22 <- resample(bat_beagle_proj_again, bat_IBCSO, method = "ngb")

# Create the mosaic
r3 <- mosaic(bat_IBCSO, r22, fun = mean)

# Save the mosaic
writeRaster(r3, filename = "bat_mosaico.nc", filetype = "CDF", overwrite = TRUE)
writeRaster(r3, filename = "Bat_mosaico.tif", filetype = "GTiff", overwrite = TRUE)

# New bathymetry for each point
# Load lat and long of each station
stations <- read.csv("Nueva_prof/Bat_mosaico_global_total.csv")

# Create an object to project the points
crs.global <- "+proj=longlat +datum=WGS84 +towgs84=0,0,0"

# Plot points
plot_locations <- st_as_sf(stations, coords = c("LONG", "LAT"), crs = crs.global)

ggplot() +
  geom_sf(data = plot_locations)

# Boundary object
boundary <- st_read("Mapa_base/mapa_base.shp")

# Visualization of points and base map
ggplot() +
  geom_sf(data = boundary) +
  geom_sf(data = plot_locations)

st_write(plot_locations, "Shape_estaciones/estaciones_total.shp", driver = "ESRI Shapefile")

# Extract values in a buffer around each point
estaciones_0.1km <- extract(r3, plot_locations, fun = mean, buffer = 10, df = TRUE)
#save
write.csv(estaciones_0.1km, "Nueva_prof/Nueva_prof_global_total.csv", row.names = FALSE)

#### Extract regional bathymetry #####

bat_mosaico <- rast("batimetria/Bat_mosaico.tif")

# Project raster
bat_mosaico_proj <- project(bat_mosaico, crs = crs.global)

# Load lat and long of each station
ptos_regional_SPOM <- read.csv("Nueva_prof/lat-long_regional_SPOM.csv")

# Create an object to project the points
plot_locations <- st_as_sf(ptos_regional_SPOM, coords = c("LONG", "LAT"), crs = crs.global)

# Save 
st_write(plot_locations, "Shape_estaciones/Regional_SPOM.shp", driver = "ESRI Shapefile")

# Extract values in a buffer around each point
regional <- extract(bat_mosaico_proj, plot_locations, fun = mean, buffer = 10, df = TRUE)

# Save CSV
write.csv(regional, "Nueva_prof/Nueva_prof_regional_SPOM.csv", row.names = FALSE)
