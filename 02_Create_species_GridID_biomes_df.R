# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)

# data --------------------------------------------------------------------
# 1. Presence of species in cells
spOcc<-readRDS("./data/base/2018_02_08_BIEN_OccurData.RData")

# 2. Shapefiles
biome_shp <- shapefile("./data/maps/Olson_processed/Biomes_olson_projected.shp")

# 3. Using Richness raster as template raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")


# Cleaning Occurrence data ------------------------------------------------
spOcc <- spOcc %>% 
  select(scrubbed_species_binomial, longitude, latitude) %>% 
  filter ((longitude >= -180 & longitude <= 180) & (latitude >= -90 & latitude <= 90))

## Get rid of duplicate occurrences
dups <- duplicated(spOcc)
spOcc2 <- spOcc[!dups, ] # Around 12807047 duplicates

spOcc_geo <- spOcc2
coordinates(spOcc_geo) <- c("longitude", "latitude")
crs(spOcc_geo) <- CRS("+proj=longlat +ellps=WGS84") # WGS 84

## Rasterize Biomes shapefile 100 km^2
r_ref <- r_Total_Rich
r_ref[] <- 1:ncell(r_ref)

## Rasterize the biomes
# https://cran.r-project.org/web/packages/fasterize/vignettes/using-fasterize.html
biomes_ras <- rasterize(biome_shp, r_ref)

## Make sure coordinates have the same projection
spOcc_geo <- spTransform(spOcc_geo, crs(r_ref))


# Create species, grid ID, biome dataframe --------------------------------
system.time({
  spOcc_geo$grid_id <- raster::extract(r_ref,spOcc_geo,factors=TRUE)
  spOcc_geo$Biomes <- raster::extract(biomes_ras,spOcc_geo,factors=TRUE)
})

## Eliminate the occurrences out of the extent of the biomes raster
spOcc_geo <- spOcc_geo[!is.na(spOcc_geo$grid_id),]

# then take the raster value with lowest distance to point AND non-NA value in the raster
# Classify grid without biome class

cell_nas <- which(is.na(spOcc_geo$Biomes))


cell_nearest <- function(xy) {
  
  biomes_ras@data@values[which.min(replace(distanceFromPoints(biomes_ras, xy), is.na(biomes_ras), NA))]
  
}

# This takes around two hours to finish
system.time({
  sampled = apply(X = coordinates(spOcc_geo[cell_nas,]), MARGIN = 1, FUN = cell_nearest)
})

# Replace the Biomes NA with the ones estimated previously
spOcc_geo$Biomes[cell_nas] <- sampled


## Replace biomes numbers for names
spOcc_geo$Biomes<-biome_shp$biomes[spOcc_geo$Biomes]

# Write dataframe into a rds file
write_rds(spOcc_geo,"./outputs/02_Species_grid_id_biomes_df.rds")
