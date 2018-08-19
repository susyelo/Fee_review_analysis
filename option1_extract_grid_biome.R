## Option 1. Using the raster package
biome_shp <- shapefile("./data/maps/Olson_processed/Biomes_olson_projected.shp")


## Using Richness raster as template raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")


## Occurrence data
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

## Create species, grid ID, biome dataframe
system.time({
spOcc_geo$grid_id <- raster::extract(r_ref,spOcc_geo,factors=TRUE)
spOcc_geo$Biomes <- raster::extract(biomes_ras,spOcc_geo,factors=TRUE)
})

cell_nas <- which(is.na(spOcc_geo$Biomes))


# then take the raster value with lowest distance to point AND non-NA value in the raster
cell_nearest <- function(xy) {
  
  biomes_ras@data@values[which.min(replace(distanceFromPoints(biomes_ras, xy), is.na(biomes_ras), NA))]

}

# This takes around two hours to finish
system.time({
sampled = apply(X = coordinates(spOcc_geo[cell_nas,]), MARGIN = 1, FUN = cell_nearest)
})

spOcc_geo$Biomes[cell_nas] <- sampled

write_rds(spOcc_geo,"./outputs/02_Species_grid_id_biomes_df.rds")

### Example with the sf package
## Example 
grid_tmp <- st_make_grid(biome_shp, 
                       cellsize = c(1000000, 1000000)) %>%
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid_tmp) %>% cbind(st_coordinates(.))  

geo_tmp <- spOcc_geo[grep("Quercus",spOcc_geo$scrubbed_species_binomial),]

ggplot() + 
  geom_sf(data = biome_shp, fill = 'white', lwd = 0.05) + 
  geom_sf(data = grid_tmp, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = geo_tmp, color = "red") +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")

# which grid square is each point in?
sp_grid_biomes = geo_tmp %>% 
  st_join(grid_tmp, join = st_intersects) %>% 
  st_join(biome_shp, join = st_intersects)
