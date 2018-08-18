## Option 1. Using the raster package
biome_shp <- shapefile("./data/maps/Olson_processed/Biomes_olson_projected.shp")

## Occurrence data
spOcc <- spOcc %>% 
  select(scrubbed_species_binomial, longitude, latitude)

## Get rid of duplicate occurrences
dups <- duplicated(spOcc)
spOcc2 <- spOcc[!dups, ] # Around 12807047 duplicates

spOcc_geo <- spOcc2
coordinates(spOcc_geo) <- c("longitude", "latitude")
crs(spOcc_geo) <- CRS("+init=epsg:4326") # WGS 84

## Rasterize Biomes shapefile
r_ref <- raster(biome_shp)
res(r_ref) <- 100000
r_ref[] <- 1:ncell(r_ref)

## Rasterize the biomes
# https://cran.r-project.org/web/packages/fasterize/vignettes/using-fasterize.html
biomes_ras <- rasterize(biome_shp, r_ref)

## Make sure coordinates have the same projection
CRS.new <- CRS(" +proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
+ellps=WGS84 +towgs84=0,0,")
spOcc_geo<-spTransform(spOcc_geo, CRS.new)

## Create species, grid ID, biome dataframe
spOcc_geo$grid_id <- raster::extract(r_ref,spOcc_geo,factors=TRUE)
spOcc_geo$Biomes <- raster::extract(biomes_ras,spOcc_geo,factors=TRUE)



### Example with the sf package
## Example 
grid_tmp<-st_make_grid(biome_shp, 
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
