library(rgdal)
library(raster)
library(tidyverse)

# 1. data ------------------------------------------------------------------------
# 1.1 Presence of species in cells and biomes with occurrences without filtering
spOcc_geo<-readRDS("outputs/02_Species_grid_id_biomes_df.rds")

# 1.2. Trait species dataframe
Traits_phylo<-read.csv("./data/traits/traits_ALLMB_lambda.csv")

############################################################################
# Occurrence data patterns ------------------------------------------------
############################################################################

# 2. Create total richness dataframe -----------------------------------------
# Convert into a normal dataframe
spOcc_geo_df <- as_tibble(spOcc_geo)

## Species richnnes per biomes
spMatrix_occ <- table(spOcc_geo_df$grid_id,spOcc_geo_df$scrubbed_species_binomial)

spMatrix_occ[which(spMatrix_occ > 0)] <- 1
Sp_Rich_Occ <- rowSums(spMatrix_occ)

## Include the richnnes values to a dataframe
grid_richness_total <- data.frame(grid_id = names(Sp_Rich_Occ), Richness = Sp_Rich_Occ)
indx<-match(grid_richness_total$grid_id,spOcc_geo_df$grid_id)
grid_richness_total$Biomes <- spOcc_geo_df$Biomes[indx]
  
library(wesanderson)

grid_richness_total$Biomes<-recode(grid_richness_total$Biomes,Moist_Forest="Moist",
                                 Savannas="Savannas",
                                 Tropical_Grasslands="Trop_Grass",
                                 Dry_Forest="Dry",
                                 Xeric_Woodlands="Xeric",
                                 Mediterranean_Woodlands="Mediterranean",
                                 Temperate_Grasslands="Temp_Grass",
                                 Temperate_Mixed="Temp_Mixed",
                                 Coniferous_Forests="Coniferous",
                                 Taiga="Taiga",
                                 Tundra="Tundra")

grid_richness_total$Biomes<-factor(grid_richness_total$Biomes,
                                 levels=c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra"))


png("./figs/Total_richnnes_occ.png", width = 700)
ggplot(data=grid_richness_total,aes(x=Biomes,y=log(Richness+1))) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("Richness")))
dev.off()



# 3. Richness for grids with trait values ------------------------------------
spOcc_geo_df$Species <- gsub(" ", "_", spOcc_geo_df$scrubbed_species_binomial)

spOcc_geo_traits <-
  spOcc_geo_df %>% 
  dplyr::filter(Species%in%unique(Traits_phylo$species))

## Species richnnes per biomes
spMatrix_occ_trait <- table(spOcc_geo_traits$grid_id,spOcc_geo_traits$Species)

spMatrix_occ_trait[which(spMatrix_occ_trait > 0)] <- 1
Sp_Rich_Occ_trait<- rowSums(spMatrix_occ_trait)

## Include the richnnes values to a dataframe
grid_richness_traits <- data.frame(grid_id = names(Sp_Rich_Occ_trait), Richness = Sp_Rich_Occ_trait)
indx<-match(grid_richness_traits$grid_id,spOcc_geo_traits$grid_id)
grid_richness_traits$Biomes <- spOcc_geo_traits$Biomes[indx]

grid_richness_traits$Biomes<-recode(grid_richness_traits$Biomes,Moist_Forest="Moist",
                                   Savannas="Savannas",
                                   Tropical_Grasslands="Trop_Grass",
                                   Dry_Forest="Dry",
                                   Xeric_Woodlands="Xeric",
                                   Mediterranean_Woodlands="Mediterranean",
                                   Temperate_Grasslands="Temp_Grass",
                                   Temperate_Mixed="Temp_Mixed",
                                   Coniferous_Forests="Coniferous",
                                   Taiga="Taiga",
                                   Tundra="Tundra")

grid_richness_traits$Biomes<-factor(grid_richness_traits$Biomes,
                                   levels=c("Moist","Savannas","Trop_Grass",
                                            "Dry","Xeric","Mediterranean",
                                            "Temp_Grass","Temp_Mixed","Coniferous",
                                            "Taiga","Tundra"))


png("./figs/Trait_richnnes_occ.png", width = 700)
ggplot(data=grid_richness_traits,aes(x=Biomes,y=log(Richness+1))) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("Richness")))
dev.off()

############################################################################
# Spatial models distributions ---------------------------------------------
############################################################################

# 1. Data  -----------------------------------------------------------------

# Shapefiles
biome_shp <- shapefile("./data/maps/Olson_processed/Biomes_olson_projected.shp")

# Using Richness raster as template raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

# 2. Create dataframes  -----------------------------------------------------
## Rasterize Biomes shapefile 100 km^2
r_ref <- r_Total_Rich
r_ref[] <- 1:ncell(r_ref)

## Rasterize the biomes
# https://cran.r-project.org/web/packages/fasterize/vignettes/using-fasterize.html
biomes_ras <- rasterize(biome_shp, r_ref)

boxplot(log(r_Total_Rich),biomes_ras)


