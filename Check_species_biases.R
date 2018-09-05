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

cell_N_records <- rowSums(spMatrix_occ)
undesampled_cells <- names(cell_N_records)[which(cell_N_records<100)]


ix <- which(rownames(spMatrix_occ)%in%undesampled_cells==FALSE)
spMatrix_occ <- spMatrix_occ[ix,]

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

## Since there are a lot of NAs that are on the coast I will use the dataframe used in the original analysis

spPresence_biomes <- readRDS("./data/spPresence_biomes_all.rds")

spMatrix_SDM_total <- table(spPresence_biomes$cells,spPresence_biomes$Species)

spMatrix_SDM_total[which(spMatrix_SDM_total > 0)] <- 1
Sp_Rich_SDM_total<- rowSums(spMatrix_SDM_total)

## Include the richnnes values to a dataframe
grid_richness_SDM <- data.frame(grid_id = names(Sp_Rich_SDM_total), 
                                Richness = Sp_Rich_SDM_total)

indx<-match(grid_richness_SDM$grid_id,spPresence_biomes$cells)

grid_richness_SDM$Biomes <- spPresence_biomes$biomes[indx]


grid_richness_SDM$Biomes<-recode(grid_richness_SDM$Biomes,Moist_Forest="Moist",
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

grid_richness_SDM$Biomes<-factor(grid_richness_SDM$Biomes,
                                  levels=c("Moist","Savannas","Trop_Grass",
                                           "Dry","Xeric","Mediterranean",
                                           "Temp_Grass","Temp_Mixed","Coniferous",
                                           "Taiga","Tundra"))

png("./figs/Total_richnnes_SDM.png", width = 700)
ggplot(data=grid_richness_SDM,aes(x=Biomes,y=log(Richness+1))) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("Richness")))
dev.off()

### Richness counting only those species with traits

spPresence_biomes_traits <- 
  spPresence_biomes %>% 
  dplyr::filter(Species%in%unique(Traits_phylo$species))


spMatrix_SDM_trait <- table(spPresence_biomes_traits$cells,spPresence_biomes_traits$Species)

spMatrix_SDM_trait[which(spMatrix_SDM_trait > 0)] <- 1
Sp_Rich_SDM_traits<- rowSums(spMatrix_SDM_trait)

## Include the richnnes values to a dataframe
grid_richness_SDM_trait <- data.frame(grid_id = names(Sp_Rich_SDM_traits), 
                                Richness = Sp_Rich_SDM_traits)

indx<-match(grid_richness_SDM_trait$grid_id,spPresence_biomes_traits$cells)

grid_richness_SDM_trait$Biomes <- spPresence_biomes_traits$biomes[indx]


grid_richness_SDM_trait$Biomes<-recode(grid_richness_SDM_trait$Biomes,Moist_Forest="Moist",
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

grid_richness_SDM_trait$Biomes<-factor(grid_richness_SDM_trait$Biomes,
                                 levels=c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra"))

png("./figs/Trait_richnnes_SDM.png", width = 700)
ggplot(data=grid_richness_SDM_trait,aes(x=Biomes,y=log(Richness+1))) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("Richness")))
dev.off()
