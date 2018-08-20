# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)
library(sp)
library(foreach)
library(raster)
library(viridis)
library(ggridges)
library(wesanderson)
library(lattice)
library(ggplot2)
library(ggpubr)


# Settings ----------------------------------------------------------------
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

# Functions ---------------------------------------------------------------
source("./functions/check_functions.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/trait_distribution_functions.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/traits/traits_ALLMB_lambda.csv")

#2. Richness raster and dataframe with cells and biome info
r_Total_Rich <- raster("./data/base/BIEN_2_Ranges/richness100km.tif")
cell_sp_biomes <- readRDS("./outputs/02_Species_grid_id_biomes_df.rds")

#4. Biomes shapefiles
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

# Check raster and shapefiles are in the same projection
#proj4string(r_Total_Rich)==proj4string(biome_shp)

# Presence absence matrix of species in each biome -------------------------------------------------
TraitSpecies <- unique(Traits_phylo$species)

## Filter species occurrence with those that have traits
cell_sp_biomes$Species <- gsub(" ", "_", cell_sp_biomes$scrubbed_species_binomial)

cell_sp_biomes_traits <-
  cell_sp_biomes %>% 
  filter(Species%in%TraitSpecies)

## Presence/absence of species in each biome
Biomes_pabs_cells <- table(cell_sp_biomes_traits$grid_id,cell_sp_biomes_traits$Species)

index <- match(rownames(Biomes_pabs_cells),as.character(cell_sp_biomes_traits$grid_id))

rownames(Biomes_pabs_cells) <- cell_sp_biomes_traits$Biomes[index]

biomes_names<-unique(cell_sp_biomes_traits$Biomes)

# Calculate the relative "abundance" of species in each biome
Biomes_Abun_sp<-foreach(i=1:length(biomes_names),.combine=rbind)%do%
{
  
  indx<-which(row.names(Biomes_pabs_cells)==biomes_names[i])
  biome_abun<-colSums(Biomes_pabs_cells[indx,])
  
}
row.names(Biomes_Abun_sp)<-biomes_names

#3. Final presence/absence matrix species vs biomes
Biomes_pabs_sp<-Biomes_Abun_sp
Biomes_pabs_sp[which(Biomes_pabs_sp>0)]<-1

# Compute distance matrix of trait between each pair of species  ----------
rownames(Traits_phylo)<-Traits_phylo$species

# Select traits to calculate the distances among species
Traits_phylo$logseed_mass<-log(Traits_phylo$Seed_mass)
Traits_phylo$logHeight<-log(Traits_phylo$Height)
Traits_phylo$logWoodDensity<-log(Traits_phylo$Wood_density)
Traits_phylo$sqrtSLA<-sqrt(Traits_phylo$SLA)

traits<-c("logWoodDensity","Leaf_N","sqrtSLA","logseed_mass","logHeight","Leaf_P")

Dist_matrix<-compute_dist_matrix(Traits_phylo[,traits],metric="euclidean",
                                 center = TRUE,
                                 scale = TRUE) ## This can take a while



# Compute functional distinctiveness per biome ----------------------------

# 1. Calculating distinctiveness
# Calculating relative abundance
Biome_relAbun<-make_relative(Biomes_Abun_sp)
Biomes_di_abun = distinctiveness(Biome_relAbun, Dist_matrix)

Biomes_di_abun_clean<-
  Biomes_di_abun%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_di_abun))%>%
  gather(key="Biome",value="Di",-species) %>% 
  filter(!is.na(Di))


# 2. Calculating distinctiveness
Biomes_di = distinctiveness(Biomes_pabs_sp, Dist_matrix)

Biomes_di_clean<-
  Biomes_di%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_di))%>%
  gather(key="Biome",value="Di",-species) %>% 
  filter(!is.na(Di))

indx<-match(Biomes_di_clean$species,Biomes_di_abun_clean$species)
Biomes_di_clean$Di_abun<-Biomes_di_abun_clean$Di[indx]

# Compute functional Uniqueness per biome ----------------------------
Biomes_ui = apply(Biomes_pabs_sp, 1,
               function(x, dist_m) {
                 smaller_com = x[x > 0 & !is.na(x)]
                 uniqueness(t(as.matrix(smaller_com)), dist_m)
               }, dist_m = Dist_matrix)

Biomes_ui_clean<-do.call(rbind.data.frame, Biomes_ui)
Biomes_ui_clean$Biome<-gsub('[0-9]+', '', rownames(Biomes_ui_clean))
Biomes_ui_clean$Biome<-gsub('\\.', '', Biomes_ui_clean$Biome)

# Compute functional restrictiness per biome ----------------------------

## A 0 value indicates that the focal species is present in all the sites.
rest_species<-
  foreach(i=1:length(biomes_names), .combine = rbind)%do%
  {
    indx<-which(rownames(Biomes_pabs_cells)==biomes_names[i])
    biome_PAbs_tmp<-Biomes_pabs_cells[indx,]
    rest_species<-restrictedness(pres_matrix = biome_PAbs_tmp)
    rest_species$Biome<-biomes_names[i]
    
    rest_species
  }

## Exclude all that are not present in any biome cell (which are not in Biome_Di)
# Usually species restricted to just one cell which is not classified in any of the biomes
rest_species<-
  rest_species %>% 
  filter(species%in%unique(Biomes_di_clean$species))


## Scaling Di values per biome
Biomes_di_clean<-
  Biomes_di_clean %>% 
  group_by(Biome) %>% 
  mutate(DiScale=rescaleRas01(Di))

## Scaling Ui values per biome
Biomes_ui_clean<-
  Biomes_ui_clean %>% 
  group_by(Biome) %>% 
  mutate(UiScale=rescaleRas01(Ui))

## Merge Di and Ri
Biome_Di_Ri<-merge(Biomes_di_clean, rest_species)
Biome_Di_Ri<-merge(Biome_Di_Ri,Biomes_ui_clean)

write.csv(Biome_Di_Ri, "./outputs/04_Biome_Di_Ri_phylo.csv")

Biome_Di_Ri$Widespread<-1-Biome_Di_Ri$Ri

## Heatmaps
# Total headmap
Biome_Di_Ri$Biome<-as.factor(Biome_Di_Ri$Biome)

Biome_Di_Ri$Biome<-factor(Biome_Di_Ri$Biome, levels=c("Moist_Forest","Savannas","Tropical_Grasslands",
                                                       "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                                       "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                                       "Taiga","Tundra"))

my_palette <-colorRampPalette(c("#02401b","#02401b","#32806e","white"))(n = 10)



# Hexagonal binning
pdf("./figs/04_All_biomes_heatmap_logTraits2.pdf")
Biome_Di_Ri %>% 
  ggplot(aes(Widespread, DiScale)) +
  stat_binhex(bins=20,aes(fill=log(..count..)))+
  scale_fill_gradientn(colours=rev(c("#02401B","#02401B","#2C7C4C","#2C7C4C","#69A682","#ACCCB9")),name = "log(Richness)")  +
  theme_minimal()+
  facet_wrap( ~ Biome, ncol = 3)+
  ylab("Distinctiveness")+
  xlab("Widespreadness")
dev.off()