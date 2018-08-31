# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(foreach)
library(raster)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/Biomes_hypervolumes_fun.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/traits/traits_ALLMB_lambda.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/04_Biome_Di_Ri_phylo.csv", row.names = 1)

# 5. Presence matrix of species
species_cell_biomes <- readRDS("./outputs/02_Species_grid_id_biomes_df.rds")

## Remove species without traits
species_cell_biomes$Species <- gsub(" ", "_", species_cell_biomes$scrubbed_species_binomial)

species_cell_biomes<- species_cell_biomes %>%
  dplyr::filter(Species%in%unique(Traits_phylo$species))


spMatrix_sub <- table(species_cell_biomes$grid_id,species_cell_biomes$Species)


  # Data manipulation -------------------------------------------------------
# 1. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)


## 3. Rename and order biomes
Traits_Biome_Di_Ri$Biome<-recode(Traits_Biome_Di_Ri$Biome,Moist_Forest="Moist",
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

Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome,
                                 levels=c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra"))

# Calculate hypervolumes -------------------------------

# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density)
Traits_Biome_Di_Ri$sqrtSLA<-sqrt(Traits_Biome_Di_Ri$SLA)

#Selecting and Scalling variables
Traits_Biome_Di_Ri<-
  Traits_Biome_Di_Ri %>%
  mutate(Scaled_logSeed_mass=as.numeric(scale(logseed_mass)),
         Scaled_logHeight=as.numeric((logHeight)),
         Scaled_SLA=as.numeric(scale(sqrtSLA)),
         Scaled_logWood_density=as.numeric(scale(logWoodDensity)),
         Scaled_Leaf_N=as.numeric(scale(Leaf_N)),
         Scaled_Leaf_P=as.numeric(scale(Leaf_P))
  )


## Hypervolumes using a list of species function

Trait_df<-
  Traits_Biome_Di_Ri %>%
  dplyr::select(species,contains("Scaled"))


### Taking only the 10% of the cells per each biomes
cell_biomes <-tapply(species_cell_biomes$grid_id, species_cell_biomes$Biomes, unique)

Random_cells<-
  lapply(cell_biomes,
         function(x)
           sample(x,length(x)*.20)
  )

cells_names<-as.character(as.vector(unlist(Random_cells)))

Tmp<-NULL
count <- 0

system.time(

for (i in cells_names)
{
  print(i)
  count <- count + 1
  x <- spMatrix_sub[i,]

  print(paste("Processing",count, "out of ",length(cells_names)))

  sp_names<-names(x[x > 0 & !is.na(x)])

  #if(length(sp_names)>100){

   # sample_sp<-sample(sp_names,100)
  #}else{
   # sample_sp<-sp_names
  #}

  if (length(sp_names)>1){

    res<- tryCatch({
      cell_hyper<-Trait_df %>%
        filter(species%in%sp_names) %>%
        dplyr::select(contains("Scaled")) %>%
        hypervolume_gaussian()

      cell_hyper@Volume

    },
    error = function(cond){
      message("Species with the same trait values")
      return(NA)
    })

  }else{

    res=NA

  }

  tmp_df <- data.frame(cell = i, vol = res)

  Tmp<-rbind(Tmp,tmp_df)

  write_rds(Tmp, "06_Hypervolume_sp_sample_box_occurrences.rds")
}
)

cell_hyper_df <- Tmp

indx<-match(cell_hyper_df$cell,species_cell_biomes$grid_id)
cell_hyper_df$biomes<-species_cell_biomes$Biomes[indx]

cell_hyper_df$biomes<-factor(cell_hyper_df$biomes,
                             levels=c("Moist_Forest","Savannas","Tropical_Grasslands",
                                      "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                      "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                      "Taiga","Tundra"))

cell_hyper_df$biomes<-recode(cell_hyper_df$biomes,Moist_Forest="Moist",
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

library(wesanderson)

pdf("./figs/06_Hypervolume_cells_sp.pdf", width=10)
ggplot(data=cell_hyper_df,aes(x=biomes,y=vol)) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("SD"^"6")))
dev.off()



## Ignoring the cells that have less than 100 records

cell_N_records <- rowSums(spMatrix_sub)
undesampled_cells <- names(cell_N_records)[which(cell_N_records<100)]
ix <- which(cell_hyper_df$cell%in%undesampled_cells==FALSE)

pdf("./figs/06_Hypervolume_cells_occ_NoUndersampled.pdf", width=10)
ggplot(data=cell_hyper_df[ix,],aes(x=biomes,y=vol)) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("SD"^"6")))
dev.off()

### Richness vs hypervolumes plots

spMatrix_sub_copy <- spMatrix_sub

spMatrix_sub_copy[which(spMatrix_sub_copy>0)]<-1

cell_richness <- rowSums(spMatrix_sub_copy)

indx<-match(cell_hyper_df$cell,names(cell_richness))
cell_hyper_df$Richness<-cell_richness[indx]


## Ignoring the cells with hypervolumes more than 250
cell_hyper_df %>% 
  filter(vol<350) %>% 
  with(plot(log(Richness),vol))


# Use box plot as marginal plots
library(ggpmisc)
library(ggpubr)

biomes_to_plot <- c("Mediterranean","Dry","Xeric")

cell_hyper_df$logRich <- log(cell_hyper_df$Richness)
cell_hyper_df$sqrtVol <- sqrt(cell_hyper_df$vol)

tmp_df <- cell_hyper_df %>% 
  filter(biomes%in%biomes_to_plot)

ggscatterhist(data = tmp_df, x = "logRich", y = "vol",
              color = "biomes", size = 3, alpha = 0.6,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              margin.plot = "boxplot",
              ggtheme = theme_bw())


## Ignoring undersampled grids
ix <- which(tmp_df$cell%in%undesampled_cells==FALSE)
ggscatterhist(data = tmp_df[ix,], x = "logRich", y = "vol",
              color = "biomes", size = 3, alpha = 0.6,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              margin.plot = "boxplot",
              ggtheme = theme_bw())