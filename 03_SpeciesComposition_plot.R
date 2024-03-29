# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(fuzzySim)
library(circlize)
library(RColorBrewer)
library(wesanderson)
library(dendextend)
library(sf)
library(gplots)
library(ggplot2)
library(data.table)
library(foreach)

# data --------------------------------------------------------------------
# 1. Presence of species in cells and biomes
spOcc_geo<-readRDS("outputs/02_Species_grid_id_biomes_df.rds")

# Convert into a normal dataframe
spOcc_geo_df <- as_tibble(spOcc_geo)


## Number of cells per species in each biome
cells_in_sp<-spOcc_geo_df %>% 
  group_by(scrubbed_species_binomial,Biomes) %>% 
  summarise(N_cells=n_distinct(grid_id)) %>% 
  group_by(scrubbed_species_binomial) %>% 
  mutate(Total_cells=sum(N_cells), prop_cells=N_cells/sum(N_cells)) %>% 
  mutate(max_prop=max(prop_cells))

cells_in_sp$Biomes<-recode(cells_in_sp$Biomes,Moist_Forest="Moist",
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

# 2. Species list for each biome ------------------------------------------
cells_in_sp$Species <- cells_in_sp$scrubbed_species_binomial

# 2.1 Total numbr of species
Total_sp_list<-tapply(cells_in_sp$Species,cells_in_sp$Biomes,unique)

# 2.2 Species with highest proportion of their ranges in each biome
Wides_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==max_prop)

Wides_sp_list<-tapply(Wides_sp$Species,Wides_sp$Biomes,unique)

# 2.3 Endemics for each biome
Endemics_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==1)

Endemics_sp_list<-tapply(Endemics_sp$Species,Endemics_sp$Biomes,unique)

# 2.4 Proportion of endemics in each biome
total_n<-unlist(lapply(Total_sp_list,length))
endemics_n<-unlist(lapply(Endemics_sp_list,length))
prop_endemics<-round(endemics_n/total_n,3)*100


# 3. Create similarity matrix ---------------------------------------------
## Create a loop to calculate the similarity (number of species shared among biomes)
biome_richness<-Total_sp_list

spSimilarity<-foreach(i=1:length(biome_richness), .combine='cbind') %:%
  foreach(j=1:length(biome_richness), .combine='c') %do% {
    length(intersect(biome_richness[[i]],biome_richness[[j]]))
  }

colnames(spSimilarity)<-names(biome_richness)
rownames(spSimilarity)<-names(biome_richness)

## Double check that the numbers are correct
#diag(spSimilarity)==unlist(lapply(biome_richness, n_distinct))
#spSimilarity/diag(spSimilarity)

# 4. Chordplot of similarities --------------------------------------------

# 4.1 Species composition among biomes using all the species

col=c(wes_palette("Darjeeling1",6,type="continuous"),
      wes_palette("Cavalcanti1",5,type="continuous"))

diag(spSimilarity)<-0

indx<-match(colnames(spSimilarity),names(prop_endemics))
colnames(spSimilarity)<-paste(colnames(spSimilarity),", ", prop_endemics[indx],"%", sep="")

rownames(spSimilarity)<-colnames(spSimilarity)

pdf("./figs/species_composition/Total_similarity_biomes_withEndemics.pdf")
par(mar=c(0, 0, 0, 0))
chordDiagram(spSimilarity, annotationTrack = "grid", preAllocateTracks = 1, grid.col =col,symmetric = TRUE,
             column.col = col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .2, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.7)
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


# 4.2 Species composition using Dominant species (those that occupy most of their ranges in each biome)
# This is to give some direction of shared species among biome 
spSimilarity_Wides<-foreach(i=1:length(Wides_sp_list), .combine='rbind') %:%
  foreach(j=1:length(Total_sp_list), .combine='rbind') %do% {
    
    N_sp=length(intersect(Wides_sp_list[[i]],Total_sp_list[[j]]))
    df<-data.frame(from=names(Wides_sp_list)[i], 
                   to=names(Total_sp_list)[j],
                   Sp_shared=N_sp)
    df
  }

spSimilarity_ma<-matrix(data = spSimilarity_Wides$Sp_shared, nrow = n_distinct(spSimilarity_Wides$from), ncol = n_distinct(spSimilarity_Wides$to), byrow = FALSE,
            dimnames = NULL)

rownames(spSimilarity_ma)<-unique(spSimilarity_Wides$from)
colnames(spSimilarity_ma)<-unique(spSimilarity_Wides$to)

biome_order<-c("Moist","Savannas","Dry",
               "Xeric","Trop_Grass",
               "Coniferous","Temp_Mixed","Temp_Grass",
               "Mediterranean","Taiga","Tundra")

spSimilarity_ma<-spSimilarity_ma[biome_order,biome_order]

## Print file to include into the supplementary information
diag(spSimilarity_ma)<-endemics_n[biome_order]
write.csv(spSimilarity_ma,"./supp_info/Shared_species_matrix.csv")

diag(spSimilarity_ma)<-0
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))


## Proportion of widespread species
total_n<-unlist(lapply(Total_sp_list,length))
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))
prop_widespread<-round(Wides_sp_total/total_n,2)*100


colnames(spSimilarity_ma)<-paste(colnames(spSimilarity_ma),", ", prop_widespread[biome_order],"%", sep="")
rownames(spSimilarity_ma)<-colnames(spSimilarity_ma)

pdf("./figs/03_Total_similarity_biomes_DominantSp_Occ.pdf",width = 8, height = 8)
par(mar=c(0, 0, 0, 0))
chordDiagram(spSimilarity_ma,column.col = col,
             grid.col =col, directional = -1, 
             direction.type = c("diffHeight", "arrows"),link.largest.ontop=TRUE,
             annotationTrack = c("grid"),link.arr.length = 0.2,link.arr.type = "big.arrow",
             preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.8)
}, bg.border = NA)
dev.off()


# 5. Species dissimilarity among biomes --------------------------------------
Similarity_sp_biomes<-function(sp_list){
  
  biome_similarity <- foreach(i = 1:length(sp_list), .combine='rbind') %:%
    foreach(j = 1:length(Total_sp_list), .combine='rbind') %do% {
      
      sp_intersection = length(intersect(sp_list[[i]],sp_list[[j]]))
      biome1 <- length(sp_list[[i]])
      biome2 <- length(sp_list[[j]])
      
      sorensen <- 2 * sp_intersection/(biome1 + biome2)
      sorensen_df<-data.frame(sorensen=sorensen, 
                             from=names(sp_list)[i], 
                             to=names(sp_list)[j])
      sorensen_df
    }
  
  Similarity_ma<-matrix(data = biome_similarity$sorensen, 
                        nrow = n_distinct(biome_similarity$from), 
                        ncol = n_distinct(biome_similarity$to), byrow = FALSE,
                          dimnames = NULL)
  
  rownames(Similarity_ma)<-unique(biome_similarity$from)
  colnames(Similarity_ma)<-unique(biome_similarity$to)
  
  Similarity_ma
}
  
## Dissimilarity all species  
Total_similarity<-Similarity_sp_biomes(Total_sp_list)

fit_total_sim <-hclust(as.dist(1-Total_similarity))

labels(fit_total_sim)<-c("Trop_Grass", "Moist","Savannas", "Dry","Xeric", 
                         "Taiga","Tundra","Mediterranean","Coniferous","Temp_Grass","Temp_Mixed")

dend_total<-
  fit_total_sim %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti1")[2]) %>% 
  dendextend::set("branches_lwd", 4)  %>% 
  dendextend::set("labels_cex", 1.5)

circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)

### Heatmaps

my_palette <- colorRampPalette(c(wes_palette("Cavalcanti1")[2],
                                 wes_palette("Cavalcanti1")[c(2,4)],
                                 "white"))(n = 100)

my_palette <-rev(colorRampPalette(c('#ffffcc','#c2e699','#78c679','#31a354','#006837','#006837'))(n = 100))


my_palette <-rev(colorRampPalette(c('#ffffcc','#c2e699','#78c679','#31a354','#006837','#006837'))(n = 100))

col_breaks<-seq(0,1,by=0.01)

pdf("./figs/03_species_composition_heatmap.pdf", width = 10)
heatmap.2(as.matrix(1-Total_similarity), symm = TRUE,
          distfun = function(x) as.dist(x),dendrogram = "both",margins = c(12,10),
          revC = TRUE,
          cexRow=1.5,cexCol=1.5,
          trace = "none", density.info = "none",keysize = 1.3,
          key.title = "",
          key.xlab = "",
          col=my_palette,
          breaks=col_breaks)
dev.off()
