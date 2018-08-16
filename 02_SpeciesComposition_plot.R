# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(foreach)
library(fuzzySim)
library(circlize)
library(RColorBrewer)
library(wesanderson)
library(dendextend)
library(sf)
library(gplots)
library(ggplot2)

# data --------------------------------------------------------------------
# 1. Presence of species in cells
spOcc<-readRDS("./data/base/2018_02_08_BIEN_OccurData.RData")

# 2. Shapefiles
biome_shp<-st_read("./data/maps/Olson_processed/Biomes_olson_projected.shp")


###############
# PROCEDURE
###############

### define coordinates and convert to SpatialPointsDataFrame
## Working with a subset
occ = spOcc[sample(nrow(spOcc), 1000), ]

spOcc_geo<-
  st_as_sf(occ, coords = c("longitude","latitude"),
           crs = 4326)


spOcc_geo<-st_transform(spOcc_geo,"+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# create 100km grid with grid_id as the name of the biomes


grid_biome <- function(x){
  st_make_grid(biome_shp[biome_shp$biomes==x,], cellsize = c(100000, 100000)) %>%
    st_sf(grid_id = rep(x,length(.)))
}
  
grid_100_biomes<-lapply(biome_shp$biomes,grid_biome) %>% rbindlist() %>% st_sf()

tmp<-st_sf(rbindlist(grid_100_biomes))


grid_100 <- st_make_grid(biome_shp, cellsize = c(100000, 100000)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid_100) %>% cbind(st_coordinates(.))     
                  
               
# view the sampled points, polygons and grid
ggplot() +
  geom_sf(data = biome_shp[biome_shp$biomes=="Moist_Forest",], fill = 'white', lwd = 0.05) +
  geom_sf(data = spOcc_geo[1,], color = 'red', size = 1.7) + 
  geom_sf(data = grid_100, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")


ggplot() +
  geom_sf(data = grid_100, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2)

# which grid square is each point in?
spOcc_geo[1,] %>% st_join(grid_100_biomes, join = st_intersects) %>% as.data.frame


# Ref: https://gis.stackexchange.com/questions/282750/identify-polygon-containing-point-with-r-sf-package
# intersect and extract state name
region <- apply(st_intersects(spOcc_geo[1:100,],biome_shp, sparse = FALSE), 1, 
                     function(col) { 
                       biome_shp[which(col), ]$biomes
                     })



## Starts from here!!


## Extract the biome classification for each grid cell
# Also extract the area cover for each biome in each grid cell
Cells_biomes2<-
  foreach(i=1:length(biome_shp$biomes), .combine = rbind)%do%
  
  {
    print(paste("Extract cells from",biome_shp$biomes[i]))
    biome_tmp<-biome_shp[i,]
    int <- as_tibble(st_intersection(st_buffer(biome_tmp, 0),ms))
    int$areaBiome <- st_area(int$geometry)
    
    int
  }

## Area of a pixel
area_ref<-st_area(p[1,])

# Calculate the proportion of pixel area in each biome
tb_biome <- 
  Cells_biomes %>%
  group_by(cell) %>%
  mutate(areaProp = (areaBiome*100)/area_ref) %>% 
  mutate(maxArea=max(areaProp))

## Select the biome that have the highest proportion of land of the pixel. 
tb_biome <- 
  tb_biome %>% 
  filter(areaProp==maxArea) %>% 
  dplyr::select(biomes,cell)

spPresence_biome<-merge(spPresence, tb_biome, by.x="cells", by.y="cell")
saveRDS(spPresence_biome, file="./outputs/spPresence_biomes_all.rds")

## Number of cells per species in each biome
cells_in_sp<-spPresence_biome %>% 
  group_by(Species,biomes) %>% 
  summarise(N_cells=n_distinct(cells)) %>% 
  group_by(Species) %>% 
  mutate(Total_cells=sum(N_cells), prop_cells=N_cells/sum(N_cells)) %>% 
  mutate(max_prop=max(prop_cells))

saveRDS(cells_in_sp, file="./outputs/spPresence_cell_prop_biomes_all.rds")

cells_in_sp$biomes<-recode(cells_in_sp$biomes,Moist_Forest="Moist",
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

# 2.1 Total numbr of species
Total_sp_list<-tapply(cells_in_sp$Species,cells_in_sp$biomes,unique)

# 2.2 Species with highest proportion of their ranges in each biome
Wides_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==max_prop)

Wides_sp_list<-tapply(Wides_sp$Species,Wides_sp$biomes,unique)

# 2.3 Endemics for each biome
Endemics_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==1)

Endemics_sp_list<-tapply(Endemics_sp$Species,Endemics_sp$biomes,unique)

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

biome_order<-c("Moist","Savannas","Dry",
               "Xeric","Trop_Grass",
               "Coniferous","Temp_Mixed","Temp_Grass",
               "Mediterranean","Taiga","Tundra")

spSimilarity_1<-spSimilarity[biome_order,biome_order]


col=c(wes_palette("Darjeeling1",6,type="continuous"),
      wes_palette("Cavalcanti1",5,type="continuous"))

diag(spSimilarity_1)<-0
colnames(spSimilarity_1)<-c("Moist","Savannas","Dry",
                            "Xeric","Trop_Grass",
                            "Coniferous","Temp_Mixed","Temp_Grass",
                            "Mediterranean","Taiga","Tundra")

indx<-match(colnames(spSimilarity_1),names(prop_endemics))
colnames(spSimilarity_1)<-paste(colnames(spSimilarity_1),", ", prop_endemics[indx],"%", sep="")

rownames(spSimilarity_1)<-colnames(spSimilarity_1)

pdf("./figs/species_composition/Total_similarity_biomes_withEndemics.pdf")
par(mar=c(0, 0, 0, 0))
chordDiagram(spSimilarity_1, annotationTrack = "grid", preAllocateTracks = 1, grid.col =col,symmetric = TRUE,
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

spSimilarity_ma<-spSimilarity_ma[biome_order,biome_order]

## Print file to include into the supplementary information
diag(spSimilarity_ma)<-endemics_n[biome_order]
write.csv(spSimilarity_ma,"./supp_info/Shared_species_matrix.csv")

diag(spSimilarity_ma)<-0
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))

# Rename biomes
colnames(spSimilarity_ma)<-c("Moist","Savannas","Dry",
                             "Xeric","Trop Grass","Coniferous","Temp Mixed",
                             "Temp Grass","Mediterranean","Taiga","Tundra")


## Proportion of widespread species
total_n<-unlist(lapply(Total_sp_list,length))
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))
prop_widespread<-round(Wides_sp_total/total_n,2)*100


colnames(spSimilarity_ma)<-paste(colnames(spSimilarity_ma),", ", prop_widespread[biome_order],"%", sep="")
rownames(spSimilarity_ma)<-colnames(spSimilarity_ma)

pdf("./figs/species_composition/Total_similarity_biomes_DominantSp.pdf",width = 8, height = 8)
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
  color_branches(1,col=wes_palette("Cavalcanti")[2]) %>% 
  set("branches_lwd", 4)  %>% 
  set("labels_cex", 1.5)

pdf("./figs/species_composition/species_composition_cluster_allsp.pdf", height = 9.7, width = 9.6)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

### Heatmaps

my_palette <- colorRampPalette(c(wes_palette("Cavalcanti1")[2],
                                 wes_palette("Cavalcanti1")[c(2,4)],
                                 "white"))(n = 100)

my_palette <-rev(colorRampPalette(c('#ffffcc','#c2e699','#78c679','#31a354','#006837','#006837'))(n = 100))


my_palette <-colorRampPalette(c("#02401b","#02401b","#32806e","white"))(n = 100)


col_breaks<-seq(0,1,by=0.01)

pdf("./figs/species_composition/species_composition_heatmap.pdf", width = 10)
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


## Dissimilarity with dominant species
Dominant_similarity<-Similarity_sp_biomes(Wides_sp_list)
fit_Dominant_similarity <-hclust(as.dist(1-Dominant_similarity))

## Check the order first
#labels(fit_Dominant_similarity)
labels(fit_Dominant_similarity)<-c("Taiga","Tundra","Mediterranean", "Trop grass", 
                                   "Trop Dry", "Xeric","Moist","Savannas", 
                                   "Coniferous", 
                                   "Temp Grass","Temp Mixed")
  
dend_dom<-
  fit_Dominant_similarity %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[3]) %>% 
  set("branches_lwd", 4) %>% 
  set("labels_cex", 1.5)

pdf("./figs/species_composition/species_composition_cluster_Dominant_sp.pdf", height = 10, width = 9.1)
circlize_dendrogram(dend_dom,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()
