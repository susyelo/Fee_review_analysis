# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)
library(factoextra)
library(dendextend)
library(wesanderson)
library(gplots)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/Biomes_hypervolumes_fun.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/traits/traits_ALLMB_lambda.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/04_Biome_Di_Ri_phylo.csv", row.names = 1)


# Data manipulation -------------------------------------------------------
# 1. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)


## 2. Rename and order biomes
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
         Scaled_logHeight=as.numeric(scale(logHeight)),
         Scaled_SLA=as.numeric(scale(sqrtSLA)),
         Scaled_logWood_density=as.numeric(scale(logWoodDensity)),
         Scaled_Leaf_N=as.numeric(scale(Leaf_N)),
         Scaled_Leaf_P=as.numeric(scale(Leaf_P))
  )


# Total species hypervolumes ---------------------------
biome_names<-levels(Traits_Biome_Di_Ri$Biome)

## Hypervolumes for all species
Total_hypervol<-
  Traits_Biome_Di_Ri %>%
  dplyr::select(Biome,contains("Scaled")) %>%
  Biomes_hypervolume(., biome_names)

saveRDS(Total_hypervol, "./outputs/Total_hypervolumes.rds")


pdf("./figs/hypervolumes_clusters/Total_Moist_Temperated_Mixed_Taiga.pdf", width = 12)
plot(
  hypervolume_join(
    Total_hypervol$Moist,
    Total_hypervol$Temp_Mixed,
    Total_hypervol$Taiga
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Temperate mixed", "Taiga"),
       text.col=brewer.pal(n=3,"Set1"),bty="n",cex=1.5,text.font=2)
dev.off()

pdf("./figs/hypervolumes_clusters/Total_Moist_Dry_Savanna_Xeric.pdf", width = 12)
plot(
  hypervolume_join(
    Total_hypervol$Moist,
    Total_hypervol$Dry,
    Total_hypervol$Savannas,
    Total_hypervol$Xeric
  ),
  contour.lwd=1.5,
  colors=c("#FF0000","#91a737","#32806e","#F49c00"),
  cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Dry", "Savannas", "Xeric"),
       text.col=c("#FF0000","#91a737","#32806e","#F49c00"),bty="n",cex=1.5,text.font=2.5)
dev.off()

png("./figs/hypervolumes_clusters/Total_Moist_Dry_Savanna.png", width = 600)
plot(
  hypervolume_join(
    Total_hypervol$Moist,
    Total_hypervol$Dry,
    Total_hypervol$Savanna
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  cex.data=2,cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Dry", "Savanna"),
       text.col=brewer.pal(n=3,"Set1"),bty="n",cex=1.5,text.font=2)
dev.off()


## Similarity hypervolumes with the total species
Total_Sim<-similarity_hypervol(Total_hypervol)
#saveRDS(Total_Sim, "./outputs/05_Total_similarity_hypervolumes.rds")
fit_total <-hclust(as.dist(1-Total_Sim))

dend_total<-
  fit_total %>%
  as.dendrogram() %>%
  color_branches(1,col=wes_palette("Cavalcanti1")[1]) %>%
  set("branches_lwd", 4) %>%
  set("labels_cex", 1.5)

pdf("./figs/hypervolumes_clusters/Total_Sorensen.pdf",height = 11, width = 11.5)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()


my_palette <-rev(colorRampPalette(c('#ffffcc','#c2e699','#78c679','#31a354','#006837','#006837'))(n = 100))


col_breaks<-seq(0,1,by=0.01)

Total_Sim[upper.tri(Total_Sim)]<-t(Total_Sim)[upper.tri(Total_Sim)]
pdf("./figs/05_Total_hypervolumes_similarity_heatmap.pdf", width = 10)
heatmap.2(as.matrix(1-Total_Sim), symm = TRUE,
          distfun = function(x) as.dist(x),dendrogram = "both",margins = c(12,10),
          revC = TRUE,
          cexRow=1.5,cexCol=1.5,
          trace = "none", density.info = "none",keysize = 1.3,
          key.title = "",
          key.xlab = "",
          col=my_palette,
          breaks=col_breaks)
dev.off()


# Redundant and widespread species hypervolumes ---------------------------

## Hypervolumes for widespread and redundant species
Redun_Wides_hypervol<-
  Traits_Biome_Di_Ri %>%
  dplyr::filter(DiScale < 0.25) %>%
  dplyr::select(Biome,contains("Scaled")) %>%
  Biomes_hypervolume(biome_names)

saveRDS(Redun_Wides_hypervol, "./outputs/ReduntWides_hypervolumes_Wides0.5_0.25.rds")

pdf("./figs/hypervolumes_clusters/Redun_Moist_Temperated_Mixed_Taiga.pdf", width = 12)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist,
    Redun_Wides_hypervol$Temp_Mixed,
    Redun_Wides_hypervol$Taiga
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Temperate mixed", "Taiga"),
       text.col=brewer.pal(n=3,"Set1"),bty="n",cex=1.5,text.font=2)
dev.off()

pdf("./figs/hypervolumes_clusters/Redun_Moist_Dry_Savanna_Xeric.pdf", width = 12)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist,
    Redun_Wides_hypervol$Dry,
    Redun_Wides_hypervol$Savannas,
    Redun_Wides_hypervol$Xeric
  ),
  contour.lwd=1.5,
  colors=c("#FF0000","#91a737","#32806e","#F49c00"),
  cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Dry", "Savannas", "Xeric"),
       text.col=c("#FF0000","#91a737","#32806e","#F49c00"),bty="n",cex=1.5,text.font=2.5)
dev.off()


png("./figs/hypervolumes_clusters/Redundant_Moist_Dry_Savanna.png", width = 600)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist,
    Redun_Wides_hypervol$Dry,
    Redun_Wides_hypervol$Savanna
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  cex.data=2,cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Moist","Dry", "Savanna"),
       text.col=brewer.pal(n=3,"Set1"),bty="n",cex=1.5,text.font=2)
dev.off()

## Calculate Hypervolume similarity using Sorense's index -----
## With Redundant species
redun_Sim<-similarity_hypervol(Redun_Wides_hypervol)
#saveRDS(redun_Sim, "./outputs/Redunt_similarity_hypervolumes.rds")
fit_red <-hclust(as.dist(1-redun_Sim))

dend_red<-
  fit_red %>%
  as.dendrogram() %>%
  color_branches(1,col=wes_palette("Cavalcanti1")[3]) %>%
  set("branches_lwd", 4) %>%
  set("labels_cex", 1.5)


#dir.create("./figs/hypervolumes_clusters")
pdf("./figs/hypervolumes_clusters/Redundant_Sorensen0.5_0.25.pdf" ,height = 11, width = 11.5)
circlize_dendrogram(dend_red,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

redun_Sim[upper.tri(redun_Sim)]<-t(redun_Sim)[upper.tri(redun_Sim)]
pdf("./figs/hypervolumes_clusters/Redundant_hypervolumes_similarity_heatmap.pdf", width = 10)
heatmap.2(as.matrix(1-redun_Sim), symm = TRUE,
          distfun = function(x) as.dist(x),dendrogram = "both",margins = c(12,10),
          revC = TRUE,
          cexRow=1.5,cexCol=1.5,
          trace = "none", density.info = "none",keysize = 1.3,
          key.title = "",
          key.xlab = "",
          col=my_palette,
          breaks=col_breaks)
dev.off()


# Hypervolumes for climatic categories ------------------------------------
tropical<-c("Moist","Dry","Trop_Grass","Savannas")
temperate<-c("Temp_Mixed","Coniferous","Temp_Grass","Mediterranean")
cold<-c("Taiga","Tundra")

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%tropical,"Tropical",NA)
Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%temperate,"Temperate",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%cold,"Cold",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes[which(Traits_Biome_Di_Ri$Biome=="Xeric")]<-"Xeric"


## The function needs a variable called "Biome"
## so I created a temporal dataframe with the Biome variables using the climate classification
Traits_Biome_Di_Ri_tmp<-Traits_Biome_Di_Ri
Traits_Biome_Di_Ri_tmp$Biome<-Traits_Biome_Di_Ri_tmp$ClimBiomes

biome_names<-unique(Traits_Biome_Di_Ri_tmp$Biome)


Climatic_hypervol<-
  Traits_Biome_Di_Ri_tmp %>%
  dplyr::select(Biome,contains("Scaled")) %>%
  Biomes_hypervolume(biome_names)

saveRDS(Climatic_hypervol,"./outputs/Climatic_hypervolumes.rds")


Climatic_Sim<-similarity_hypervol(Climatic_hypervol)
fit_Climatic <-hclust(as.dist(1-Climatic_Sim))

dend_climatic<-
  fit_Climatic %>%
  as.dendrogram() %>%
  color_branches(1,col=wes_palette("Cavalcanti1")[1]) %>%
  set("branches_lwd", 4)



pdf("./figs/hypervolumes_clusters/Climatic_hypervolumes_total.pdf",width = 12)
plot(
  hypervolume_join(
    Climatic_hypervol$Tropical,
    Climatic_hypervol$Xeric,
    Climatic_hypervol$Temperate,
    Climatic_hypervol$Cold
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  cex.axis=1,cex.names=1.5,
  show.legend=FALSE,
  names=c("log(seed mass)","log(Height)", "sqrt(SLA)", "log(wood density)", "Leaf N","Leaf P")
)
legend("bottomleft",legend = c("Tropical","Xeric", "Temperate", "Polar"),
       text.col=brewer.pal(n=4,"Set1"),bty="n",cex=1.5,text.font=2)
dev.off()

### Hypervolumes for Dominant growth forms -----

# Extract dominant growth forms per biome ---------------------------------
Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome,
                                 levels=rev(c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra")))

biome_name<-unique(Traits_Biome_Di_Ri$Biome)

Total_GF<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Traits_Biome_Di_Ri %>%
      filter(Biome==biome_name[i]) %>%
      group_by(GROWTHFORM_STD) %>%
      dplyr::summarise(N_sp=length(species)) %>%
      mutate(Dist="Total_prop",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,1))
  }

#Using the ScaleUi values produced the same results as the ScaleDi
Traits_Biome_Di_Ri$DiScale<-round(Traits_Biome_Di_Ri$DiScale,1)
Traits_Biome_Di_Ri$Widespread<-round(Traits_Biome_Di_Ri$Widespread,1)

RedWides_GF<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Traits_Biome_Di_Ri %>%
      filter(Biome==biome_name[i]) %>%
      filter(DiScale < 0.25 & Widespread > 0.5) %>%
      group_by(GROWTHFORM_STD) %>%
      dplyr::summarise(N_sp=length(species)) %>%
      mutate(Dist="Redun_wides",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,1))
  }

Total_GF$Tmnt<-"Total"
RedWides_GF$Tmnt<-"RedWid"

new_df<-rbind(Total_GF,RedWides_GF)
col_GF<-c(wes_palette("Cavalcanti1")[c(2:4,1)],"grey")

pdf("./figs/Growth_forms/Total_vs_redundant_species_Wid50_Dist25.pdf", width = 12, height = 8)
ggplot(data = new_df,
       mapping = aes(x = Biome, fill = GROWTHFORM_STD,
                     y = ifelse(test = Tmnt == "Total",
                                yes = -prop, no = prop))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values=col_GF,name = "Growth form")+
  labs(y = "%",x="") +
  coord_flip()+
  geom_hline(yintercept=0)+
  theme(axis.text.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)))
dev.off()


## Hypervolumes for dominant forms
Traits_Biome_Di_Ri<-as.data.frame(Traits_Biome_Di_Ri)

Dom_growth_forms<-list(Moist=c("Tree"),
                       Savannas=c("Tree","Grass"),
                       Trop_Grass=c("Grass"),
                       Dry=c("Tree"),
                       Xeric=c("Shrub","Tree"),
                       Temp_Mixed=c("Tree"),
                       Coniferous=c("Tree"),
                       Temp_Grass=c("Grass"),
                       Mediterranean=c("Shrub","Tree"),
                       Taiga=c("Tree"),
                       Tundra=c("Herb","Grass"))


Dom_GF_hypervolumes<-
  foreach(i=1:length(names(Dom_growth_forms)))%do%{

    print(names(Dom_growth_forms)[i])
    Traits_Biome_Di_Ri %>%
      filter(Biome==names(Dom_growth_forms)[i] & GROWTHFORM_STD%in%Dom_growth_forms[[i]]) %>%
      dplyr::select(contains("Scaled")) %>%
      hypervolume_box(name = names(Dom_growth_forms)[i])

  }

names(Dom_GF_hypervolumes)<-names(Dom_growth_forms)
saveRDS(Dom_GF_hypervolumes, "./outputs/Dom_growth_forms.rds")


plot(
  hypervolume_join(
    Dom_GF_hypervolumes$Temp_Grass,
    Dom_GF_hypervolumes$Trop_Grass,
    Dom_GF_hypervolumes$Tundra
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE,
  show.3d=TRUE
)

plot(
  hypervolume_join(
    Dom_GF_hypervolumes$Moist,
    Dom_GF_hypervolumes$Coniferous,
    Dom_GF_hypervolumes$Temp_Mixed,
    Dom_GF_hypervolumes$Taiga
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE
)

plot(
  hypervolume_join(
    Dom_GF_hypervolumes$Moist,
    Dom_GF_hypervolumes$Coniferous,
    Dom_GF_hypervolumes$Temp_Mixed,
    Dom_GF_hypervolumes$Taiga
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE,
  show.3d=TRUE
)


Dom_GF_Sim<-similarity_hypervol(Dom_GF_hypervolumes)
#saveRDS(redun_Sim, "./outputs/Redunt_similarity_hypervolumes.rds")
fit_DomGF <-hclust(as.dist(1-Dom_GF_Sim))


dend_DomGF<-
  fit_DomGF %>%
  as.dendrogram() %>%
  color_branches(1,col=wes_palette("Cavalcanti1")[3]) %>%
  set("branches_lwd", 4) %>%
  set("labels_cex", 1.5)

circlize_dendrogram(dend_DomGF,dend_track_height = 0.7,labels_track_height = 0.2)


## Hypervolumes per growth forms
biomes_names<-levels(Traits_Biome_Di_Ri$Biome)

Tree_hypervolumes<-
  foreach(i=1:length(biomes_names))%do%{

    print(biomes_names[i])
    Traits_Biome_Di_Ri %>%
      filter(Biome==biomes_names[i] & GROWTHFORM_STD=="Tree") %>%
      dplyr::select(contains("Scaled")) %>%
      hypervolume_box(name = biomes_names[i])

  }

names(Tree_hypervolumes)<-biomes_names

Shrub_hypervolumes<-
  foreach(i=1:length(biomes_names))%do%{

    print(biomes_names[i])
    Traits_Biome_Di_Ri %>%
      filter(Biome==biomes_names[i] & GROWTHFORM_STD=="Shrub") %>%
      dplyr::select(contains("Scaled")) %>%
      hypervolume_box(name = biomes_names[i])

  }
names(Shrub_hypervolumes)<-biomes_names


Herb_hypervolumes<-
  foreach(i=1:length(biomes_names))%do%{

    print(biomes_names[i])
    Traits_Biome_Di_Ri %>%
      filter(Biome==biomes_names[i] & GROWTHFORM_STD=="Herb") %>%
      dplyr::select(contains("Scaled"),-Scaled_logWood_density) %>%
      hypervolume_box(name = biomes_names[i])

  }
names(Herb_hypervolumes)<-biomes_names


Grass_hypervolumes<-
  foreach(i=1:length(biomes_names))%do%{

    print(biomes_names[i])
    Traits_Biome_Di_Ri %>%
      filter(Biome==biomes_names[i] & GROWTHFORM_STD=="Grass") %>%
      dplyr::select(contains("Scaled"),-Scaled_logWood_density) %>%
      hypervolume_box(name = biomes_names[i])

  }
names(Grass_hypervolumes)<-biomes_names

plot(
  hypervolume_join(
    Tree_hypervolumes$Moist,
    Tree_hypervolumes$Coniferous,
    Tree_hypervolumes$Temp_Mixed,
    Tree_hypervolumes$Taiga
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE
)

plot(
  hypervolume_join(
    Shrub_hypervolumes$Xeric,
    Shrub_hypervolumes$Mediterranean,
    Shrub_hypervolumes$Tundra
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE
)

plot(
  hypervolume_join(
    Shrub_hypervolumes$Dry,
    Shrub_hypervolumes$Xeric,
    Shrub_hypervolumes$Mediterranean
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE
)


plot(
  hypervolume_join(
    Grass_hypervolumes$Trop_Grass,
    Grass_hypervolumes$Temp_Grass,
    Grass_hypervolumes$Tundra
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  show.random=FALSE
)


biomes<-c("Moist","Temp_Mixed")
