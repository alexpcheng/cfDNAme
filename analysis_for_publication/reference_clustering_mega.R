#!/usr/bin/env RScript
# Title: Unsupervised clustering of reference methylomes and sample projection
# Authors: Alexandre Pellan Cheng
# Brief description: Generates figures X and Y for [publication]

# Assume setwd = project_folder/Bin
#Initialize console ----------------------------------------------------------
rm(list=ls())

#libraries --------------------------------------------------------------------
library(data.table)
library(parallel)
library(ggplot2)
library(umap)
library(dplyr)
library(ggpubr)
library(pals)
library(scales)
library(ggtern)
source('~/theme_alex.R')
#Colors ------------------------------------------------------------------------
as.vector(pals::alphabet(19))
custom_palette <- c("#0075dc", "#0007dc", "#00dcd5",
                    "#4C005C", "#191919",
                    "#005C31",
                    "#ffcc99", "#ff4da6", "#ff9999",
                    "#94FFB5", "#8F7C00", "#9DCC00",
                    "#C20088", 
                    "#003380",
                    "#FFA405",
                    "#808080", #"#FFA8BB",
                    "#426600",
                    "#FF0010",
                    "#993F00")



#Functions --------------------------------------------------------------------
# Creates a list with each element containing the reference matrix or a sample.
# NANs are removed from dataframes and sample methylation percentages are 
# scaled to a fraction to be consistent with the reference matrix.

process_files<-function(filename, references.path){
  col<-fread("../references/reference_methylomes/MEGAKA_Refs/reference_methylomes_megaka.tsv", header=TRUE)
  tiled<-fread(paste0(references.path, filename))
  tiled<-tiled[complete.cases(tiled),]
  colnames(tiled)<-c("chr", "start", "end", col$tissue_name)
  tiled <- tiled[, !c('colon_5', 'hsc_5', 'hsc_2')]
  return(tiled)
}

get_UMAP <- function(references.path, random_state, n_neighbors, min_dist){
  references_list<-list.files(path=references.path, pattern="MethylMatrix_binned")
  #Select common regions to references and samples ------------------------------
  refs <- process_files(references_list[[1]], references.path)
  num_references<-ncol(refs)-3
  refs.features<-refs[,4:ncol(refs)]
  ref <- data.frame(t(refs.features))
  
  
  m<- melt(refs.features)
  
  ggplot(data=m)+geom_density(aes(x=value, color=variable))+theme(legend.position="none")
  
  # Perform unsupervised clustering on references -------------------------------
  # Kmeans
  pca <- prcomp(t(refs.features))
  
  # UMAP
  UMAP=umap(t(refs.features), #pca$x, 
            random_state=random_state, #3
            n_neighbors= n_neighbors, #15,
            min_dist= min_dist, #0.1,
            metric = 'euclidean') #random state set for consistency
  UMAP_dims<-data.frame(UMAP$layout)
  UMAP_dims$sample<-factor(colnames(refs.features), levels=colnames(refs.features))
  UMAP_dims$tissue<-gsub("_.*","", UMAP_dims$sample)
  
  return(UMAP_dims)
}

UMAP_dims_golden <- get_UMAP(references.path = "../references/reference_methylomes/MEGAKA_Refs/",
                             random_state = 42,
                             n_neighbors = 15,
                             min_dist=0.1)

UMAP_dims_golden$cluster <- kmeans(UMAP_dims_golden[, c('X1', 'X2')], centers=10)$cluster
UMAP_dims_golden$tissue <- factor(UMAP_dims_golden$tissue, levels=c('macrophage', 'monocyte', 'dendritic',
                                                                    'neutrophil', 'eosinophil', 
                                                                    'tcell', 'nkcell', 'bcell',
                                                                    'erythroblast', 'hsc', 'progenitor',
                                                                    'spleen', 'heart', 'skin','kidney', 'lung',
                                                                    'liver', 
                                                                    'pancreas',
                                                                    'colon', 'megakaryocyte'))

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UMAP_polychrome_mega.pdf",
                width=47.5/25.4, height=47.5/25.4, paper="special", bg="transparent",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=UMAP_dims_golden)+
  #stat_ellipse(geom="polygon",alpha=0.5, aes(x=X1, y=X2, group=cluster), color='black', fill=NA)+
  geom_point(aes(x=X1, y=X2, fill=(tissue)), size=3, pch=21, color='black')+
  scale_fill_manual(values=as.vector(pals::polychrome(30))[c(1:3,9,5:8,4,10,20,12:19,28)])+
  ggrepel::geom_text_repel(aes(x=X1, y=X2, label=factor(tissue)))
  #scale_fill_manual(values = as.vector(custom_palette))+
  #scale_color_gradient2(low='red', mid='purple', high='green', midpoint=5)+
  theme_alex()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t=0, r=0, b=0, l=0)),
        axis.title.x = element_text(margin = margin(t=0, r=0, b=0, l=0)),
        #axis.title = element_blank(),
        plot.margin=grid::unit(c(1,1,1,1), "pt"))
if(save_fig)dev.off()
# Coloring by UMAP coordinates -----
  if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UMAP_polychrome_mega_zoom.pdf",
                  width=22/25.4, height=22/25.4, paper="special", bg="transparent",
                  fonts="Helvetica", colormodel = "cmyk", pointsize=6)
  ggplot(data=UMAP_dims_golden[grepl("kidney|eryth|megaka|heart|skin|spleen|lung", UMAP_dims_golden$tissue), ])+
    #stat_ellipse(geom="polygon",alpha=0.5, aes(x=X1, y=X2, group=cluster), color='black', fill=NA)+
    geom_point(aes(x=X1, y=X2, fill=(tissue)), size=3, pch=21, color='black')+
    scale_fill_manual(values=as.vector(pals::polychrome(30))[c(4,12,13,14,15,16,28)])+
    ggrepel::geom_text_repel(aes(x=X1, y=X2, label=factor(tissue)))
  #scale_fill_manual(values = as.vector(custom_palette))+
  #scale_color_gradient2(low='red', mid='purple', high='green', midpoint=5)+
  theme_alex()+
    xlab("UMAP1")+
    ylab("UMAP2")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          #axis.title = element_blank(),
          plot.margin=grid::unit(c(1,1,1,1), "pt"))
  if(save_fig)dev.off()
  
color_coder_outliers<-function(x,y, x_min, x_max, y_min, y_max, bl, br, tl, tr){
  # Take in hex colors
  bottom_left <- bl
  bottom_right <- br
  top_left <- tl
  top_right <- tr
  
  x_dim_avg <- scales::rescale(c(x_min, x, x_max))[2]
  y_dim_avg <- scales::rescale(c(y_min, y, y_max))[2]
  
  left_right_interpolation_bottom <- colorRamp(colors = c(bottom_left, bottom_right), interpolate = "linear")(x_dim_avg)
  left_right_interpolation_bottom <-
    rgb(left_right_interpolation_bottom[1,1],
        left_right_interpolation_bottom[1,2],
        left_right_interpolation_bottom[1,3],
        maxColorValue = 255)
  
  left_right_interpolation_top <- colorRamp(colors = c(top_left, top_right), interpolate = "linear")(x_dim_avg)
  left_right_interpolation_top <-
    rgb(left_right_interpolation_top[1,1],
        left_right_interpolation_top[1,2],
        left_right_interpolation_top[1,3],
        maxColorValue = 255)
  
  
  final_color <- colorRamp(colors=c(left_right_interpolation_bottom,left_right_interpolation_top), interpolate = "linear")(y_dim_avg)
  
  point_r <- final_color[1,1]
  point_g <- final_color[1,2]
  point_b <- final_color[1,3]
  final_col <- rgb(point_r, point_g, point_b, maxColorValue = 255)
  
  return(final_col)
}





UMAP_dims_golden <- get_UMAP(references.path = "../references/reference_methylomes/processed_refs/",
                             random_state = 3,
                             n_neighbors = 15,
                             min_dist=0.1)
UMAP_dims_golden$hex <- NA
UMAP_dims_golden$n <- NA

for (i in 1:nrow(UMAP_dims_golden)){
  UMAP_dims_golden$n[i]<-as.character(i)
  UMAP_dims_golden$hex[i] <- color_coder_outliers(x = UMAP_dims_golden$X1[i],
                                                  y = UMAP_dims_golden$X2[i],
                                                  x_min = min(UMAP_dims_golden$X1),
                                                  x_max = max(UMAP_dims_golden$X1),
                                                  y_min = min(UMAP_dims_golden$X2),
                                                  y_max = max(UMAP_dims_golden$X2),
                                                  bl = 'darkgrey',
                                                  br = 'orange',
                                                  tl = 'black',
                                                  tr = 'magenta')
}



inter_cluster_distance<-function(df){
  a = max(dist(df[, c('X1', 'X2')]))
  return(a)
}
tt = UMAP_dims_golden
clusters <- c(3)
UMAP_dims_golden=tt

color_count = 1
color_pal = c("white", "orange", "darkblue", "purple", "red", "blue", "green", "yellow")

#color_pal = c("white", "orange", "darkblue", "forestgreen", "red", "yellow", "yellow", "yellow")


for (c in clusters){
  set.seed(4)
  UMAP_dims_golden <- cbind(UMAP_dims_golden, kmeans(UMAP_dims_golden[, c('X1', 'X2')], centers=c, iter.max = 10**5)$cluster)
  index <- ncol(UMAP_dims_golden)
  colnames(UMAP_dims_golden)[[index]]<- 'centers'
  
  
  for (center in seq(c)){
    df <- UMAP_dims_golden[UMAP_dims_golden[, index]==center, ]
    UMAP_dims_golden =  UMAP_dims_golden[!(UMAP_dims_golden[, index]==center), ]
    print(inter_cluster_distance(df))
    print(df)
    if (inter_cluster_distance(df) > 5){
      for (i in 1:nrow(df)){
        df$n[i]<-as.character(i)
        df$hex[i] <- color_coder_outliers(x = df$X1[i],
                                          y = df$X2[i],
                                          x_min = min(df$X1),
                                          x_max = max(df$X1),
                                          y_min = min(df$X2),
                                          y_max = max(df$X2),
                                          bl = color_pal[color_count],
                                          br = color_pal[color_count+1],
                                          tl = color_pal[color_count+2],
                                          tr = color_pal[color_count+3])
        print(color_pal[color_count])
        print(color_pal[color_count+1])
        print(color_pal[color_count+2])
        print(color_pal[color_count+3])
        print('====================')
      }
      color_count= color_count+4
    }
    UMAP_dims_golden = rbind(UMAP_dims_golden, df)
  }
}

UMAP_dims_golden <- UMAP_dims_golden[gtools::mixedorder(UMAP_dims_golden$sample), ]

if(save_fig)pdf(file="../Figures/UMAP_refs_golden_markers_colors_new.pdf",
                width=3, height=3, paper="special", bg="transparent",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=UMAP_dims_golden)+
  #stat_ellipse(geom="polygon",alpha=0.5, aes(x=X1, y=X2, group=centers), color='black', fill=NA)+
  ggrepel::geom_text_repel(aes(x=X1,y=X2, label=sample))+
  geom_point(aes(x=X1, y=X2, fill=(sample)), size=3, pch=21, color='black')+
  scale_fill_manual(values = UMAP_dims_golden$hex)+
  #geom_point(aes(x=X1, y=X2, fill=factor(cluster1)), size=3, pch=21, color='black')+
  #scale_fill_manual(values=pals::kelly(19))+
  theme_alex()+xlab("UMAP1")+ylab("UMAP2")+
  theme(panel.background = element_blank(),
        #panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
#axis.title = element_blank())
if(save_fig)dev.off()

# UMAP_dims_golden_mains <- UMAP_dims_golden[UMAP_dims_golden$X1> -2.5, ]
# 
# for (i in 1:nrow(UMAP_dims_golden_mains)){
#   UMAP_dims_golden_mains$hex[i] <- color_coder_mains(x = UMAP_dims_golden_mains$X1[i],
#                                                    y = UMAP_dims_golden_mains$X2[i],
#                                                    x_min = min(UMAP_dims_golden_mains$X1),
#                                                    x_max = max(UMAP_dims_golden_mains$X1),
#                                                    y_min = min(UMAP_dims_golden_mains$X2),
#                                                    y_max = max(UMAP_dims_golden_mains$X2))
# }
# 
# 
# UMAP_dims_golden$n <- factor(UMAP_dims_golden$n, levels=UMAP_dims_golden$n)
# 
# for (i in 1:nrow(UMAP_dims_golden)){
#   t <- rownames(UMAP_dims_golden)[i]
#   if (t %in% rownames(UMAP_dims_golden_mains)){
#     UMAP_dims_golden$hex[i] <- UMAP_dims_golden_mains[rownames(UMAP_dims_golden_mains)==t, ]$hex
#   }
# }

if(save_fig)pdf(file="../Figures/UMAP_refs_golden_markers_colors_new_zoom.pdf",
                width=1.5, height=1.5, paper="special", bg="transparent",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=UMAP_dims_golden[UMAP_dims_golden$X1>2 & UMAP_dims_golden$X1<5 & UMAP_dims_golden$X2>0 & UMAP_dims_golden$X2<3, ])+
  #ggrepel::geom_text_repel(aes(x=X1,y=X2, label=tissue))+
  geom_point(aes(x=X1, y=X2, fill=(n)), size=3, pch=21, color='black')+
  scale_fill_manual(values = UMAP_dims_golden[UMAP_dims_golden$X1>2 & UMAP_dims_golden$X1<5 & UMAP_dims_golden$X2>0 & UMAP_dims_golden$X2<3, ]$hex)+
  theme_alex()+xlab("UMAP1")+ylab("UMAP2")+
  theme(panel.background = element_blank(),
        #panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
#axis.title = element_blank())
if(save_fig)dev.off()
# Aggregate and save new color scheme for later ----

color_aggregate <- function(coll){
  rgbs <- data.frame(col2rgb(coll))
  means <- rowMeans(rgbs)
  
  hex <- rgb(means[1],
             means[2],
             means[3],
             maxColorValue = 255)
  return(hex)
}

tissues_and_colors <- UMAP_dims_golden[, c("tissue", "hex")]
tissues_and_colors <- aggregate(.~tissue, tissues_and_colors, color_aggregate)

fwrite(file = "tables/tissue_hex_new.txt", x=tissues_and_colors, quote = FALSE, col.names = TRUE, sep = '\t')


medd <- fread('/workdir/apc88/GVHD/Alignment_analysis/sample_output/V1/binned_samples/golden_markers/BRIP01')
ggplot(data=dd)+geom_density(aes(x=dd$V4))




get_PCA <- function(references.path){
  references_list<-list.files(path=references.path, pattern="MethylMatrix_binned")
  #Select common regions to references and samples ------------------------------
  refs <- process_files(references_list[[1]], references.path)
  num_references<-ncol(refs)-3
  refs.features<-refs[,4:ncol(refs)]
  ref <- data.frame(t(refs.features))
  
  m<- melt(refs.features)
  
  # Perform unsupervised clustering on references -------------------------------
  # Kmeans
  pca <- prcomp(t(refs.features))
  pca <- data.frame(pca$x)
  PCA_dims<-data.frame(pca[, c('PC1', 'PC2')])
  PCA_dims$sample<-factor(colnames(refs.features), levels=colnames(refs.features))
  PCA_dims$tissue<-apply(X = PCA_dims, MARGIN = 1, FUN = tissue_names, ncol(PCA_dims))
  PCA_dims$tissue<- factor(PCA_dims$tissue, levels = c("macrophage", "monocyte", "dendritic",
                                                       "eosonophil", "neutrophil",
                                                       "erythroblast",
                                                       "Bcell", "NKell", "Tcell",
                                                       "hema. stem cell/progenitor", "lym. progenitor", "mye. progenitor",
                                                       "spleen",
                                                       "bladder",
                                                       "skin",
                                                       "kidney",
                                                       "liver",
                                                       "pancreas",
                                                       "intestine"))
  return(PCA_dims)
}

PCA_dims_golden <- get_PCA(references.path = "../Reference_Methylomes/MethylMatrix/golden_markers/")

# Coloring by UMAP coordinates -----

PCA_dims_golden$hex <- UMAP_dims_golden$hex
PCA_dims_golden$n <- UMAP_dims_golden$n

if(save_fig)pdf(file="../Figures/PCA_refs_golden_markers_colors_new.pdf",
                width=45/25.4, height=45/25.4, paper="special", bg="transparent",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=PCA_dims_golden)+
  geom_point(aes(x=PC1, y=PC2, fill=(n)), size=3, pch=21, color='black')+
  #geom_text(aes(x=X1,y=X2, label=tissue))+
  scale_fill_manual(values = PCA_dims_golden$hex)+
  theme_alex()+xlab("PC1")+ylab("PC2")+
  theme(panel.background = element_blank(),
        #panel.border = element_blank(),
        panel.grid = element_blank())
#axis.ticks=element_blank(),
#axis.text=element_blank(),
#axis.title = element_blank())
if(save_fig)dev.off()
