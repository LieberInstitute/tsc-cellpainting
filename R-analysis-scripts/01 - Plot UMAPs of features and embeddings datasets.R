install.packages("dplyr","umap","ggplot2","magrittr"))
library(dplyr)
library(umap)
library(ggplot2)
library(magrittr)

c25 <- c(
  "dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise",
  "green1", "yellow4", "yellow3","darkorange4", "brown"
)

#Import data
embeddings = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\Spring analysis\\combined embeddings\\combined embeddings - all oxygen - well - normalized by plate.csv")
embeddings <- embeddings[embeddings$Metadata_Cell_line != 'CT29' & embeddings$Metadata_Cell_line != 'CT30',]
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\CellProfiler analysis\\all-oxygen_triplicate_feature-selected-manually-curated.csv")
features <- features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]

set.seed(120)

#Embeddings UMAPs
#Unfiltered
embeddings_umap_data <- embeddings[,11:650] %>% umap()
embeddings_umap_labels <-embeddings[1:7]
embeddings_umap <- as.data.frame(embeddings_umap_data$layout) %>% cbind(.,embeddings_umap_labels)
colnames(embeddings_umap) <- c("UMAP1","UMAP2","Plate","Well","Cell_Line","Sex","Dx","Age","Oxygen")

#Normoxia
embeddings_umap_data_normoxia <- embeddings %>% filter(Metadata_Oxygen=="Normoxia") %>% .[,11:650] %>% umap()
embeddings_umap_labels_normoxia <-embeddings %>% filter(Metadata_Oxygen=="Normoxia") %>% .[1:7]
embeddings_umap_normoxia <- as.data.frame(embeddings_umap_data_normoxia$layout) %>% cbind(.,embeddings_umap_labels_normoxia)
colnames(embeddings_umap_normoxia) <- c("UMAP1","UMAP2","Plate","Well","Cell_Line","Sex","Dx","Age","Oxygen")

#Hypoxia
embeddings_umap_data_hypoxia <- embeddings %>% filter(Metadata_Oxygen=="Hypoxia") %>% .[,11:650] %>% umap()
embeddings_umap_labels_hypoxia <-embeddings %>% filter(Metadata_Oxygen=="Hypoxia") %>% .[1:7]
embeddings_umap_hypoxia <- as.data.frame(embeddings_umap_data_hypoxia$layout) %>% cbind(.,embeddings_umap_labels_hypoxia)
colnames(embeddings_umap_hypoxia) <- c("UMAP1","UMAP2","Plate","Well","Cell_Line","Sex","Dx","Age","Oxygen")

#Features UMAPs
#Unfiltered
features_umap_data <- features[,10:449] %>% umap()
features_umap_labels <-features[1:6]
features_umap <- as.data.frame(features_umap_data$layout) %>% cbind(.,features_umap_labels)
colnames(features_umap) <- c("UMAP1","UMAP2","Plate","Cell_Line","Sex","Dx","Oxygen","Well")

#Normoxia
features_umap_data_normoxia <- features %>% filter(Metadata_Oxygen=="Normoxia") %>% .[,10:449] %>% umap()
features_umap_labels_normoxia <-features %>% filter(Metadata_Oxygen=="Normoxia") %>% .[1:6]
features_umap_normoxia <- as.data.frame(features_umap_data_normoxia$layout) %>% cbind(.,features_umap_labels_normoxia)
colnames(features_umap_normoxia) <- c("UMAP1","UMAP2","Plate","Cell_Line","Sex","Dx","Oxygen","Well")

#Hypoxia
features_umap_data_hypoxia <- features %>% filter(Metadata_Oxygen=="Hypoxia") %>% .[,10:449] %>% umap()
features_umap_labels_hypoxia <-features %>% filter(Metadata_Oxygen=="Hypoxia") %>% .[1:6]
features_umap_hypoxia <- as.data.frame(features_umap_data_hypoxia$layout) %>% cbind(.,features_umap_labels_hypoxia)
colnames(features_umap_hypoxia) <- c("UMAP1","UMAP2","Plate","Cell_Line","Sex","Dx","Oxygen","Well")

#Plot UMAPs
embeddingsNormoxiaSex <- ggplot(embeddings_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Sex)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("palegreen3","darkgreen")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.91),
        legend.background = element_rect(color='black')
  )

embeddingsHypoxiaSex <- ggplot(embeddings_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Sex)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("palegreen3","darkgreen")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.9),
        legend.background = element_rect(color='black')
  )

embeddingsNormoxiaDx <- ggplot(embeddings_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Dx)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("violet","magenta4")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.87,0.91),
        legend.background = element_rect(color='black')
  )

embeddingsHypoxiaDx <- ggplot(embeddings_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Dx)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("violet","magenta4")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.87,0.9),
        legend.background = element_rect(color='black')
  )

embeddingsLinesNormoxia <- ggplot(embeddings_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Cell_Line)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Cell Line") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

embeddingsLinesHypoxia <- ggplot(embeddings_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Cell_Line)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Cell Line") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

platesEmbeddings <- ggplot(embeddings_umap, aes(x = UMAP1,y = UMAP2, color = Plate)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Plate") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

featuresNormoxiaSex <- ggplot(features_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Sex)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("palegreen3","darkgreen")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.16,0.91),
        legend.background = element_rect(color='black')
  )

featuresHypoxiaSex <- ggplot(features_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Sex)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("palegreen3","darkgreen")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.9),
        legend.background = element_rect(color='black')
  )

featuresNormoxiaDx <- ggplot(features_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Dx)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("violet","magenta4")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.11,0.91),
        legend.background = element_rect(color='black')
  )

featuresHypoxiaDx <- ggplot(features_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Dx)) +
  theme_bw() +
  geom_point(size = 2.5) +
  scale_color_manual(values=c("violet","magenta4")) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.13,0.9),
        legend.background = element_rect(color='black')
  )

featuresLinesNormoxia <- ggplot(features_umap_normoxia, aes(x = UMAP1,y = UMAP2, color = Cell_Line)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Cell Line") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

featuresLinesHypoxia <- ggplot(features_umap_hypoxia, aes(x = UMAP1,y = UMAP2, color = Cell_Line)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Cell Line") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

platesFeatures <- ggplot(features_umap, aes(x = UMAP1,y = UMAP2, color = Plate)) +
  theme_bw() +
  geom_point(size = 2.5) +
  labs(color = "Plate") +
  scale_color_manual(values=c25) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold",hjust = 0.5),
        legend.position = "right",
        legend.background = element_rect(color='black')
  )

#Eport plots
ggsave("normoxia_sexClass.tiff", embeddingsNormoxiaSex, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("hypoxia_sexClass.tiff", embeddingsHypoxiaSex, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("normoxia_dxClass.tiff", embeddingsNormoxiaDx, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("hypoxia_dxClass.tiff", embeddingsHypoxiaDx, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("cellLines_normoxia.tiff", embeddingsLinesNormoxia, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("cellLines_hypoxia.tiff", embeddingsLinesHypoxia, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("platesEmbeddings.tiff", platesEmbeddings, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("normoxia_sexClass_features.tiff", featuresNormoxiaSex, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("hypoxia_sexClass_features.tiff", featuresHypoxiaSex, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("normoxia_dxClass_features.tiff", featuresNormoxiaDx, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("hypoxia_dxClass_features.tiff", featuresHypoxiaDx, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("cellLines_features_normoxia.tiff", featuresLinesNormoxia, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("cellLines_features_hypoxia.tiff", featuresLinesHypoxia, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")
ggsave("plates_features.tiff", platesFeatures, width = 8, height = 6, units = 'in', dpi = 800, compression = "lzw")