install.packages(c("tidyverse","magrittr"))
library(tidyverse)
library(magrittr)

#Find all importance .csv files
setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\kfold-CV-importance-values\\test\\")
files <-list.files(pattern = '\\.csv$')

#Define groupings
channels <- c("AGP","DNA","ER","Mito","RNA","Other")
features <- c("Correlation","Granularity","Intensity","RadialDistribution","Texture","AreaShape","Other")

#Group features and plot
for (i in 1:length(files)) {
  importance = read.csv(files[[i]])
  colnames(importance) <- c("Feature","Gini")
  
  new_search <- c("Cells", "Cytoplasm","Nuclei") %>% paste(collapse = "|")
  new <- importance %>% mutate(Compartment = str_extract(paste(Feature), new_search))
  
  new_search2 <- c("DNA", "ER","AGP","RNA","Mito","mito") %>% paste(collapse = "|")
  new2 <- new %>% mutate(Channel = str_extract(paste(Feature), new_search2))
  new2[is.na(new2)] <-"Other"
  
  new_search3 <- c("AreaShape", "Granularity","Intensity","Correlation","Neighbors","RadialDistribution","Texture","ObjectSkeleton","Parent") %>% paste(collapse = "|")
  new3 <- new2 %>% mutate(Description = str_extract(paste(Feature), new_search3)) %>% .[order(.$Gini, decreasing = TRUE),]
  new3$Channel <- gsub("mito","Mito",new3$Channel)
  
  new3$Description <- gsub("Neighbors","Other",new3$Description)
  new3$Description <- gsub("ObjectSkeleton","Other",new3$Description)
  new3$Description <- gsub("Parent","Other",new3$Description)
  
  nuclei <- new3[new3$Compartment == "Nuclei",]
  cells <- new3[new3$Compartment == "Cells",]
  cyto <- new3[new3$Compartment == "Cytoplasm",]
  
  # All compartments
  plot <-
    ggplot(new3) +
    theme_bw() +
    geom_bar(aes(x = factor(Channel,channels), y = Gini, fill = Channel),
             position = "stack",
             stat = "identity",
             color = "black",
             linewidth = .3) +
    ylab("Mean Decrease Gini") +
    scale_fill_manual(values = c("#339933","#0066CC","#FFFF66","#CC0000","grey","#66CCCC")) +
    coord_flip() +
    facet_grid(factor(Description,features) ~ Compartment,
               scales = "free_y",
               space = "free_y") +
    scale_x_discrete(expand = c(0,0)) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(size = 6),
          panel.spacing.y = unit(1,"lines"))
  
  filename <- paste0(gsub(".csv","",files[[i]]),".tiff")
  ggsave(filename,plot,width = 5, height = 7, units = 'in', dpi = 800, compression = "lzw")
}
