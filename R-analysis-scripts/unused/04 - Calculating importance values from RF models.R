install.packages("randomForest","magrittr")
library(randomForest)
library(magrittr)

setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\kfold-CV-importance-values\\test\\")

#Calculate and export importance values for Dx classification
DxImportance <- function(data,folds,file_prefix) {
  
  for (i in 1:nrow(folds)) {
    #Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% folds[i,],] %>% .[,-c(1:3,5:9)]
    train <-data[!(data$Metadata_Cell_line %in% folds[i,]),] %>% .[,-c(1:3,5:9)]
    
    #Random Forest
    set.seed(120)
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Dx,
                                 ntree = 2500)
    
    importanceClass <- importance(classifier_RF)
    filename <- paste0(file_prefix,"_",folds[i,1],"_",folds[i,2],"_",folds[i,3],".csv")
    write.csv(importanceClass,filename,row.names=TRUE)
  }
}

#Calculate and export importance values for Sex classification
SexImportance <- function(data,file_prefix) {
  
  for (i in 1:10) {
    #Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:2,4:9)]
    train <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:2,4:9)]
    
    #Random Forest
    set.seed(120)
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Sex,
                                 ntree = 2500)
    
    importanceClass <- importance(classifier_RF)
    filename <- paste0(file_prefix,"_",foldsSex[i,1],"_",foldsSex[i,2],"_",foldsSex[i,3],".csv")
    write.csv(importanceClass,filename,row.names=TRUE)
  }
}

#Import data
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\CellProfiler analysis\\all-oxygen_triplicate_feature-selected-manually-curated.csv")
features <-features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]
features$Metadata_Sex <- as.factor(features$Metadata_Sex)
features$Metadata_Dx <- as.factor(features$Metadata_Dx)

#Data subsets
featuresNormoxia <- features[features$Metadata_Oxygen == 'Normoxia',]
featuresNormoxiaM <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Male',]
featuresNormoxiaF <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Female',]
featuresHypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]
featuresHypoxiaM <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Male',]
featuresHypoxiaF <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Female',]

#Define folds for cross-validation
foldsDx <- data.frame(
  Scz1 = c("17A","20A","11A","1A","2A","1002.01","1006.02","1013.04","19A","13A"),
  Ctrl1 = c("16B","9B","8B","3B","4B","1009.04","1016.02","15B","9c1","7c6"),
  Scz2 = c(NA,NA,NA,NA,NA,"7A","6A","12A",NA,NA)
)

foldsDxF <- foldsDx[1:5,]
foldsDxM <- foldsDx[6:10,]

foldsSex <- data.frame(
  Male1 = c("19A","9c1","6A","7c6","1016.02","1002.01","7A","1009.04","15B","1013.04"),
  Female1 = c("9B","8B","11A","20A","1A","3B","2A","4B","16B","17A"),
  Male2 = c(NA,NA,NA,NA,"13A","1006.02","12A",NA,NA,NA)
)

#Run RF modeling
DxImportance(featuresNormoxia,foldsDx,"normoxiaDx")
DxImportance(featuresNormoxiaM,foldsDxM,"normoxiaDxM")
DxImportance(featuresNormoxiaF,foldsDxF,"normoxiaDxF")
DxImportance(featuresHypoxia,foldsDx,"hypoxiaDx")
DxImportance(featuresHypoxiaM,foldsDxM,"hypoxiaDxM")
DxImportance(featuresHypoxiaF,foldsDxF,"hypoxiaDxF")
SexImportance(featuresNormoxia,"normoxiaSex")
SexImportance(featuresHypoxia,"hypoxiaSex")
