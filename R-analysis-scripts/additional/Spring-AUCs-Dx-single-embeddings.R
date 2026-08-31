setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\formal-Dx-reanalysis\\AUCs\\single-channels\\")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------Features----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#Normoxia
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\Spring analysis\\combined embeddings\\combined embeddings - all oxygen - well - normalized by plate.csv")
features <- features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]
features[1:10] <- lapply(features[1:10],as.factor)

features_normoxia <- features[features$Metadata_Oxygen == 'Normoxia',]
featuresNormoxiaF <- features_normoxia[features_normoxia$Metadata_Sex == 'Female',]

trainingDxF <- featuresNormoxiaF[,-c(1:4,6:10)]
trainingDxF_AGP <- featuresNormoxiaF[,-c(1:4,6:10,139:650)]
trainingDxF_DNA <- featuresNormoxiaF[,-c(1:4,6:10,11:138,267:650)]
trainingDxF_ER <- featuresNormoxiaF[,-c(1:4,6:10,11:266,395:650)]
trainingDxF_Mito <- featuresNormoxiaF[,-c(1:4,6:10,11:394,523:650)]
trainingDxF_RNA <- featuresNormoxiaF[,-c(1:4,6:10,11:522)]

trainingDxF_full <- list(trainingDxF, trainingDxF_DNA, trainingDxF_ER, trainingDxF_RNA, trainingDxF_AGP, trainingDxF_Mito)

foldsF <- list()
foldsF[[1]] <- which(featuresNormoxiaF$Metadata_Cell_line != "16B" & featuresNormoxiaF$Metadata_Cell_line != "17A")
foldsF[[2]] <- which(featuresNormoxiaF$Metadata_Cell_line != "9B" & featuresNormoxiaF$Metadata_Cell_line != "20A")
foldsF[[3]] <- which(featuresNormoxiaF$Metadata_Cell_line != "8B" & featuresNormoxiaF$Metadata_Cell_line != "11A")
foldsF[[4]] <- which(featuresNormoxiaF$Metadata_Cell_line != "3B" & featuresNormoxiaF$Metadata_Cell_line != "1A")
foldsF[[5]] <- which(featuresNormoxiaF$Metadata_Cell_line != "4B" & featuresNormoxiaF$Metadata_Cell_line != "2A")

train_controlF <- trainControl(method = "cv",
                               index = foldsF,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

AUC_F <- data.frame(matrix(ncol=6,nrow=5))
colnames(AUC_F) <- c("All","DNA","ER","RNA","AGP","Mito")

for (i in seq_along(trainingDxF_full)) {
  set.seed(120)
  modelDxF <- train(
    Metadata_Dx ~ .,
    data = trainingDxF_full[[i]],
    method = 'ranger',
    trControl = train_controlF,
    metric = "ROC",
  )
  
  AUC_F[,i] <- modelDxF$resample$ROC
}

write.csv(AUC_F,file="caret_Dx_Normoxia_Females_Spring_singleEmbeddings_ranger.csv")

#Hypoxia
features_hypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]
featuresHypoxiaF <- features_hypoxia[features_hypoxia$Metadata_Sex == 'Female',]

trainingDxF <- featuresHypoxiaF[,-c(1:4,6:10)]
trainingDxF_AGP <- featuresHypoxiaF[,-c(1:4,6:10,139:650)]
trainingDxF_DNA <- featuresHypoxiaF[,-c(1:4,6:10,11:138,267:650)]
trainingDxF_ER <- featuresHypoxiaF[,-c(1:4,6:10,11:266,395:650)]
trainingDxF_Mito <- featuresHypoxiaF[,-c(1:4,6:10,11:394,523:650)]
trainingDxF_RNA <- featuresHypoxiaF[,-c(1:4,6:10,11:522)]

trainingDxF_full <- list(trainingDxF, trainingDxF_DNA, trainingDxF_ER, trainingDxF_RNA, trainingDxF_AGP, trainingDxF_Mito)

foldsF <- list()
foldsF[[1]] <- which(featuresHypoxiaF$Metadata_Cell_line != "16B" & featuresHypoxiaF$Metadata_Cell_line != "17A")
foldsF[[2]] <- which(featuresHypoxiaF$Metadata_Cell_line != "9B" & featuresHypoxiaF$Metadata_Cell_line != "20A")
foldsF[[3]] <- which(featuresHypoxiaF$Metadata_Cell_line != "8B" & featuresHypoxiaF$Metadata_Cell_line != "11A")
foldsF[[4]] <- which(featuresHypoxiaF$Metadata_Cell_line != "3B" & featuresHypoxiaF$Metadata_Cell_line != "1A")
foldsF[[5]] <- which(featuresHypoxiaF$Metadata_Cell_line != "4B" & featuresHypoxiaF$Metadata_Cell_line != "2A")

train_controlF <- trainControl(method = "cv",
                               index = foldsF,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

for (i in seq_along(trainingDxF_full)) {
  set.seed(120)
  modelDxF <- train(
    Metadata_Dx ~ .,
    data = trainingDxF_full[[i]],
    method = 'ranger',
    trControl = train_controlF,
    metric = "ROC",
  )
  
  AUC_F[,i] <- modelDxF$resample$ROC
}

write.csv(AUC_F,file="caret_Dx_Hypoxia_Females_Spring_singleEmbeddings_ranger.csv")