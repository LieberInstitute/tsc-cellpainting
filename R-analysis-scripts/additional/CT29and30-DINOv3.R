setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\formal-Dx-reanalysis\\AUCs\\")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------Features----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#Normoxia
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\dinov3_artifacts\\normalized_features\\libd_dinov3_plate_robustize_USE_THIS.csv")
features[1:10] <- lapply(features[1:10],as.factor)

features_normoxia <- features[features$Metadata_Oxygen == 'Normoxia',]

models <- c("ranger","glmnet","mlp")

folds <- list()
folds[[1]] <- which(features_normoxia$Metadata_Cell_line != "19A" & features_normoxia$Metadata_Cell_line != "9B")
folds[[2]] <- which(features_normoxia$Metadata_Cell_line != "9c1" & features_normoxia$Metadata_Cell_line != "8B")
folds[[3]] <- which(features_normoxia$Metadata_Cell_line != "6A" & features_normoxia$Metadata_Cell_line != "11A")
folds[[4]] <- which(features_normoxia$Metadata_Cell_line != "7c6" & features_normoxia$Metadata_Cell_line != "20A")
folds[[5]] <- which(features_normoxia$Metadata_Cell_line != "1016.02" & features_normoxia$Metadata_Cell_line != "1A" & features_normoxia$Metadata_Cell_line != "13A")
folds[[6]] <- which(features_normoxia$Metadata_Cell_line != "1002.01" & features_normoxia$Metadata_Cell_line != "3B" & features_normoxia$Metadata_Cell_line != "1006.02")
folds[[7]] <- which(features_normoxia$Metadata_Cell_line != "7A" & features_normoxia$Metadata_Cell_line != "2A" & features_normoxia$Metadata_Cell_line != "12A")
folds[[8]] <- which(features_normoxia$Metadata_Cell_line != "1009.04" & features_normoxia$Metadata_Cell_line != "4B")
folds[[9]] <- which(features_normoxia$Metadata_Cell_line != "15B" & features_normoxia$Metadata_Cell_line != "16B")
folds[[10]] <- which(features_normoxia$Metadata_Cell_line != "1013.04" & features_normoxia$Metadata_Cell_line != "17A")
folds[[11]] <- which(features_normoxia$Metadata_Cell_line != "CT29" & features_normoxia$Metadata_Cell_line != "CT30")

#Full dataset
trainingSex <- features_normoxia[,-c(1:3,5:10)]

train_control <- trainControl(method = "cv",
                              index = folds,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              savePredictions = TRUE)

AUC5 <- data.frame(matrix(ncol=3,nrow=11))
colnames(AUC5) <- c("ranger","glmnet","mlp")

for (i in seq_along(models)) {
  set.seed(120)
  modelSex <- train(
    Metadata_Sex ~ .,
    data = trainingSex,
    method = models[i],
    trControl = train_control,
    metric = "ROC",
  )
  
  AUC5[,i] <- modelSex$resample$ROC
}

write.csv(AUC,file="caret_Sex_Normoxia_DINOv3.csv")

#Hypoxia
features_hypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]

AUC6 <- data.frame(matrix(ncol=3,nrow=11))
colnames(AUC6) <- c("ranger","glmnet","mlp")

folds <- list()
folds[[1]] <- which(features_hypoxia$Metadata_Cell_line != "19A" & features_hypoxia$Metadata_Cell_line != "9B")
folds[[2]] <- which(features_hypoxia$Metadata_Cell_line != "9c1" & features_hypoxia$Metadata_Cell_line != "8B")
folds[[3]] <- which(features_hypoxia$Metadata_Cell_line != "6A" & features_hypoxia$Metadata_Cell_line != "11A")
folds[[4]] <- which(features_hypoxia$Metadata_Cell_line != "7c6" & features_hypoxia$Metadata_Cell_line != "20A")
folds[[5]] <- which(features_hypoxia$Metadata_Cell_line != "1016.02" & features_hypoxia$Metadata_Cell_line != "1A" & features_hypoxia$Metadata_Cell_line != "13A")
folds[[6]] <- which(features_hypoxia$Metadata_Cell_line != "1002.01" & features_hypoxia$Metadata_Cell_line != "3B" & features_hypoxia$Metadata_Cell_line != "1006.02")
folds[[7]] <- which(features_hypoxia$Metadata_Cell_line != "7A" & features_hypoxia$Metadata_Cell_line != "2A" & features_hypoxia$Metadata_Cell_line != "12A")
folds[[8]] <- which(features_hypoxia$Metadata_Cell_line != "1009.04" & features_hypoxia$Metadata_Cell_line != "4B")
folds[[9]] <- which(features_hypoxia$Metadata_Cell_line != "15B" & features_hypoxia$Metadata_Cell_line != "16B")
folds[[10]] <- which(features_hypoxia$Metadata_Cell_line != "1013.04" & features_hypoxia$Metadata_Cell_line != "17A")
folds[[11]] <- which(features_hypoxia$Metadata_Cell_line != "CT29" & features_hypoxia$Metadata_Cell_line != "CT30")

train_control <- trainControl(method = "cv",
                              index = folds,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              savePredictions = TRUE)

#Full dataset
trainingSex <- features_hypoxia[,-c(1:3,5:10)]

for (i in seq_along(models)) {
  set.seed(120)
  modelSex <- train(
    Metadata_Sex ~ .,
    data = trainingSex,
    method = models[i],
    trControl = train_control,
    metric = "ROC",
  )
  
  AUC6[,i] <- modelSex$resample$ROC
}

write.csv(AUC,file="caret_Sex_Hypoxia_DINOv3.csv")
