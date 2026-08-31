setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\formal-Dx-reanalysis\\importance\\")

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------Features----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#Normoxia
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\CellProfiler analysis\\all-oxygen_triplicate_feature-selected-manually-curated.csv")
features <- features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]
features[1:9] <- lapply(features[1:9],as.factor)

features_normoxia <- features[features$Metadata_Oxygen == 'Normoxia',]
featuresNormoxiaM <- features_normoxia[features_normoxia$Metadata_Sex == 'Male',]
featuresNormoxiaF <- features_normoxia[features_normoxia$Metadata_Sex == 'Female',]

models <- c("ranger")
#models <- c("glmnet","mlp")

folds <- list()
folds[[1]] <- which(features_normoxia$Metadata_Cell_line != "16B" & features_normoxia$Metadata_Cell_line != "17A")
folds[[2]] <- which(features_normoxia$Metadata_Cell_line != "9B" & features_normoxia$Metadata_Cell_line != "20A")
folds[[3]] <- which(features_normoxia$Metadata_Cell_line != "8B" & features_normoxia$Metadata_Cell_line != "11A")
folds[[4]] <- which(features_normoxia$Metadata_Cell_line != "3B" & features_normoxia$Metadata_Cell_line != "1A")
folds[[5]] <- which(features_normoxia$Metadata_Cell_line != "4B" & features_normoxia$Metadata_Cell_line != "2A")
folds[[6]] <- which(features_normoxia$Metadata_Cell_line != "1009.04" & features_normoxia$Metadata_Cell_line != "1002.01" & features_normoxia$Metadata_Cell_line != "7A")
folds[[7]] <- which(features_normoxia$Metadata_Cell_line != "1016.02" & features_normoxia$Metadata_Cell_line != "1006.02" & features_normoxia$Metadata_Cell_line != "6A")
folds[[8]] <- which(features_normoxia$Metadata_Cell_line != "15B" & features_normoxia$Metadata_Cell_line != "1013.04" & features_normoxia$Metadata_Cell_line != "12A")
folds[[9]] <- which(features_normoxia$Metadata_Cell_line != "9c1" & features_normoxia$Metadata_Cell_line != "19A")
folds[[10]] <- which(features_normoxia$Metadata_Cell_line != "7c6" & features_normoxia$Metadata_Cell_line != "13A")

#Full dataset
trainingDx <- features_normoxia[,-c(1:3,5:9)]

train_control <- trainControl(method = "cv",
                              index = folds,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              savePredictions = TRUE)

for (i in seq_along(models)) {
  set.seed(120)
  modelDx <- train(
    Metadata_Dx ~ .,
    data = trainingDx,
    method = models[i],
    trControl = train_control,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- (varImp(modelDx, scale = FALSE))
  write.csv(DxImportance$importance,file=paste0("Dx_normoxia_Importance_",models[i],".csv"))
}

trainingDxF <- featuresNormoxiaF[,-c(1:3,5:9)]
trainingDxM <- featuresNormoxiaM[,-c(1:3,5:9)]

foldsF <- list()
foldsF[[1]] <- which(featuresNormoxiaF$Metadata_Cell_line != "16B" & featuresNormoxiaF$Metadata_Cell_line != "17A")
foldsF[[2]] <- which(featuresNormoxiaF$Metadata_Cell_line != "9B" & featuresNormoxiaF$Metadata_Cell_line != "20A")
foldsF[[3]] <- which(featuresNormoxiaF$Metadata_Cell_line != "8B" & featuresNormoxiaF$Metadata_Cell_line != "11A")
foldsF[[4]] <- which(featuresNormoxiaF$Metadata_Cell_line != "3B" & featuresNormoxiaF$Metadata_Cell_line != "1A")
foldsF[[5]] <- which(featuresNormoxiaF$Metadata_Cell_line != "4B" & featuresNormoxiaF$Metadata_Cell_line != "2A")

foldsM <- list()
foldsM[[1]] <- which(featuresNormoxiaM$Metadata_Cell_line != "1009.04" & featuresNormoxiaM$Metadata_Cell_line != "1002.01" & featuresNormoxiaM$Metadata_Cell_line != "7A")
foldsM[[2]] <- which(featuresNormoxiaM$Metadata_Cell_line != "1016.02" & featuresNormoxiaM$Metadata_Cell_line != "1006.02" & featuresNormoxiaM$Metadata_Cell_line != "6A")
foldsM[[3]] <- which(featuresNormoxiaM$Metadata_Cell_line != "15B" & featuresNormoxiaM$Metadata_Cell_line != "1013.04" & featuresNormoxiaM$Metadata_Cell_line != "12A")
foldsM[[4]] <- which(featuresNormoxiaM$Metadata_Cell_line != "9c1" & featuresNormoxiaM$Metadata_Cell_line != "19A")
foldsM[[5]] <- which(featuresNormoxiaM$Metadata_Cell_line != "7c6" & featuresNormoxiaM$Metadata_Cell_line != "13A")

train_controlF <- trainControl(method = "cv",
                               index = foldsF,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

train_controlM <- trainControl(method = "cv",
                               index = foldsM,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

for (i in seq_along(models)) {
  set.seed(120)
  modelDxF <- train(
    Metadata_Dx ~ .,
    data = trainingDxF,
    method = models[i],
    trControl = train_controlF,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- varImp(modelDxF, scale = FALSE)
  write.csv(DxImportance$importance,file=paste0("Dx_Female_normoxia_Importance_",models[i],".csv"))
}

for (i in seq_along(models)) {
  set.seed(120)
  modelDxM <- train(
    Metadata_Dx ~ .,
    data = trainingDxM,
    method = models[i],
    trControl = train_controlM,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- varImp(modelDxM, scale = FALSE)
  write.csv(DxImportance$importance,file=paste0("Dx_Male_normoxia_Importance_",models[i],".csv"))
}

#Hypoxia
features_hypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]
featuresHypoxiaM <- features_hypoxia[features_hypoxia$Metadata_Sex == 'Male',]
featuresHypoxiaF <- features_hypoxia[features_hypoxia$Metadata_Sex == 'Female',]

folds <- list()
folds[[1]] <- which(features_hypoxia$Metadata_Cell_line != "16B" & features_hypoxia$Metadata_Cell_line != "17A")
folds[[2]] <- which(features_hypoxia$Metadata_Cell_line != "9B" & features_hypoxia$Metadata_Cell_line != "20A")
folds[[3]] <- which(features_hypoxia$Metadata_Cell_line != "8B" & features_hypoxia$Metadata_Cell_line != "11A")
folds[[4]] <- which(features_hypoxia$Metadata_Cell_line != "3B" & features_hypoxia$Metadata_Cell_line != "1A")
folds[[5]] <- which(features_hypoxia$Metadata_Cell_line != "4B" & features_hypoxia$Metadata_Cell_line != "2A")
folds[[6]] <- which(features_hypoxia$Metadata_Cell_line != "1009.04" & features_hypoxia$Metadata_Cell_line != "1002.01" & features_hypoxia$Metadata_Cell_line != "7A")
folds[[7]] <- which(features_hypoxia$Metadata_Cell_line != "1016.02" & features_hypoxia$Metadata_Cell_line != "1006.02" & features_hypoxia$Metadata_Cell_line != "6A")
folds[[8]] <- which(features_hypoxia$Metadata_Cell_line != "15B" & features_hypoxia$Metadata_Cell_line != "1013.04" & features_hypoxia$Metadata_Cell_line != "12A")
folds[[9]] <- which(features_hypoxia$Metadata_Cell_line != "9c1" & features_hypoxia$Metadata_Cell_line != "19A")
folds[[10]] <- which(features_hypoxia$Metadata_Cell_line != "7c6" & features_hypoxia$Metadata_Cell_line != "13A")

foldsF <- list()
foldsF[[1]] <- which(featuresHypoxiaF$Metadata_Cell_line != "16B" & featuresHypoxiaF$Metadata_Cell_line != "17A")
foldsF[[2]] <- which(featuresHypoxiaF$Metadata_Cell_line != "9B" & featuresHypoxiaF$Metadata_Cell_line != "20A")
foldsF[[3]] <- which(featuresHypoxiaF$Metadata_Cell_line != "8B" & featuresHypoxiaF$Metadata_Cell_line != "11A")
foldsF[[4]] <- which(featuresHypoxiaF$Metadata_Cell_line != "3B" & featuresHypoxiaF$Metadata_Cell_line != "1A")
foldsF[[5]] <- which(featuresHypoxiaF$Metadata_Cell_line != "4B" & featuresHypoxiaF$Metadata_Cell_line != "2A")

foldsM <- list()
foldsM[[1]] <- which(featuresHypoxiaM$Metadata_Cell_line != "1009.04" & featuresHypoxiaM$Metadata_Cell_line != "1002.01" & featuresHypoxiaM$Metadata_Cell_line != "7A")
foldsM[[2]] <- which(featuresHypoxiaM$Metadata_Cell_line != "1016.02" & featuresHypoxiaM$Metadata_Cell_line != "1006.02" & featuresHypoxiaM$Metadata_Cell_line != "6A")
foldsM[[3]] <- which(featuresHypoxiaM$Metadata_Cell_line != "15B" & featuresHypoxiaM$Metadata_Cell_line != "1013.04" & featuresHypoxiaM$Metadata_Cell_line != "12A")
foldsM[[4]] <- which(featuresHypoxiaM$Metadata_Cell_line != "9c1" & featuresHypoxiaM$Metadata_Cell_line != "19A")
foldsM[[5]] <- which(featuresHypoxiaM$Metadata_Cell_line != "7c6" & featuresHypoxiaM$Metadata_Cell_line != "13A")

train_control <- trainControl(method = "cv",
                              index = folds,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              savePredictions = TRUE)

train_controlF <- trainControl(method = "cv",
                               index = foldsF,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

train_controlM <- trainControl(method = "cv",
                               index = foldsM,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               savePredictions = TRUE)

#Full dataset
trainingDx <- features_hypoxia[,-c(1:3,5:9)]

for (i in seq_along(models)) {
  set.seed(120)
  modelDx <- train(
    Metadata_Dx ~ .,
    data = trainingDx,
    method = models[i],
    trControl = train_control,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- varImp(modelDx, scale = FALSE)
  write.csv(DxImportance$importance,file=paste0("Dx_hypoxia_Importance_",models[i],".csv"))
}

trainingDxF <- featuresHypoxiaF[,-c(1:3,5:9)]
trainingDxM <- featuresHypoxiaM[,-c(1:3,5:9)]

for (i in seq_along(models)) {
  set.seed(120)
  modelDxF <- train(
    Metadata_Dx ~ .,
    data = trainingDxF,
    method = models[i],
    trControl = train_controlF,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- varImp(modelDxF, scale = FALSE)
  write.csv(DxImportance$importance,file=paste0("Dx_Female_hypoxia_Importance_",models[i],".csv"))
}

for (i in seq_along(models)) {
  set.seed(120)
  modelDxM <- train(
    Metadata_Dx ~ .,
    data = trainingDxM,
    method = models[i],
    trControl = train_controlM,
    metric = "ROC",
    importance = "impurity"
#   importance = TRUE
  )
  
  DxImportance <- varImp(modelDxM, scale = FALSE)
  write.csv(DxImportance$importance,file=paste0("Dx_Male_hypoxia_Importance_",models[i],".csv"))
}