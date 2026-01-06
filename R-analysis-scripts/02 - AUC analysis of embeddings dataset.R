install.packages(c("tidyverse","magrittr","randomForest","pROC","tidymodels","glmnet","neuralnet"))

library(tidyverse)
library(magrittr)
library(randomForest)
library(pROC)
library(tidymodels)
library(glmnet)
library(neuralnet)

DxAnalysis <- function(data,folds) {
  
  AUC <- data.frame(matrix(ncol=3,nrow=nrow(folds)))
  colnames(AUC) <- c("RF","LR","MLP")

  for (i in 1:nrow(folds)) {
    set.seed(120)
    
    test <- data[data$Metadata_Cell_line %in% folds[i,],] %>% .[,-c(1:4,6:10)]
    train <-data[!(data$Metadata_Cell_line %in% folds[i,]),] %>% .[,-c(1:4,6:10)]
    
    #Random Forest
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index = 2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test$Metadata_Dx,y_pred$Scz)
    AUC[i,1] <- auc(RF_roc)
    
    #Logistic regression
    model <- logistic_reg(mixture = double(1),penalty = double(1)) %>%
      set_engine("glmnet") %>%
      set_mode("classification") %>%
      fit(Metadata_Dx ~ ., data = train)
    
    results <- tidy(model)
    
    pred_class <- predict(model,
                          new_data = test,
                          type = "class")
    pred_proba <- predict(model,
                          new_data = test,
                          type = "prob")
    
    results <- test %>%
      select(Metadata_Dx) %>%
      bind_cols(pred_class,pred_proba)
    
    pred_proba <- data.frame(pred_proba)
    logreg_roc <- roc(test$Metadata_Dx,pred_proba$.pred_Scz)
    AUC[i,2] <- auc(logreg_roc)
    
    #Multilayer Perceptron
    train$Metadata_Dx <- as.numeric(train$Metadata_Dx) - 1
    test$Metadata_Dx <- as.numeric(test$Metadata_Dx) - 1
    
    model <- neuralnet(Metadata_Dx~., data = train, hidden = c(5,3), linear.output = FALSE)
    predictions <- predict(model, test)
    
    MLP_roc <- roc(test$Metadata_Dx,as.vector(predictions))
    AUC[i,3] <- auc(MLP_roc)
  }
  
  return(AUC)
}

SexAnalysis <- function(data) {
  
  AUC <- data.frame(matrix(ncol=3,nrow=10))
  colnames(AUC) <- c("RF","LR","MLP")
  
  for (i in 1:10) {
    set.seed(120)
    
    #Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10)]
    train <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10)]
    
    #Random Forest
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test$Metadata_Sex,y_pred$Male)
    AUC[i,1] <- auc(RF_roc)
    
    #Logistic regression
    model <- logistic_reg(mixture = double(1),penalty = double(1)) %>%
      set_engine("glmnet") %>%
      set_mode("classification") %>%
      fit(Metadata_Sex ~ ., data = train)
    
    results <- tidy(model)
    
    pred_class <- predict(model,
                          new_data = test,
                          type = "class")
    pred_proba <- predict(model,
                          new_data = test,
                          type = "prob")
    
    results <- test %>%
      select(Metadata_Sex) %>%
      bind_cols(pred_class,pred_proba)
    
    pred_proba <- data.frame(pred_proba)
    logreg_roc <- roc(test$Metadata_Sex,pred_proba$.pred_Male)
    AUC[i,2] <- auc(logreg_roc)
    
    #Multilayer Perceptron
    train$Metadata_Sex <- as.numeric(train$Metadata_Sex) - 1
    test$Metadata_Sex <- as.numeric(test$Metadata_Sex) - 1
    
    model <- neuralnet(Metadata_Sex~., data = train, hidden = c(5,3), linear.output = FALSE)
    predictions <- predict(model, test)
    
    MLP_roc <- roc(test$Metadata_Sex,as.vector(predictions))
    AUC[i,3] <- auc(MLP_roc)
  }
  
  return(AUC)
}

singleChannelDxRF <- function(data) {
  
  AUC <- data.frame(matrix(ncol=6,nrow=10))
  colnames(AUC) <- c("All","AGP","DNA","ER","Mito","RNA")

  for (i in 1:10) {
    set.seed(120)
    
    #Create single channel datasets - Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10)]
    test_AGP <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10,139:650)]
    test_DNA <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10,11:138,267:650)]
    test_ER <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10,11:266,395:650)]
    test_Mito <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10,11:394,523:650)]
    test_RNA <-data[data$Metadata_Cell_line %in% foldsDx[i,],] %>% .[,-c(1:4,6:10,11:522)]
    
    train <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10)]
    train_AGP <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10,139:650)]
    train_DNA <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10,11:138,267:650)]
    train_ER <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10,11:266,395:650)]
    train_Mito <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10,11:394,523:650)]
    train_RNA <-data[!(data$Metadata_Cell_line %in% foldsDx[i,]),] %>% .[,-c(1:4,6:10,11:522)]
    
    #Random Forest
    #All
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test$Metadata_Dx,y_pred$Scz)
    AUC[i,1] <- auc(RF_roc)
    
    #AGP
    classifier_RF = randomForest(x = train_AGP[-1],
                                 y = train_AGP$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_AGP$Metadata_Dx,y_pred$Scz)
    AUC[i,2] <- auc(RF_roc)
    
    #DNA
    classifier_RF = randomForest(x = train_DNA[-1],
                                 y = train_DNA$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_DNA$Metadata_Dx,y_pred$Scz)
    AUC[i,3] <- auc(RF_roc)
    
    #ER
    classifier_RF = randomForest(x = train_ER[-1],
                                 y = train_ER$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_ER$Metadata_Dx,y_pred$Scz)
    AUC[i,4] <- auc(RF_roc)
    
    #Mito
    classifier_RF = randomForest(x = train_Mito[-1],
                                 y = train_Mito$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_Mito$Metadata_Dx,y_pred$Scz)
    AUC[i,5] <- auc(RF_roc)
    
    #RNA
    classifier_RF = randomForest(x = train_RNA[-1],
                                 y = train_RNA$Metadata_Dx,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_RNA$Metadata_Dx,y_pred$Scz)
    AUC[i,6] <- auc(RF_roc)
  }
  
  return(AUC)
}

singleChannelSexRF <- function(data) {
  
  AUC <- data.frame(matrix(ncol=6,nrow=10))
  colnames(AUC) <- c("All","AGP","DNA","ER","Mito","RNA")

  for (i in 1:10) {
    set.seed(120)
    
    #Create single channel datasets - Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10)]
    test_AGP <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10,139:650)]
    test_DNA <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10,11:138,267:650)]
    test_ER <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10,11:266,395:650)]
    test_Mito <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10,11:394,523:650)]
    test_RNA <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:3,5:10,11:522)]
    
    train <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10)]
    train_AGP <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10,139:650)]
    train_DNA <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10,11:138,267:650)]
    train_ER <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10,11:266,395:650)]
    train_Mito <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10,11:394,523:650)]
    train_RNA <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:3,5:10,11:522)]
    
    #Random Forest
    #All
    classifier_RF = randomForest(x = train[-1],
                                 y = train$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test$Metadata_Sex,y_pred$Male)
    AUC[i,1] <- auc(RF_roc)
    
    #AGP
    classifier_RF = randomForest(x = train_AGP[-1],
                                 y = train_AGP$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_AGP$Metadata_Sex,y_pred$Male)
    AUC[i,2] <- auc(RF_roc)
    
    #DNA
    classifier_RF = randomForest(x = train_DNA[-1],
                                 y = train_DNA$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_DNA$Metadata_Sex,y_pred$Male)
    AUC[i,3] <- auc(RF_roc)
    
    #ER
    classifier_RF = randomForest(x = train_ER[-1],
                                 y = train_ER$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_ER$Metadata_Sex,y_pred$Male)
    AUC[i,4] <- auc(RF_roc)
    
    #Mito
    classifier_RF = randomForest(x = train_Mito[-1],
                                 y = train_Mito$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_Mito$Metadata_Sex,y_pred$Male)
    AUC[i,5] <- auc(RF_roc)
    
    #RNA
    classifier_RF = randomForest(x = train_RNA[-1],
                                 y = train_RNA$Metadata_Sex,
                                 ntree = 2500)
    
    y_pred = predict(classifier_RF, newdata = test[-1],index=2, type="prob") %>% data.frame(.)
    RF_roc <- roc(test_RNA$Metadata_Sex,y_pred$Male)
    AUC[i,6] <- auc(RF_roc)
    }

  return(AUC)
}
  
#Import data
embeddings = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\Spring analysis\\combined embeddings\\combined embeddings - all oxygen - well - normalized by plate.csv")
embeddings <- embeddings[embeddings$Metadata_Cell_line != 'CT29' & embeddings$Metadata_Cell_line != 'CT30',]
embeddings[1:10] <- lapply(embeddings[1:10],as.factor)

#Data subsets
embeddingsHypoxia <- embeddings[embeddings$Metadata_Oxygen == 'Hypoxia',]
embeddingsHypoxiaM <- embeddingsHypoxia[embeddingsHypoxia$Metadata_Sex == 'Male',]
embeddingsHypoxiaF <- embeddingsHypoxia[embeddingsHypoxia$Metadata_Sex == 'Female',]
embeddingsNormoxia <- embeddings[embeddings$Metadata_Oxygen == 'Normoxia',]
embeddingsNormoxiaM <- embeddingsNormoxia[embeddingsNormoxia$Metadata_Sex == 'Male',]
embeddingsNormoxiaF <- embeddingsNormoxia[embeddingsNormoxia$Metadata_Sex == 'Female',]

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

#Randomized labels for control test
labelsDx <- c("Scz","Ctrl")
labelsSex <- c("Male","Female")
shuffledNormoxia <- embeddingsNormoxia

set.seed(120)
shuffledNormoxia$Metadata_Dx <- sample(labelsDx,size=nrow(shuffledNormoxia),replace=TRUE, prob=c(0.6,0.4)) %>% as.factor(.)
set.seed(120)
shuffledNormoxia$Metadata_Sex <- sample(labelsSex,size=nrow(shuffledNormoxia),replace=TRUE, prob=c(0.6,0.4)) %>% as.factor(.)

#Analysis
normoxiaDxAUC <- DxAnalysis(embeddingsNormoxia,foldsDx)
hypoxiaDxAUC <- DxAnalysis(embeddingsHypoxia,foldsDx)
normoxiaDxMaleAUC <- DxAnalysis(embeddingsNormoxiaM,foldsDxM)
normoxiaDxFemaleAUC <- DxAnalysis(embeddingsNormoxiaF,foldsDxF)
hypoxiaDxMaleAUC <- DxAnalysis(embeddingsHypoxiaM,foldsDxM)
hypoxiaDxFemaleAUC <- DxAnalysis(embeddingsHypoxiaF,foldsDxF)
normoxiaSexAUC <- SexAnalysis(embeddingsNormoxia)
hypoxiaSexAUC <- SexAnalysis(embeddingsHypoxia)
normoxiaShuffledDxAUC <- DxAnalysis(shuffledNormoxia,foldsDx)
normoxiaShuffledSexAUC <- SexAnalysis(shuffledNormoxia)
normoxiaDxSingleChannelAUC <- singleChannelDxRF(embeddingsNormoxia)
hypoxiaDxSingleChannelAUC <- singleChannelDxRF(embeddingsHypoxia)
normoxiaSexSingleChannelAUC <- singleChannelSexRF(embeddingsNormoxia)
hypoxiaSexSingleChannelAUC <- singleChannelSexRF(embeddingsHypoxia)

#Export results
write.csv(normoxiaDxAUC,file="aucNormoxiaDx.csv")
write.csv(hypoxiaDxAUC,file="aucHypoxiaDx.csv")
write.csv(normoxiaDxMaleAUC,file="aucNormoxiaMaleDx.csv")
write.csv(normoxiaDxFemaleAUC,file="aucNormoxiaFemaleDx.csv")
write.csv(hypoxiaDxMaleAUC,file="aucHypoxiaMaleDx.csv")
write.csv(hypoxiaDxFemaleAUC,file="aucHypoxiaFemaleDx.csv")
write.csv(normoxiaSexAUC,file="aucNormoxiaSex.csv")
write.csv(hypoxiaSexAUC,file="aucHypoxiaSex.csv")
write.csv(normoxiaShuffledDxAUC,file="aucNormoxiaDx_random.csv")
write.csv(normoxiaShuffledSexAUC,file="aucNormoxiaSex_random.csv")
write.csv(normoxiaDxSingleChannelAUC,file="aucNormoxiaDx_singlechannel.csv")
write.csv(hypoxiaDxSingleChannelAUC,file="aucHypoxiaDx_singlechannel.csv")
write.csv(normoxiaSexSingleChannelAUC,file="aucNormoxiaSex_singlechannel.csv")
write.csv(hypoxiaSexSingleChannelAUC,file="aucHypoxiaSex_singlechannel.csv")