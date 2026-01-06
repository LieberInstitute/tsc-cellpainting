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
    
    #Split data, remove irrelevant metadata columns
    test <-data[data$Metadata_Cell_line %in% folds[i,],] %>% .[,-c(1:3,5:9)]
    train <-data[!(data$Metadata_Cell_line %in% folds[i,]),] %>% .[,-c(1:3,5:9)]
    
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
    test <-data[data$Metadata_Cell_line %in% foldsSex[i,],] %>% .[,-c(1:2,4:9)]
    train <-data[!(data$Metadata_Cell_line %in% foldsSex[i,]),] %>% .[,-c(1:2,4:9)]
    
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

#Import data
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\CellProfiler analysis\\all-oxygen_triplicate_feature-selected-manually-curated.csv")
features <- features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]
features[1:9] <- lapply(features[1:9],as.factor)

#Data subsets
featuresHypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]
featuresHypoxiaM <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Male',]
featuresHypoxiaF <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Female',]
featuresNormoxia <- features[features$Metadata_Oxygen == 'Normoxia',]
featuresNormoxiaM <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Male',]
featuresNormoxiaF <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Female',]

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
shuffledNormoxia <- featuresNormoxia

set.seed(120)
shuffledNormoxia$Metadata_Dx <- sample(labelsDx,size=nrow(shuffledNormoxia),replace=TRUE, prob=c(0.6,0.4)) %>% as.factor(.)
set.seed(120)
shuffledNormoxia$Metadata_Sex <- sample(labelsSex,size=nrow(shuffledNormoxia),replace=TRUE, prob=c(0.6,0.4)) %>% as.factor(.)

#Analysis
normoxiaDxAUC <- DxAnalysis(featuresNormoxia,foldsDx)
hypoxiaDxAUC <- DxAnalysis(featuresHypoxia,foldsDx)
normoxiaDxMaleAUC <- DxAnalysis(featuresNormoxiaM,foldsDxM)
normoxiaDxFemaleAUC <- DxAnalysis(featuresNormoxiaF,foldsDxF)
hypoxiaDxMaleAUC <- DxAnalysis(featuresHypoxiaM,foldsDxM)
hypoxiaDxFemaleAUC <- DxAnalysis(featuresHypoxiaF,foldsDxF)
normoxiaSexAUC <- SexAnalysis(featuresNormoxia)
hypoxiaSexAUC <- SexAnalysis(featuresHypoxia)
normoxiaShuffledDxAUC <- DxAnalysis(shuffledNormoxia,foldsDx)
normoxiaShuffledSexAUC <- SexAnalysis(shuffledNormoxia)

#Export
write.csv(normoxiaDxAUC,file="aucNormoxiaDx_features.csv")
write.csv(hypoxiaDxAUC,file="aucHypoxiaDx_features.csv")
write.csv(normoxiaDxMaleAUC,file="aucNormoxiaMaleDx_features.csv")
write.csv(normoxiaDxFemaleAUC,file="aucNormoxiaFemaleDx_features.csv")
write.csv(hypoxiaDxMaleAUC,file="aucHypoxiaMaleDx_features.csv")
write.csv(hypoxiaDxFemaleAUC,file="aucHypoxiaFemaleDx_features.csv")
write.csv(normoxiaSexAUC,file="aucNormoxiaSex_features.csv")
write.csv(hypoxiaSexAUC,file="aucHypoxiaSex_features.csv")
write.csv(normoxiaShuffledDxAUC,file="aucNormoxiaDx_random_features.csv")
write.csv(normoxiaShuffledSexAUC,file="aucNormoxiaSex_random_features.csv")