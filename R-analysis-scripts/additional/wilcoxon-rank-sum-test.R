library(dplyr)
library(magrittr)

setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\AUC-analysis-w-direction\\Wilcoxon-rank-sum-test\\")

#Import data
features = read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\CellProfiler analysis\\all-oxygen_triplicate_feature-selected-manually-curated.csv")
features <- features[features$Metadata_Cell_line != 'CT29' & features$Metadata_Cell_line != 'CT30',]
features[1:9] <- lapply(features[1:9],as.factor)

#Data subsets
featuresHypoxia <- features[features$Metadata_Oxygen == 'Hypoxia',]
featuresNormoxia <- features[features$Metadata_Oxygen == 'Normoxia',]

featuresHypoxiaM <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Male',]
featuresHypoxiaF <- featuresHypoxia[featuresHypoxia$Metadata_Sex == 'Female',]
featuresNormoxiaM <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Male',]
featuresNormoxiaF <- featuresNormoxia[featuresNormoxia$Metadata_Sex == 'Female',]

#wilcoxon Rank-Sum Test for diagnosis - well level (normoxia)
WRST_normoxia_Dx <- featuresNormoxia[,-c(1:3,5:9)]

p_values1 <- sapply(WRST_normoxia_Dx[,2:441], function(x) {
  wilcox.test(x ~ WRST_normoxia_Dx$Metadata_Dx)$p.value
})
adjusted_p1 <- p.adjust(p_values1, method = "BH")

results1 <- data.frame(
  Variable = names(p_values1),
  Raw_P = p_values1,
  Adjusted_P = adjusted_p1
)

#wilcoxon Rank-Sum Test for sex - well level (normoxia)
WRST_normoxia_Sex <- featuresNormoxia[,-c(1:2,4:9)]

p_values3 <- sapply(WRST_normoxia_Sex[,2:441], function(x) {
  wilcox.test(x ~ WRST_normoxia_Sex$Metadata_Sex)$p.value
})
adjusted_p3 <- p.adjust(p_values3, method = "BH")

results3 <- data.frame(
  Variable = names(p_values3),
  Raw_P = p_values3,
  Adjusted_P = adjusted_p3
)

#wilcoxon Rank-Sum Test for diagnosis - well level (hypoxia)
WRST_hypoxia_Dx <- featuresHypoxia[,-c(1:3,5:9)]

p_values5 <- sapply(WRST_hypoxia_Dx[,2:441], function(x) {
  wilcox.test(x ~ WRST_hypoxia_Dx$Metadata_Dx)$p.value
})
adjusted_p5 <- p.adjust(p_values5, method = "BH")

results5 <- data.frame(
  Variable = names(p_values5),
  Raw_P = p_values5,
  Adjusted_P = adjusted_p5
)

#wilcoxon Rank-Sum Test for sex - well level (hypoxia)
WRST_hypoxia_Sex <- featuresHypoxia[,-c(1:2,4:9)]

p_values7 <- sapply(WRST_hypoxia_Sex[,2:441], function(x) {
  wilcox.test(x ~ WRST_hypoxia_Sex$Metadata_Sex)$p.value
})
adjusted_p7 <- p.adjust(p_values7, method = "BH")

results7 <- data.frame(
  Variable = names(p_values7),
  Raw_P = p_values7,
  Adjusted_P = adjusted_p7
)

write.csv(results1,file="feature-pvalues-normoxia-Dx-well.csv")
write.csv(results3,file="feature-pvalues-normoxia-Sex-well.csv")
write.csv(results5,file="feature-pvalues-hypoxia-Dx-well.csv")
write.csv(results7,file="feature-pvalues-hypoxia-Sex-well.csv")


#wilcoxon Rank-Sum Test for diagnosis - well level (normoxia, MALE)
WRST_normoxia_DxM <- featuresNormoxiaM[,-c(1:3,5:9)]

p_values9 <- sapply(WRST_normoxia_DxM[,2:441], function(x) {
  wilcox.test(x ~ WRST_normoxia_DxM$Metadata_Dx)$p.value
})
adjusted_p9 <- p.adjust(p_values9, method = "BH")

results9 <- data.frame(
  Variable = names(p_values9),
  Raw_P = p_values9,
  Adjusted_P = adjusted_p9
)

#wilcoxon Rank-Sum Test for diagnosis - well level (normoxia, FEMALE)
WRST_normoxia_DxF <- featuresNormoxiaF[,-c(1:3,5:9)]

p_values10 <- sapply(WRST_normoxia_DxF[,2:441], function(x) {
  wilcox.test(x ~ WRST_normoxia_DxF$Metadata_Dx)$p.value
})
adjusted_p10 <- p.adjust(p_values10, method = "BH")

results10 <- data.frame(
  Variable = names(p_values10),
  Raw_P = p_values10,
  Adjusted_P = adjusted_p10
)

#wilcoxon Rank-Sum Test for diagnosis - well level (hypoxia, MALE)
WRST_hypoxia_DxM <- featuresHypoxiaM[,-c(1:3,5:9)]

p_values11 <- sapply(WRST_hypoxia_DxM[,2:441], function(x) {
  wilcox.test(x ~ WRST_hypoxia_DxM$Metadata_Dx)$p.value
})
adjusted_p11 <- p.adjust(p_values11, method = "BH")

results11 <- data.frame(
  Variable = names(p_values11),
  Raw_P = p_values11,
  Adjusted_P = adjusted_p11
)

#wilcoxon Rank-Sum Test for diagnosis - well level (hypoxia, MALE)
WRST_hypoxia_DxF <- featuresHypoxiaF[,-c(1:3,5:9)]

p_values12 <- sapply(WRST_hypoxia_DxF[,2:441], function(x) {
  wilcox.test(x ~ WRST_hypoxia_DxF$Metadata_Dx)$p.value
})
adjusted_p12 <- p.adjust(p_values12, method = "BH")

results12 <- data.frame(
  Variable = names(p_values12),
  Raw_P = p_values12,
  Adjusted_P = adjusted_p12
)

write.csv(results9,file="feature-pvalues-normoxia-DxM-well.csv")
write.csv(results10,file="feature-pvalues-normoxia-DxF-well.csv")
write.csv(results11,file="feature-pvalues-hypoxia-DxM-well.csv")
write.csv(results12,file="feature-pvalues-hypoxia-DxF-well.csv")
