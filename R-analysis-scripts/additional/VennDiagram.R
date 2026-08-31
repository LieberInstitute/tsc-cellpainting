library(VennDiagram)
library(magrittr)

windowsFonts(Arial = windowsFont("Arial"))

#Mann-Whitney features
setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\formal-Dx-reanalysis\\Wilcoxon-rank-sum-test")

Male_hypoxia = read.csv("feature-pvalues-hypoxia-DxM-well.csv") %>% .[,c(2,4)] %>%  .[.$Adjusted_P < 0.05,]
Female_hypoxia = read.csv("feature-pvalues-hypoxia-DxF-well.csv") %>% .[,c(2,4)] %>% .[.$Adjusted_P < 0.05,]
Male_normoxia = read.csv("feature-pvalues-normoxia-DxM-well.csv") %>% .[,c(2,4)] %>% .[.$Adjusted_P < 0.05,]
Female_normoxia = read.csv("feature-pvalues-normoxia-DxF-well.csv") %>% .[,c(2,4)] %>% .[.$Adjusted_P < 0.05,]

Sex_hypoxia = read.csv("feature-pvalues-hypoxia-Sex-well.csv") %>% .[,c(2,4)] %>%  .[.$Adjusted_P < 0.05,]
Sex_normoxia = read.csv("feature-pvalues-normoxia-Sex-well.csv") %>% .[,c(2,4)] %>%  .[.$Adjusted_P < 0.05,]

feature_list <- list(
  MaleHyp = Male_hypoxia[,1],
  FemaleHyp = Female_hypoxia[,1],
  MaleNorm = Male_normoxia[,1],
  FemaleNorm = Female_normoxia[,1],
  SexHyp = Sex_hypoxia[,1],
  SexNorm = Sex_normoxia[,1]
)

overlap <- calculate.overlap(feature_list[c(1:4)])
all_overlap <- overlap[[1]]

diagram <- venn.diagram(
  x = feature_list[c(1:4)],
  category.names = c("Male\nHypoxia", "Female\nHypoxia", "Male\nNormoxia", "Female\nNormoxia"),
  filename = "venndiagram_MannWhitney_Dx.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56", "#d6584f"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 1.5,
  cat.cex = 1.3
)

diagram2 <- venn.diagram(
  x = feature_list[c(3,4,6)],
  category.names = c("Dx (male)", "Dx (female)", "Sex"),
  filename = "venndiagram_MannWhitney_DxSex_norm.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2
)

diagram3 <- venn.diagram(
  x = feature_list[c(1,2,5)],
  category.names = c("Dx (male)", "Dx (female)", "Sex"),
  filename = "venndiagram_MannWhitney_DxSex_hyp.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2
)

#Importance values
setwd("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\formal-Dx-reanalysis\\importance")

Male_hypoxia2 = read.csv("Dx_Male_hypoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:50,]
Female_hypoxia2 = read.csv("Dx_Female_hypoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:150,]
Male_normoxia2 = read.csv("Dx_Male_normoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:200,]
Female_normoxia2 = read.csv("Dx_Female_normoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:150,]

Sex_hypoxia2 = read.csv("Sex\\Sex_hypoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:250,]
Sex_normoxia2 = read.csv("Sex\\Sex_normoxia_Importance_ranger.csv") %>% .[order(-.$Overall),] %>% .[1:250,]


feature_list2 <- list(
  MaleHyp = Male_hypoxia2[,1],
  FemaleHyp = Female_hypoxia2[,1],
  MaleNorm = Male_normoxia2[,1],
  FemaleNorm = Female_normoxia2[,1],
  SexHyp = Sex_hypoxia2[,1],
  SexNorm = Sex_normoxia2[,1]
)

overlap <- calculate.overlap(feature_list[c(1:4)])
all_overlap <- overlap[[1]]

diagram4 <- venn.diagram(
  x = feature_list2[c(1:4)],
  category.names = c("Male\nHypoxia", "Female\nHypoxia", "Male\nNormoxia", "Female\nNormoxia"),
  filename = "venndiagram_importance_Dx.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56", "#d6584f"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 1.5,
  cat.cex = 1.2
)

diagram5 <- venn.diagram(
  x = feature_list2[c(3,4,6)],
  category.names = c("Dx (male)", "Dx (female)", "Sex"),
  filename = "venndiagram_importance_DxSex_norm.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2
)

diagram6 <- venn.diagram(
  x = feature_list2[c(1,2,5)],
  category.names = c("Dx (male)", "Dx (female)", "Sex"),
  filename = "venndiagram_importance_DxSex_hyp.png",
  fill = c("#3a86d4", "#f4a259", "#5fad56"),
  alpha = 0.5,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  fontface = "bold",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2
)