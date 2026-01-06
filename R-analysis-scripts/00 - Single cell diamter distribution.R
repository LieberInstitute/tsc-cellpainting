install.packages("ggplot2","scales")
library(ggplot2)
library(scales)

#Import data
rna <- read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\rna_diameters.csv")
dna <- read.csv("C:\\Users\\frank.piscotta\\OneDrive - Lieber Institute for Brain Development\\Documents\\Cell Painting\\TSC Project\\TSC paper\\nucleus_diameters.csv")

p <- ggplot() +
  geom_histogram(data = dna, aes(x = diameter), fill = "blue", color = "black", binwidth = 1) +
  geom_histogram(data = rna, aes(x = diameter), fill = "yellow", color = "black", binwidth = 1) +
  theme_bw() +
  scale_y_continuous(name = "Frequency", breaks = seq(0,600000,100000), limits = c(0,600000), labels = label_comma(), expand = expansion(add = 0)) +
  scale_x_continuous(name = "Diameter (px)", breaks = seq(0,140,20), limits = c(0,140), expand = expansion(add = 0)) +
  geom_vline(xintercept = 82, linetype = "dashed", color = "red",size = 1) +
  theme(axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=16,face="bold")
  )

ggsave("diameter_distribution2.tiff", p, width = 6, height = 6, units = 'in', dpi = 800, compression = "lzw")