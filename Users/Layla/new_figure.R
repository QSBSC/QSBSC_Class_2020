library(gplots) 
library(ggfortify)
library(ggplot2)

# load data
d <- read.csv('../../Data/Porpiglia/FR-FCM-ZY3K_gated.csv', as.is=TRUE)

ggplot2::autoplot(stats::prcomp(d[,c(4,11,16,17,18,19)], scale=TRUE), 
                  label = FALSE, loadings.label = TRUE)
