## Jason Code

library(FlowSOM)
library(flowCore)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(grid)
library(gridExtra)
library(ggridges)
library("ggpubr")
#file <- "EP2012-07-11-WT-8 weeks_4.fcs"
#flowO <- flowCore::read.FCS(file,transformation=FALSE)
#file <- "EP20150815-InjuryTimeCourse-Day6_3.fcs"
file <- "EP20150513-InjuryTimeCourse-Day6_1.fcs"
file <- "EP20150513-InjuryTimeCourse-Day6_1.fcs"
flowO <- flowCore::read.FCS(file,transformation=FALSE)
data <- flowO@exprs
#(flowO)
attributes <-flowO@parameters@data
# list all cells
list <- 1:dim(data)[1]
set.seed(42)
# random permutation of cells
newOrder = sample(list,dim(data)[1],replace = FALSE, prob = NULL)

#generating all data except four markers to visualize
listElse <- 1
listElse[1:6] <- 1:6
listElse[7:10] <- 8:12
listElse[11:19] <- 14:23
listElse[20:23] <- 25:29

CD9 <- data[,28]
CD104 <- data[,5]

Pax7 <- data[,13]
Myf5 <- data[,24]
MyoD <- data[,7]
Myogenin <- data[,11]

scale <- "lnear"
if(scale=="linear")
{
  CD9 <- CD9
  CD10 <- CD10
  Pax7 <- Pax7
  Myf5 <- Myf5
  MyoD <- MyoD
  Myogenin <- Myogenin
}
if (scale=="asinh")
{
  CD9 <- asinh(CD9)
  CD10 <- asin(CD10)
  Pax7 <- asinh(Pax7)
  Myf5 <- asinh(Myf5)
  MyoD <- asinh(MyoD)
  Myogenin <- asinh(Myogenin)
}
if (scale=="log")
{
  CD9 <- log10(CD9+1)
  CD10 <- log10(CD10+1)
  Pax7 <- log10(Pax7+1)
  Myf5 <- log10(Myf5+1)
  MyoD <- log10(MyoD+1)
  Myogenin <- log10(Myogenin)
}

mycol <- colorpanel(100,"blue","yellow","red")

set.seed(43)
# Run UMAP on all scaled surface markers
data2 <- data[,listElse]
myumap <-umap(data2,n_neighbors = 15,n_threads = 1, verbose = TRUE)
umap.data = as.data.frame(myumap)

plot.new()
frame()
p1 <- qplot(CD104,CD9, col = Pax7)
p2 <- qplot(CD104,CD9, col = Myf5)
p3 <- qplot(CD104,CD9, col = MyoD)
p4 <- qplot(CD104,CD9, col = Myogenin)
grid.arrange(p1, p2, p3, p4, nrow = 1)


set.seed(43)
# Run UMAP on all scaled surface markers
data2 <- data[,listElse]
if(scale=="linear")
{
  data2 <- data2
}
if (scale=="asinh")
{
  data2 <- asinh(data2)
}
if (scale=="log")
{
  data2 <- log10(data2+1)
}
myumap <-umap(data2,n_neighbors = 15,n_threads = 1, verbose = TRUE)
umap.data = as.data.frame(myumap)

plot.new()
frame()
p1 <- plot(umap.data[,1],umap.data[,2], col = Pax7)
p2 <- plot(umap.data[,1],umap.data[,2], col = Myf5)
p3 <- plot(umap.data[,1],umap.data[,2], col = MyoD)
p4 <- plot(umap.data[,1],umap.data[,2], col = Myogenin)

grid.arrange(p1, p2, p3, p4,labels = c("Pax7", "Myf5", "MyoD","Myogenin"),nrow = 1)
  