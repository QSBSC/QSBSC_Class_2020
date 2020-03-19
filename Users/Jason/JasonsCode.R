---
output:
  pdf_document:
  latex_engine: xelatex
html_document:
  df_print: paged
editor_options:
  chunk_output_type: inline
---
library(FlowSOM)
library(flowCore)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(grid)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(umap)

#file <- "EP2012-07-11-WT-8 weeks_4.fcs"
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

scale <- "log"
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
  CD104 <- log10(CD104+1)
  Pax7 <- log10(Pax7+1)
  Myf5 <- log10(Myf5+1)
  MyoD <- log10(MyoD+1)
  Myogenin <- log10(Myogenin+1)
}

mycol <- colorpanel(100,"blue","yellow","red")



plot.new()
frame()

getCol <- function(var)
{
  i <- ifelse(100*var/max(var) > 100, 100, round(100*var/max(var),1))
  col[i] 
}
var <- Pax7
Pax7C <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- Myf5
Myf5C <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- MyoD
MyoDC <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- Myogenin
MyogeninC <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]

par(mfrow = c(1, 4))
plot(CD104,CD9, main = "Pax7",col = Pax7C)
plot(CD104,CD9, main = "Myf5",col = Myf5C)
plot(CD104,CD9, main = "MyoD",col = MyoDC)
plot(CD104,CD9, main = "Myogenin",col = MyogeninC)

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
attributes
UMAP <- myumap[["layout"]]

plot.new()
frame()
var <- Pax7
Pax7C <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- Myf5
Myf5C <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- MyoD
MyoDC <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
var <- Myogenic
MyogenicC <-mycol[ifelse(100*(var+1)/max(var) > 100, 100, round(100*(var+1)/max(var),1))]
par(mfrow = c(1, 4))
plot(UMAP[,1],UMAP[,2], main = "Pax7", col = Pax7C)
plot(UMAP[,1],UMAP[,2],  main = "Myf5", col = Myf5C)
plot(UMAP[,1],UMAP[,2], main = "Myo", col = MyoDC)
plot(UMAP[,1],UMAP[,2],  main = "Myogenin", col = MyogeninC)

legend(1, 95, legend=c("UMAP Pax7", "UMAP Myf5","UMAP MyoD", "UMAP Myogenin"),
       col=c("red", "blue","orange","green"), lty=1:2, cex=0.8)