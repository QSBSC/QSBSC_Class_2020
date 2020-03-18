library(gplots) # needed for color palette

# load data
d <- read.csv('../../Data/Porpiglia/FR-FCM-ZY3K_gated.csv', as.is=TRUE)

# subset day0 only
d0 <- d[d$sample == "Day0",]

mycol <- bluered(100)


# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d0$CD104), asinh(d0$CD9), 
     col=mycol[100*d0$IdU/max(d$IdU)], pch=19)

