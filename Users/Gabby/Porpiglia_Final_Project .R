library(gplots) # needed for color palette

# load data
d <- read.csv('../../Data/Porpiglia/FR-FCM-ZY3K_gated.csv', as.is=TRUE)
#########
#DAY 0 ##
########
# subset day0 only
d0 <- d[d$sample == "Day0",]

mycol <- bluered(100)

par(mfrow=c(2,3))
# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d0$CD104), asinh(d0$CD9), 
     col=mycol[100*d0$IdU/max(d$IdU)], pch=19)

#DAY 3 
# subset day3 only
d3 <- d[d$sample == "Day3",]

mycol <- bluered(100)


# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d3$CD104), asinh(d3$CD9), 
     col=mycol[100*d3$IdU/max(d$IdU)], pch=19)

#DAY6
# subset day6 only
d6 <- d[d$sample == "Day6",]

mycol <- bluered(100)


# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d6$CD104), asinh(d6$CD9), 
     col=mycol[100*d6$IdU/max(d$IdU)], pch=19)

#MYoD 
## plot points by CD104 and CD9 values colored by MyoD value scaled to max value in full data

## Day 0 
plot(asinh(d0$CD104), asinh(d0$CD9), 
     col=mycol[100*d0$MyoD/max(d$MyoD)], pch=19)

## Day 3 
# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d3$CD104), asinh(d3$CD9), 
     col=mycol[100*d3$MyoD/max(d$MyoD)], pch=19)

##Day 6 
# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d6$CD104), asinh(d6$CD9), 
     col=mycol[100*d6$MyoD/max(d$MyoD)], pch=19)
