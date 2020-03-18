library(gplots) # needed for color palette

# load data
d <- read.csv('../../Data/Porpiglia/FR-FCM-ZY3K_gated.csv', as.is=TRUE)
#########
#DAY 0 ##
########
# subset day0 only
d0 <- d[d$sample == "Day0",]

mycol <- colorpanel(100,"blue","yellow","red")

getColor <- function(val, maxval=max(val), col=mycol, scale=c("linear","asinh")[2])
{
        if(scale=="linear")
        {
                v <- val
                mv <- max(maxval)
                i <- ifelse(100*v/mv > 100, 100, round(100*v/mv,1))
        } else if (scale=="asinh")
        {
                v <- asinh(val)
                mv <- asinh(max(maxval))
                i <- ifelse(100*v/mv > 100, 100, round(100*v/mv,1))
        }
        col[i]
}

par(mfcol=c(2,3), oma=c(4,2,2,4))
# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d0$CD104), asinh(d0$CD9),
     col=getColor(d0$IdU,d$IdU), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste("Day0","IdU"))

plot(asinh(d0$CD104), asinh(d0$CD9),
     col=getColor(d0$MyoD,d$MyoD), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste("Day0","MyoD"))

#DAY 3 
# subset day3 only
d3 <- d[d$sample == "Day3",]

# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d3$CD104), asinh(d3$CD9),
     col=getColor(d3$IdU,d$IdU), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste("Day3","IdU"))

plot(asinh(d3$CD104), asinh(d3$CD9),
     col=getColor(d3$MyoD,d$MyoD), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste("Day3","MyoD"))

#DAY6
# subset day6 only
d6 <- d[d$sample == "Day6",]

# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d6$CD104), asinh(d6$CD9),
     col=getColor(d6$IdU,d$IdU), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste("Day6","IdU"))

plot(asinh(d6$CD104), asinh(d6$CD9),
     col=getColor(d$MyoD,d$MyoD), pch=19, cex=0.5,
     ylim=c(3,10), xlim=c(0,10), main=paste('Day6',"MyoD"))



