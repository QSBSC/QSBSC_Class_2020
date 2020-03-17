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
     col=mycol[100*d0$IdU/max(d$IdU)], pch=19)

#DAY 3 
# subset day3 only
d3 <- d[d$sample == "Day3",]

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

# plot points by CD104 and CD9 values colored by IdU value scaled to max value in full data
plot(asinh(d3$CD104), asinh(d3$CD9), 
     col=mycol[100*d3$IdU/max(d$IdU)], pch=19)

#DAY6
# subset day6 only
d6 <- d[d$sample == "Day6",]

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
