library(gplots) # needed for color palette

# load data
d <- read.csv('/Users/vagan/Desktop/FR-FCM-ZY3K_gated.csv', as.is=TRUE)

# subset day0 only
d0 <- d[d$sample == "Day0",]

# mycol <- bluered(100)
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

invisible(lapply(sort(unique(d$sample)), function(x)
{
    ss <- d[d$sample==x,]
    # shuffle order of cells to ensure no plotting bias (overlaid information)
    set.seed(123)
    ss <- ss[sample(rownames(ss), nrow(ss)),]
    plot(asinh(ss$CD104), asinh(ss$CD9),
         col=getColor(ss$IdU,d$IdU), pch=19, cex=0.5,
         ylim=c(3,10), xlim=c(0,10), main=paste(x,"IdU"))
    plot(asinh(ss$CD104), asinh(ss$CD9),
         col=getColor(ss$MyoD,d$MyoD), pch=19, cex=0.5,
         ylim=c(3,10), xlim=c(0,10), main=paste(x,"MyoD"))
}))
