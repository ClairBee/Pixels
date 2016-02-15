
library("IO.Pixels")
if (getwd() !=  "~Pixels") {setwd("~/Pixels")}

d <- 1996
#####################################################################################################

# Add a function to import all channels for one day, as well as all days for one channel?

# function to split data out into individual panels (compare distribution across panels)
# come up with a colour scheme based on distance from mean (above or below) & assign using function

#####################################################################################################

# Import daily snapshot
b.150828 <- load.daily(150828, "black")

# pixelwise mean & SD of all 20 acquisitions
pw.mean.150828 <- pixelwise.mean(b.150828)
pw.sd.150828 <- pixelwise.sd(b.150828)

# global mean & SD across all pixels & acquisitions
mean.150828 <- mean(b.150828)
sd.150828 <- sd(b.150828)

# histograms of pixelwise means & SDs with delineation of 'extreme' values marked
{
    par(mfrow = c(2,2))
    hist(pw.mean.150828, breaks = "FD", main = "Histogram of pixelwise means")
    abline(v = mean.150828+(6*c(-sd.150828, sd.150828)), col = "red")
    
    hist(pw.mean.150828, breaks = "FD", ylim = c(0,100), main = "Cropped histogram showing detail")
    abline(v = mean.150828+(6*c(-sd.150828, sd.150828)), col = "red")
    
    hist(pw.sd.150828, breaks = "FD", main = "Histogram of pixelwise SDs")
    
    hist(pw.sd.150828, breaks = "FD", ylim = c(0,100), main = "Cropped histogram showing detail")
    par(mfrow = c(1,1))
}

# image plot
col.bands <- c(0, mean.150828 + (c(-3,-2,0,2,3) * sd.150828), (1 + mean.150828 + (3*sd.150828))/2, 1)
filled.contour(x = c(1:d), y = c(1:d), pw.mean.150828, levels = col.bands, col = topo.colors(10), 
               plot.axes = {points(ppp.high, pch = 1, col = "white", cex = 0.5); axis(1); axis(2)},
               main = "Mean black value: 15-08-28")

# try plotting points where mean value is > global mean + 6 sd
ind.high <- which(pw.mean.150828 == 1.0, arr.ind = T)
ppp.high <- ppp(ind.high[,1], ind.high[,2], c(1,d), c(1,d))
plot.ppp(ppp.high, pch=0, legend=FALSE, cols = "blue", main = "Points with mean value 1 ('hot pixels')")

