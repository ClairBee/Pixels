library(tiff)
library(plyr)

if (getwd() !=  "~Pixels") {setwd("~/Pixels")}

# Functions to create:
#  - import all images from single day & channel - UPDATE FUNCTION TO ROTATE MATRIX
#  - pixelwise mean
#  - pixelwise SD
#  - cropped histograms

##################################################################################
# Try using a tree to categorise pixels according to ?unspecified similarities
# May have to define manual algorithm, classifying by stepping through each sequence

# Label axes of daily matrix?

# should deal with fact that not always 20 snapshots - at some point, at least!
##################################################################################

# correct orientation of matrix data to match original display direction
z <- readTIFF("sample-tif.tiff")[,,1]
zz <- t(apply(z,2,rev))     # correct matrix to be correctly oriented
t(mat[nrow(mat):1,,drop=FALSE])     # alternative approach: need to test which is better

# difference in time elapsed is v small for 1996x1996 matrix (such as those we are importing)
system.time(image(t(apply(mean.150828,2,rev))))                           # 8.18 elapsed
system.time(image(t(mean.150828[nrow(mean.150828):1,,drop=FALSE])))       # 8.09 elapsed    

##################################################################################
# import daily snapshots
b.150828 <- load.daily(150828, "black")

# daily mean & SD
mean.150828 <- apply(b.150828, c(1,2), mean)
sd.150828 <- apply(b.150828, c(1,2), sd)

# histogram of daily mean & SD

# to plot bad pixels, use 'which' with arr.ind = T to extract indices of points to plot (Sept28-377)
lower.lim <- 0.05
upper.lim <- 0.807

ind.low <- which(mean.150828 < lower.lim, arr.ind = T)
ind.high <- which(mean.150828 > upper.lim, arr.ind = T)

nrow(ind.low); nrow(ind.high)       # no. of 'bad' pixels by value only

ppp.low <-ppp(ind.low[,2], ind.low[,1], c(1,d), c(1,d), marks=rep("l",length(ind.low[,2]))) 
ppp.high <-ppp(ind.high[,2], ind.high[,1], c(1,d), c(1,d), marks=rep("h",length(ind.high[,2])))

ppp.bad <- superimpose(ppp.low, ppp.high)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(0,0,0,0))
plot.ppp(ppp.low, pch=0, legend=FALSE, cols = "blue")
points(ppp.high, pch=1, col = "red")

# need to check orientation: is this the same way up as the original .tif?

# add K-function (need to read up on this!)
plot(envelope(ppp.low, Kest))
plot(envelope(ppp.high, Kest))

plot(density(ppp.low))
points(ppp.low, pch = 0, col = "white")
plot(density(ppp.high))
points(ppp.high, pch = 0, col = "white")

##################################################################################

# function to import multiple days' daily scans for a single channel
load.channel <- function(channel, x = 1996, y = 1996, z = 20, progress = F, fpath = "./Image-data/") {
    
    dates <- list.files(fpath)
    # prob. need to add some validation on folder names as well
    
    d <- length(dates)
    
    # add a 'from date and 'to date' rather than just picking everything up, modify d
    
    m <- array(dim = c(x,y,z,d))

    if (progress) {pb <- txtProgressBar(max = d, style = 3)}
    
    for (j in 1:d) {
        m[,,,j] <- load.daily(dates[j], "black")
        if (progress) {setTxtProgressBar(pb,j)}
    }
    
    return(m)
}

tmp <- load.channel("black", progress = T)


# get pixelwise mean & SD
daily.w.mean <- apply(daily.w, c(1,2), mean)
daily.w.sd <- apply(daily.w, c(1,2), sd)

# get pixelwise mean & SD



globalmean <- mean(m)
globalSD <- sd(m)

# cropped histogram
hist(mean.150828, breaks = "FD", ylim = c(0,100))

cropped.hist <- function(data, ymax = 100) {
    hist(data, breaks = "FD", ylim = c(0,ymax))
}

cropped.hist(m[,,1])

# try handling TIFF image as a raster
r <- readTIFF("./Image-data/150828/Black/black_280815_20.tif", native = T)

p <- as.ppp(m[,,20])

########################### PITIFUL ATTEMPT TO FIX GRAPHICS DEVICE ################
sessionInfo()
capabilities()
png(); dev.off()
x11()
Sys.getenv("DISPLAY")
###################################################################################

# try making larger matrix, giving multiple days

