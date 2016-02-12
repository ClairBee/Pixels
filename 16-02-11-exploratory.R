library(tiff)
library(plyr)

if (getwd() !=  "~Pixels") {setwd("~/Pixels")}

# Functions to create:
#  - import all images from single day & channel
#  - pixelwise mean
#  - pixelwise SD
#  - cropped histograms

##################################################################################
# Try using a tree to categorise pixels according to ?unspecified similarities
# May have to define manual algorithm, classifying by stepping through each sequence

##################################################################################
# array to contain pixelwise data
m <- array(dim=c(1996, 1996, 20))

load.daily <- function(img.date, channel, x = 1996, y = 1996, z = 20, progress = F, path = "./Image-data/") {
    
    img.date <- toString(img.date)
    rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
    m <- array(dim = c(x, y, z))
    
    # show progress bar or not?
    if (progress) {pb <- txtProgressBar(max = z, style = 3)}
    
    filenm <- paste(path,img.date,"/",channel,"/",channel,"_",rev.date,"_",i,".tif", sep = "")
    
    for (i in 1:z) {
        m[,,i] <- readTIFF(filenm)
        if (progress) {setTxtProgressBar(pb,i)}
    } 
    return(m)
}

system.time({
for (i in 1:20) {
    # need to rotate/flip matrix to correct orientation
    m[,,i] <- readTIFF(paste("./Image-data/150828/Black/black_280815_",i,".tif", sep = ""))
    setTxtProgressBar(pb, i)
}   # 4 elapsed
})

# get pixelwise mean & SD

system.time({
    PWmean <- apply(m, c(1,2), mean)        # 240 elapsed
    PWsd <- apply(m, c(1,2), sd)
})


globalmean <- mean(m)
globalSD <- sd(m)

# cropped histogram
hist(m[,,20], breaks = "FD", ylim = c(0,100))

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
getwd()
