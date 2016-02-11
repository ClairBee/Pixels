library(tiff)
library(plyr)

# Functions to create:
#  - import all images from single day & channel
#  - pixelwise mean
#  - pixelwise SD


##################################################################################
# array to contain pixelwise data
m <- array(dim=c(1996, 1996, 20))

# import images for single day
pb <- txtProgressBar(max=20, style = 3)

for (i in 1:20) {
    m[,,i] <- readTIFF(paste("./Image-data/141009-black/black_091014_",i,".tif", sep = ""))
    setTxtProgressBar(pb, i)
}


# get pixelwise mean & SD
PWmean <- array(dim=c(1996,1996))
PWsd <- array(dim=c(1996,1996))

