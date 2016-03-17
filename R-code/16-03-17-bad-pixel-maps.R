library("IO.Pixels")

load.images(160314, "black")

hist(b.160314)

z <- unique(which(b.160314 > 9360), arr.ind = T)

bpm.img <- readTIFF("./Other-data/Other-images/BadPixelMap_160314/BadPixelMapBlack.tif", as.is = T)

fpath <- "./Other-data/Other-images/BadPixelMap_160314/"


# import bad pixel coordinates
bpm <- as.data.frame(t(as.matrix(as.data.frame(
    xmlToList(xmlParse(paste0(fpath, "BadPixelMap.bpm.xml")))$BadPixels))),
    stringsAsFactors = F)

for (i in 1:3) {
    bpm[,i] <- as.integer(bpm[,i])
}

pw.m.b <- pixelwise.mean(b.160314)
pixel.image(pw.m.b)
points(bpm$X, bpm$Y)

