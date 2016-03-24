
library("IO.Pixels")
pw.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
image <- pw.b[,,"160314"]

c <- 974; r <- c(1, 992)
transect <- pw.b[c, r[1]:r[2], "150828"]
sm <- lowess(transect, f = 1/15)$y
res <- transect - sm
abline(h = mad(res) * c(1,-1), col = "red")

o.plot(res)
illness <- res / mad(res)

which(illness > 2)
o.plot(transect)
hist(illness, breaks = "fd", xlim = c(-10,10))

points(cbind(which(illness > 2), transect[which(illness > 2)]), col = "purple")
points(cbind(which(illness > 10), transect[which(illness > 10)]), col = "red")

abline(v = c(992 - (1024 * c(0.25, 0.5, 0.75))), col = "blue")

abline(h = mean(transect), col = "green3")

###################################################################################################

bp <- get.bad.pixels(pw.b[,,"160314"])
bp <- get.bad.pixels(pw.b[,,"160314"], method = "mad")

bp <- get.bad.pixels(pw.b[,,"160314"], n.sigma = 6)
bp <- get.bad.pixels(pw.b[,,"160314"], n.sigma = 6, method = "mad")

bp <- get.bad.pixels(pw.b[,,"160314"], n.sigma = 10, method = "mad")
bp <- get.bad.pixels(pw.b[,,"160314"], n.sigma = 30)

points(which(pw.b[,,"160314"] == 65535, arr.ind = T), col = "gold", pch = 20)       # always on
points(which(pw.b[,,"160314"] == 0, arr.ind = T), col = "black", pch = 20)          # always off

# load 'official' bad pixel maps to compare
{
    fpath <- "./Other-data/Other-images/BadPixelMap_160314/"
    
    # import bad pixel coordinates
    bpm <- as.data.frame(t(as.matrix(as.data.frame(
        xmlToList(xmlParse(paste0(fpath, "BadPixelMap.bpm.xml")))$BadPixels))),
        stringsAsFactors = F)
    
    for (i in 1:3) bpm[,i] <- as.integer(bpm[,i])
    
    # Y-coordinates are measured from top axis: invert
    bpm$Y <- 1996 - bpm$Y
    
    bpm.params <- lapply(lapply(do.call("c", c(xmlToList(xmlParse(paste0(fpath, 
                                                                         "CalibrationParameters.xml"))),
                                               list(recursive = TRUE))), FUN = unlist), as.integer)
}

points(bpm$X, bpm$Y, col = "red")
summary(bpm$Y)

b.dead <- which(pw.b[,,"160314"] == 0, arr.ind = T)
b.dim <- which((pw.b[,,"160314"] < bpm.params$BlackMinThreshold) & (pw.b[,,"160314"] > 0), arr.ind = T)
b.warm <- which((pw.b[,,"160314"] > bpm.params$BlackMaxThreshold) & (pw.b[,,"160314"] < 65535), arr.ind = T)
b.hot <- which(pw.b[,,"160314"] == 65535, arr.ind = T)

plot(bpm$X, bpm$Y, col = "grey", pch = 16, "Bad pixels per calibration file")
draw.panels()

points(b.dead, pch = 20, col = "black", xlim = c(1,1996), ylim = c(1,1996))
points(b.dim, pch = 20, col = "blue")
points(b.warm, pch = 20, col = "orange")
points(b.hot, pch = 20, col = "red")

hist(pw.b[,,"160314"], breaks = "fd", xlim = c(3000, 10000))
abline(v = c(bpm.params$BlackMinThreshold, bpm.params$BlackMaxThreshold), col = "red")

hist(pw.b[,,"160314"], breaks = "fd", xlim = c(3000, 10000), ylim = c(0,100))
abline(v = c(bpm.params$BlackMinThreshold, bpm.params$BlackMaxThreshold), col = "red")
