
library("IO.Pixels"); library("spdep")

# replicating 'official' classification of underperforming pixels
# establish hierarchy so that each bad pixel has only one classification?

# check which outputs are still integer-valued and specify data type accordingly
# (converted to numeric on arithmetic operation)

########################################################################################################
# load day's data for all channels

z <- load.daily(150828)

bad.pixels <- list()

# bright offset corrected image: X_i - Y_b
offset.img <- pixelwise.mean(z[,,,1])
gain.img <- pixelwise.mean(z[,,,3])
bright.offset.corrected <- sweep(z[,,,3], c(1,2), offset.img, "-")

# underperforming bright pixel: value > 1.5*median
median.bright <- median(bright.offset.corrected)
unique(which(b.150828 < ll.b.150828, arr.ind = T)[,c(1:2)])

bad.pixels[[1]] <- unique(which(bright.offset.corrected > (1.5 * median.bright), arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[1]])


# no gain: dark pixel, no light response (pixelwise mean value == 0)
bad.pixels[[2]] <- unique(which(gain.img == 0, arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[2]])

ppp.off <- ppp(bad.pixels[[2]][,2], bad.pixels[[2]][,1], c(1,1996), c(1,1996))
plot.ppp(ppp.off, pch=0, legend=FALSE, cols = "blue")


# underperforming dark pixel: value < 0.45*median
bad.pixels[[3]] <- unique(which(bright.offset.corrected < (0.45 * median.bright), arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[3]])

ppp.dark <- ppp(bad.pixels[[3]][,2], bad.pixels[[3]][,1], c(1,1996), c(1,1996))
points(ppp.dark, pch = 1, col = "red")


# bright noise: in white channel, pixelwise sigma > 6*median sigma
pw.sd.bright <- pixelwise.sd(z[,,,3])
med.sd.bright <- median(pw.sd.bright)

bad.pixels[[4]] <- unique(which(pw.sd.bright > (6 * med.sd.bright), arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[4]])

ppp.bright.noise <- ppp(bad.pixels[[4]][,2], bad.pixels[[4]][,1], c(1,1996), c(1,1996))
points(ppp.bright.noise, pch = 16, col = "green")


# dark noise: in black channel, pixelwise sigma > 6*median sigma
pw.sd.dark <- pixelwise.sd(z[,,,1])
med.sd.dark <- median(pw.sd.dark)

bad.pixels[[5]] <- unique(which(pw.sd.dark > (6 * med.sd.dark), arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[5]])

ppp.dark.noise <- ppp(bad.pixels[[5]][,2], bad.pixels[[5]][,1], c(1,1996), c(1,1996))
points(ppp.dark.noise, pch = 16, col = "black")

# offset & gain corrected image: (X_i - Y_b) / (Y_w - Y_b) where X_i is grey image
grey.offset <- sweep(z[,,,2], c(1,2), offset.img, "-")
gain.offset <- gain.img - offset.img
grey.offset.gain.correction <- sweep(grey.offset, c(1,2), gain.offset, "/")
# grey.offset.gain.correction[which(is.na(grey.offset.gain.correction))] <- 0


# globally non-uniform pixel: pixel value more than +- 2% from median in offset correction image
med.grey <- median(grey.offset.gain.correction, na.rm = T)
bad.pixels[[6]] <- unique(which(grey.offset.gain.correction > (1.02 * med.grey), arr.ind = T)[,c(1:2)])
nrow(bad.pixels[[6]])

# locally non-uniform pixel: pixel value more than +- 1% of median of 9x9 neighbours
tmp <- z[, , 1, 1]
local.med <- as.matrix(focal(raster(tmp), w=matrix(c(rep(1, 40),NA,rep(1, 40)), nrow = 9), fun = median, na.rm = T))

local.upper <- local.med*1.01
local.lower <- local.med*0.99

local.high <- which(tmp > local.upper, arr.ind = T)[, c(1:2)]
local.low <- which(tmp < local.lower, arr.ind = T)[, c(1:2)]
nrow(local.high)
nrow(local.low)

########################################################################################################

# once pixels have been classified by behaviour, classify by shape (lines/clusters etc)