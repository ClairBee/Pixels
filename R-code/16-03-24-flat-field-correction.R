
# flat-field correction

pw.w <- readRDS("./Other-data/Pixelwise-means-white.rds")
pw.g <- readRDS("./Other-data/Pixelwise-means-grey.rds")
pw.b <- readRDS("./Other-data/Pixelwise-means-black.rds")

# choose an image with no spots: 150828
R <- pw.g[,,"150828"]
D <- pw.b[,,"150828"]
FF <- pw.w[,,"150828"]
m <- mean(FF - D)

corr <- m * (R - D) / (FF - D)
corr[is.na(corr)] <- 0      # otherwise get NA where FF == D

pixel.image(corr)
s.hist(corr)
abline(v = mean(corr), col = "red")
abline(v = mean(corr) + 1.96 * sd(corr), col = "red")

bp <- which(abs(corr - mean(corr)) > (1.96 * sd(corr)), arr.ind = T)
points(bp)

# per-panel linear correction to account for offset & readout dark
panel.lm <- fit.panel.lm(corr)
pixel.image(panel.lm$fitted.values)
panel.res <- corr - panel.lm$fitted.values
pixel.image(panel.res)

hist(panel.res, breaks = "fd", prob = T, xlim = c(-200, 200))
lines(c(-400:400), dnorm(c(-400:400), mean = mean(panel.res), sd = sd(panel.res)), lwd = 3, col = "cornflowerblue")
lines(c(-400:400), dnorm(c(-400:400), mean = mean(panel.res), sd = mad(panel.res)), lwd = 3, col = "red")
# normalise using mad, not SD - too susceptible to outliers

sickness <- panel.res / mad(panel.res)
image(c(1:1996), c(1:1996), sickness, 
      breaks = c(floor(min(sickness)), -5, -2, 2, 5, ceiling(max(sickness))), 
      col = c("black", "blue", "white", "blue", "black"))
hist(sickness, breaks = "fd", xlim = c(-20, 20), col = "black", ylim = c(0,20))
hist(sickness[abs(sickness) > 4], breaks = "fd", add = T, col = "red")

hot <- which(R == 65535, arr.ind = T)
warm <- which(sickness > 4, arr.ind = T)
cool <- which(sickness < -4, arr.ind = T)
dead <- which(R == 0, arr.ind = T)
    
bp <- rbind(cbind(which(R == 65535, arr.ind = T), type = 1),
            cbind(which(R == 0, arr.ind = T), 4),
            cbind(which(sickness > 4, arr.ind = T), 2),
            cbind(which(sickness < -4, arr.ind = T), 3))
bp <- bp[!duplicated(bp[,1:2]),]

plot(bp, pch = 20, col = c("greenyellow", "red", "blue", "black")[bp[,3]])

###################################################################################################
cols <- c("blue", "skyblue", "greenyellow", "yellow", "gold", "orange",
          "red3", "violetred", "purple", "black")

# find saturated points in each channel
mx.b <- which(apply(pw.b,c(1,2), max) == 65535, arr.ind = T)
mx.g <- which(apply(pw.g,c(1,2), max) == 65535, arr.ind = T)
mx.w <- which(apply(pw.w,c(1,2), max) == 65535, arr.ind = T)

mx <- rbind(mx.b, mx.g, mx.w)
mx <- mx[!duplicated(mx),]
mx <- mx[order(mx[,2], mx[,1]),]

# get all values of those points
hot.b <- apply(pw.b, 3, "[", mx)
hot.g <- apply(pw.g, 3, "[", mx)
hot.w <- apply(pw.w, 3, "[", mx)

# filter out any points that are always at 65535 - not worth plotting
hot.b.moved <- hot.b[which(apply(hot.b, 1, min) != 65535),]
hot.g.moved <- hot.g[which(apply(hot.g, 1, min) != 65535),]
hot.w.moved <- hot.w[which(apply(hot.w, 1, min) != 65535),]

# plot progression of pixels in black channel
{
plot(hot.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
     xlab = "", ylab = "", main = "Progression of black pixels")
for (i in 2:11) {
    points(hot.b[,i], pch = 20, col = adjustcolor("blue", alpha = i/11))
}
abline(h = mean(pw.b) + (1.96*sd(pw.b)), col = "blue")
abline(h = mean(pw.b), lwd = 2, col = "blue")
}

# plot progression of pixels in grey channel
{
plot(hot.g[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("green3", alpha = 1/11),
     xlab = "", ylab = "", main = "Progression of grey pixels")
for (i in 2:11) {
    points(hot.g[,i], pch = 20, col = adjustcolor("green3", alpha = i/11))
}
abline(h = mean(pw.g) + (1.96*sd(pw.g)), col = "green3")
abline(h = mean(pw.g), lwd = 2, col = "green3")
}

# plot progression of pixels in white channel
{
plot(hot.w[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11),
     xlab = "", ylab = "", main = "Progression of white pixels")
for (i in 2:11) {
    points(hot.w[,i], pch = 20, col = adjustcolor("red", alpha = i/11))
}
abline(h = mean(pw.w) + (1.96*sd(pw.w)), col = "red")
abline(h = mean(pw.w), lwd = 2, col = "red")
}

# try normalising black pixel values by subtracting per-image mean value
m.black <- matrix(rep(apply(pw.b, 3, mean), nrow(hot.b)), nrow = nrow(hot.b), byrow = T)
hot.b.adj <- hot.b - (m.black - min(m.black))
hot.b.adj[hot.b == 65535] <- 65535  # reset absolute values

# replot
{
    plot(hot.b.adj[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of black pixels")
    for (i in 2:11) {
        points(hot.b.adj[,i], pch = 20, col = adjustcolor("blue", alpha = i/11))
    }
}

# same for grey & white images
m.grey <- matrix(rep(apply(pw.g, 3, mean), nrow(hot.g)), nrow = nrow(hot.g), byrow = T)
hot.g.adj <- hot.g - (m.grey - min(m.grey))
hot.g.adj[hot.g == 65535] <- 65535  # reset absolute values

m.white <- matrix(rep(apply(pw.w, 3, mean), nrow(hot.w)), nrow = nrow(hot.w), byrow = T)
hot.w.adj <- hot.w - (m.white - min(m.white))
hot.w.adj[hot.w == 65535] <- 65535  # reset absolute values

{
    plot(hot.g.adj[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("green3", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of grey pixels")
    for (i in 2:11) {
        points(hot.g.adj[,i], pch = 20, col = adjustcolor("green3", alpha = i/11))
    }
}

{
    plot(hot.w.adj[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of white pixels")
    for (i in 2:11) {
        points(hot.w.adj[,i], pch = 20, col = adjustcolor("red", alpha = i/11))
    }
}

###################################################################################################
# 'OFFICIAL' BAD PIXEL MAP

# load bad pixel map provided by Jay to compare
fpath <- "./Other-data/Other-images/BadPixelMap_160314/"

# import bad pixel coordinates
bpm <- as.data.frame(t(as.matrix(as.data.frame(
    xmlToList(xmlParse(paste0(fpath, "BadPixelMap.bpm.xml")))$BadPixels))),
    stringsAsFactors = F)

for (i in 1:3) bpm[,i] <- as.integer(bpm[,i])

bpm$Y <- 1996 - bpm$Y       # Y-coordinates are measured from top axis: invert
bpm$X <- bpm$X + 1          # X-coordinates are measured from 0: add 1

write.csv(bpm, "./Other-data/BadPixelMap-160314.csv", row.names = F)
bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)

bpm.b <- apply(pw.b, 3, "[", as.matrix(bpm[,c("X", "Y")]))
bpm.g <- apply(pw.g, 3, "[", as.matrix(bpm[,c("X", "Y")]))
bpm.w <- apply(pw.w, 3, "[", as.matrix(bpm[,c("X", "Y")]))

# plot in black channel
{
    plot(bpm.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of black pixels")
    for (i in 2:11) {
        points(bpm.b[,i], pch = 20, col = adjustcolor("blue", alpha = i/11))
    }
}

# plot in grey channel
{
    plot(bpm.g[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("green3", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of grey pixels")
    for (i in 2:11) {
        points(bpm.g[,i], pch = 20, col = adjustcolor("green3", alpha = i/11))
    }
}

# plot in white channel
{
    plot(bpm.w[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11),
         xlab = "", ylab = "", main = "Progression of white pixels")
    for (i in 2:11) {
        points(bpm.w[,i], pch = 20, col = adjustcolor("red", alpha = i/11))
    }
}

# does overplotting make any sense?
plot(bpm.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
     xlab = "", ylab = "", main = "Progression of pixels", xlim = c(500,1000))
for (i in 2:11) {
    points(bpm.b[,i], pch = 20, col = adjustcolor("blue", alpha = 1/11))
}
for (i in 1:11) {
    points(bpm.g[,i], pch = 20, col = adjustcolor("green3", alpha = 1/11))
}
for (i in 1:11) {
    points(bpm.w[,i], pch = 20, col = adjustcolor("red", alpha = 1/11))
}

###################################################################################################
# plot non-bad pixels: what should 'normal' behaviour be?
all <- merge(merge(x = c(1:1996), y = c(1:1996)), bpm, by.x = c(1:2), by.y = c(1:2), all.x = T)
good.pix <- all[is.na(all$Mask),]

samp <- good.pix[sample(1:nrow(good.pix), nrow(bpm), replace = F),]

gpm.b <- apply(pw.b, 3, "[", as.matrix(samp[, c("x", "y")]))
gpm.g <- apply(pw.g, 3, "[", as.matrix(samp[, c("x", "y")]))
gpm.w <- apply(pw.w, 3, "[", as.matrix(samp[, c("x", "y")]))

# plot sample of 'normal' pixels (plotting all pixels takes too long)
    plot(gpm.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
         xlab = "", ylab = "", main = "Sample of 'healthy' pixels")
    points(gpm.g[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("green3", alpha = 1/11))
    points(gpm.w[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("gold", alpha = 1/11))

# overlay sample of 'normal' pixels with 'bad' pixels
    plot(gpm.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("blue", alpha = 1/11),
         xlab = "", ylab = "", main = "Healthy vs 'bad' pixels - black images")
    points(bpm.b[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11))
    
    plot(gpm.g[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("green3", alpha = 1/11),
         xlab = "", ylab = "", main = "Healthy vs 'bad' pixels - grey images")
    points(bpm.g[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11))
    
    plot(gpm.w[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("gold", alpha = 1/11),
         xlab = "", ylab = "", main = "Healthy vs 'bad' pixels - white images")
    points(bpm.w[,1], pch = 20, ylim = c(0, 65535), col = adjustcolor("red", alpha = 1/11))
    
###################################################################################################

o.plot(bpm.g[550,], ylim = c(0,65535), col = adjustcolor("green3", alpha = 0.3))

for (i in 0:30) {
    o.plot(bpm.g[550 + i,], add = T, col = adjustcolor("green3", alpha = 0.3))
}



###################################################################################################
# distinct type of 'badness' on LHS and RHS.

plot(as.matrix(bpm[,c("X", "Y")]))
points(as.matrix(bpm[1:400,c("X", "Y")]), pch = 20, col = "red")
points(as.matrix(bpm[2500:2842,c("X", "Y")]), pch = 20, col = "red")
points(as.matrix(bpm[500:900,c("X", "Y")]), pch = 20, col = "red")
