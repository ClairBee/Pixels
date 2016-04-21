
# identifying blobs on screen based on morphology

library("IO.Pixels")

library(mmand)      # needed for morphological operations

fpath <- "./Notes/Spots-on-screen/fig/"
load.pixel.means()
ff <- readRDS("./Other-data/Flat-field-corrected.rds")

####################################################################################################
# spots on beryllium screen manifest as blobs of dim pixels in white & grey images
# exclude these clusters when fitting distribution for bright/dim thresholding
# still look for bright pixels within these blobs - (individual pixel response is too high)
# not sure whether to look for dim pixels within these blobs. Maybe exclude.

####################################################################################################
# WHITE IMAGE ONLY                                                                              ####

# we know there's a panelwise gradient, so fit this to remove panel edge effects
panels <- panel.lm(pw.m[,,"white", "160314"], "x + y")
panel.res <- pw.m[,,"white", "160314"] - panels$fitted.values
pixel.image(panel.res)

# Gaussian blur of panel residuals
{
    pixel.image(blur(as.im(panel.res), sigma = 1), title = "Sigma: 1")
    pixel.image(blur(as.im(panel.res), sigma = 2), title = "Sigma: 2")
    blur.20 <- blur(as.im(panel.res), sigma = 20)
    pixel.image(blur.20, title = "Sigma: 20")
    
    im <- matrix(blur(as.im(panel.res), sigma = 1)$v, ncol = 1996)
    pixel.image(im)
    
    o.plot(panel.res[50,1:992])
    o.plot(matrix(blur.20$v, ncol = 1996)[50,1:992], col = "blue", add = T)
    o.plot(lowess(panel.res[50,1:992], f = 1/5), col = "magenta3", add = T)
}

# looks like lowess might actually be a better way to go than convolution...
# (unnecessary to fit panels first)
smoo <- lowess.per.column(pw.m[,,"white", "160314"], span = 1/5)

loess.res <- pw.m[,,"white", "160314"] - smoo
pixel.image(loess.res)
s.hist(loess.res)

o.plot(loess.res[1000,])
abline(h = -sd(abs(loess.res)), col = "red")

# TEST CASES
pixel.image(loess.res, xlim = c(770, 860), ylim = c(850, 930))
s.hist(loess.res[770:860, 850:930])

# small cluster: c. 20px across, 15 high
pixel.image(loess.res, xlim = c(1400, 1450), ylim = c(900, 950))
s.hist(loess.res[1400:1450, 900:950])

th.large <- threshold(loess.res[770:860, 850:930], method = "kmeans")
pixel.image(th.large)
th.small <- threshold(loess.res[1400:1450, 900:950], method = "kmeans")
pixel.image(th.small)

# works fine on small areas, not over full image spread
th <- threshold(loess.res, method = "kmeans")
pixel.image(th)

# need to perform opening/closing first
display(loess.res)

# kernel size 5 circular should do it
eroded <- erode(loess.res, shapeKernel(c(5, 5), type = "disc"))
pixel.image(loess.res)
pixel.image(eroded, title = "Eroded Loess residuals")
th.eroded <- threshold(eroded, method = "kmeans")
pixel.image(th.eroded)

dilated <- dilate(loess.res, shapeKernel(c(5, 5), type = "disc"))
pixel.image(dilated, title = "Dilated Loess residuals")
th.dilated <- threshold(dilated, method = "kmeans")
pixel.image(th.dilated)
hist(dilated, breaks = "fd", prob = T)
lines(c(-6000: 26000), dnorm(c(-6000: 26000), mean(dilated), mad(dilated)), col = "cyan3", lwd = 2)

edges <- dilated - eroded
pixel.image(edges)

display(edges)
# not particularly helpful.

mf <- medianFilter(loess.res, shapeKernel(c(5, 5), type = "disc"))
pixel.image(mf, title = "Median-filtered Loess residuals")
gf <- gaussianSmooth(loess.res, c(5, 5))
pixel.image(gf, title = "Gaussian smoothed Loess residuals")

erosion.gf <- gaussianSmooth(eroded, c(5,5))
pixel.image(erosion.gf, title = "Gaussian-smoothed erosion of Loess residuals")
# Gaussian smoothing just enlarges the small clusters of dim pixels. Not helpful. 

# Gaussian convolution over residuals?
{
    blur.1 <- blur(as.im(loess.res), sigma = 1)
    blur.2 <- blur(as.im(loess.res), sigma = 2)
    blur.5 <- blur(as.im(loess.res), sigma = 5)
    blur.10 <- blur(as.im(loess.res), sigma = 10)
    
    
    o.plot(loess.res[1000,])
    o.plot(matrix(blur.1$v, ncol = 1996)[1000,], add = T, col = adjustcolor("blue", alpha = 0.4))
    o.plot(matrix(blur.2$v, ncol = 1996)[1000,], add = T, col = adjustcolor("green3", alpha = 0.4))
    o.plot(matrix(blur.5$v, ncol = 1996)[1000,], add = T, col = adjustcolor("red", alpha = 0.4))
    o.plot(matrix(blur.10$v, ncol = 1996)[1000,], add = T, col = adjustcolor("gold", alpha = 0.4))
    
    # still requires threshold setting
}

# or just move straight to morph. closing?

# TRUNCATE DATA                                                                                 ####

m.mar <- c(3,2,1,1)
# raw data plot
{
    jpeg(paste0(fpath, "white-image.jpg"))
        par(mar = m.mar)
        pixel.image(pw.m[,,"white", "160314"])
    dev.off()
}

# apply loess smoothing
smoo <- lowess.per.column(pw.m[,,"white", "160314"], span = 1/5)
{
    jpeg(paste0(fpath, "loess-smoothing.jpg"))
        par(mar = m.mar)
        pixel.image(smoo, break.levels = sd.levels(pw.m[,,"white", "160314"]))
    dev.off()
}

loess.res <- pw.m[,,"white", "160314"] - smoo
{
    jpeg(paste0(fpath, "loess-residuals.jpg"))
        par(mar = m.mar)
        pixel.image(loess.res)
    dev.off()
}

# save sd levels for clearer plots
org.breaks <- sd.levels(loess.res)

# set all values > mean (since we're only interested in dim pixels)

trunc <- loess.res
trunc[trunc > mean(trunc)] <- mean(trunc)
{
    jpeg(paste0(fpath, "loess-truncated.jpg"))
        par(mar = m.mar)
        pixel.image(trunc, break.levels = org.breaks)
    dev.off()
}

trunc.eroded <- erode(trunc, shapeKernel(c(5, 5), type = "disc"))
pixel.image(trunc.eroded, break.levels = org.breaks, title = "Erosion of truncated Loess residuals")
# wrong way

trunc.dilated <- dilate(trunc, shapeKernel(c(5, 5), type = "disc"))
{
    jpeg(paste0(fpath, "dilated.jpg"))
        par(mar = m.mar)
        pixel.image(trunc.dilated, break.levels = org.breaks)
    dev.off()
}

th <- threshold(trunc.dilated, method = "kmeans")
pixel.image(th, title = "Thresholded dilation of truncated Loess residuals")

# perform erosion to return shapes to approx. original size
blobs <- erode(trunc.dilated, shapeKernel(c(5, 5), type = "disc"))
{
    jpeg(paste0(fpath, "eroded.jpg"))
        par(mar = m.mar)
        pixel.image(blobs, break.levels = org.breaks)
    dev.off()
}
pixel.image(blobs, break.levels = org.breaks, title = "blobs identified")

blobs.th <- threshold(blobs, method = "kmeans")
pixel.image(blobs.th, title = "Final classification into dim spots/normal regions")

table(c(blobs.th))
# convert this to a raster to identify clumps?