
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

blobs.th <- array(1, dim = c(1996, 1996)) - threshold(blobs, method = "kmeans")
pixel.image(blobs.th, title = "Final classification into dim spots/normal regions")

{
    jpeg(paste0(fpath, "thresholded.jpg"))
        par(mar = m.mar)
        image(blobs.th, asp = T, col = c("white", "blue"))
    dev.off()
}

table(c(blobs))

hist(blobs, breaks =  c(floor(min(blobs)):ceiling(max(blobs))), ylim = c(0,1000))
hist(blobs[blobs.th == 1], breaks = c(floor(min(blobs)):ceiling(max(blobs))), add = T, 
     col = "blue", border = "blue", ylim = c(0,1000))

hist(blobs[blobs.th == 0], breaks = c(floor(min(blobs)):ceiling(max(blobs))),
     add = T, col = "red", border = "red", ylim = c(0,1000))

# convert this to a raster to identify clumps?

# FLAT-FIELD CORRECTED IMAGE                                                                    ####

ff <- (pw.m[,,"grey", "160314"] - pw.m[,,"black", "160314"]) / (pw.m[,,"white", "160314"] - pw.m[,,"black", "160314"])
ff[is.na(ff)] <- 0

pixel.image(ff)
br <- sd.levels(ff)

# truncate
ff[ff > mean(ff)] <- mean(ff)

ff.d <- dilate(ff, shapeKernel(c(5, 5), type = "disc"))
ff.e <- erode(ff.d, shapeKernel(c(5, 5), type = "disc"))

pixel.image(ff.e, break.levels = br)

# resulting shapes are much  smaller. Flat-field correction has removed the less extreme values.

# COMPARE JOHNSON DISTRIBUTION WITH & WITHOUT DIM PATCHES                                       ####

zz <- screen.spots(160314)

# summarise spots
sc <- ddply(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))),
                       id = getValues(zz)[!is.na(getValues(zz))]),
            .(id), summarise, 
            xm = round(mean(x),0), ym = round(mean(y),0),
            xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y),
            count = length(x))

sum(sc$count) / (1996 * 1996) * 100

# fit parametric model to white image
spot <- spot.lm(pw.m[,,"white", "160314"], o = 2)
spot.res <- matrix(spot$residuals, ncol = 1996)
panels <- panel.lm(spot.res, "x + y")
res <- spot.res - panels$fitted.values

# plot raw residuals with Johnson fit
hist(res, breaks = "fd", prob = T, ylim = c(0, 2e-06), xlim = c(-10000,10000), main = "Johnson fit to all residuals")
lines(c(-50000:22900), dJohnson(c(-50000:22900), JohnsonFit(res)), col = "red", lwd = 2)

# plot residuals + Johnson fit without dim spots
res.adj <- res[xyFromCell(zz, which(is.na(getValues(zz))))]
hist(res.adj, breaks = "fd", prob = T, ylim = c(0, 2e-06), xlim = c(-10000,10000), main = "Johnson fit without dim spots")
lines(c(-50000:22900), dJohnson(c(-50000:22900), JohnsonFit(res.adj)), col = "red", lwd = 2)

# plot residuals + Johnson fit without dim spots, without edge boundary
res.adj.c <- res[xyFromCell(crop(zz, c(10.5, 1985.5, 10.5, 1985.5)), 
                            which(is.na(getValues(crop(zz, c(10.5, 1985.5, 10.5, 1985.5))))))]
hist(res.adj.c, breaks = "fd", prob = T, ylim = c(0, 2e-06), xlim = c(-10000,10000), main = "Johnson fit without dim spots")
lines(c(-50000:22900), dJohnson(c(-50000:22900), JohnsonFit(res.adj.c)), col = "red", lwd = 2)


# q-q plots to compare
{
    plot(qJohnson(c(1:999)/1000, JohnsonFit(res)), 
         quantile(res, c(1:999)/1000),  col = adjustcolor("green3", alpha = 0.4),
         pch = 20, cex = 0.7, asp = T, ylab = "Observed quantile", xlab = "Johnson quantile")
    points(qJohnson(c(1:999)/1000, JohnsonFit(res.adj)), 
           quantile(res.adj, c(1:999)/1000), pch = 20, cex = 0.7, col = adjustcolor("blue", alpha = 0.4))
    points(qJohnson(c(1:999)/1000, JohnsonFit(res.adj.c)), 
           quantile(res.adj.c, c(1:999)/1000), pch = 20, cex = 0.7, col = adjustcolor("red", alpha = 0.4))
    legend("topleft", pch = 20, cex = 0.8, col = adjustcolor(c("green3", "blue", "red"), alpha = 0.4),
           legend = c("All residuals", "Dim spots removed", "Dim spots removed, cropped"), bty = "n")
    title("Q-Q plot of Johnson distribution against observed data")
    
    abline(0, 1, col = "red", lty = 2)
    abline(h = quantile(res, c(0.01, 0.99)), col = "skyblue", lty = 2)
    abline(v = qJohnson(c(0.01, 0.99), JohnsonFit(res)), col = adjustcolor("green3", alpha = 0.4), lty = 2)
    abline(v = qJohnson(c(0.01, 0.99), JohnsonFit(res.adj)), col = adjustcolor("blue", alpha = 0.4), lty = 2)
    abline(v = qJohnson(c(0.01, 0.99), JohnsonFit(res.adj.c)), col = adjustcolor("red", alpha = 0.4), lty = 2)

}

# Kullback-Leibler divergence to compare distances
library(FNN)

mean(KL.divergence(res, rJohnson(10000, JohnsonFit(res))))
mean(KL.divergence(sample(res, 10000), rJohnson(10000, JohnsonFit(res.adj))))


# Kolmogorov-Smirnov distances
ks.test(res, pJohnson, JohnsonFit(res))                 # D = 0.0059649
ks.test(res.adj, pJohnson, JohnsonFit(res.adj))         # D = 0.005435
ks.test(res.adj.c, pJohnson, JohnsonFit(res.adj.c))     # D = 0.0032077

ks.test(res, pJohnson, JohnsonFit(res.adj))             # D = 0.007422
ks.test(res, pJohnson, JohnsonFit(res.adj.c))           # D = 0.014879



# use MSE to assess fit?

# TRACK DEVELOPMENT OF DIM PATCHES                                                              ####
# takes around 4.5 minutes to run over all 11 images
spots <- lapply(dimnames(pw.m)[[4]], screen.spots)

cols <- c("orangered2", "orange", "gold", "green", "chartreuse3", "cyan3", "black", "black",
          "skyblue", "blue", "magenta2")

pdf(paste0(fpath, "All-screen-spots.pdf"))
image(c(1:1996), c(1:1996), t(as.matrix(spots[[1]], ncol = 1996)[1996:1,]),
      col = adjustcolor(cols[1], alpha = 0.5), asp = T, xlab = "", ylab = "")
for (i in 2:11) {
    image(c(1:1996), c(1:1996), t(as.matrix(spots[[i]], ncol = 1996)[1996:1,]), 
          col = adjustcolor(cols[i], alpha = 0.5), add = T)
}
dev.off()

# split images according to screen replacement
# Jan-26-2016 and previous
pdf(paste0(fpath, "Screen-spots-pre-Jan.pdf"))
image(c(1:1996), c(1:1996), t(as.matrix(spots[[1]], ncol = 1996)[1996:1,]),
      col = adjustcolor(cols[1], alpha = 0.5), asp = T, xlab = "", ylab = "")
for (i in 2:6) {
image(c(1:1996), c(1:1996), t(as.matrix(spots[[i]], ncol = 1996)[1996:1,]), 
          col = adjustcolor(cols[i], alpha = 0.5), add = T)
}
dev.off()

pdf(paste0(fpath, "Screen-spots-post-Jan.pdf"))
image(c(1:1996), c(1:1996), t(as.matrix(spots[[9]], ncol = 1996)[1996:1,]),
      col = adjustcolor(cols[9], alpha = 0.9), asp = T, xlab = "", ylab = "")
for (i in 10:11) {
    image(c(1:1996), c(1:1996), t(as.matrix(spots[[i]], ncol = 1996)[1996:1,]), 
          col = adjustcolor(cols[i], alpha = 0.5), add = T)
}
dev.off()

sc <- ddply(data.frame(xyFromCell(spots[[i]], which(!is.na(getValues(spots[[i]])))),
                       id = getValues(spots[[i]])[!is.na(getValues(spots[[i]]))]),
            .(id), summarise, 
            xm = round(mean(x),0), ym = round(mean(y),0),
            xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y),
            count = length(x))

zz <- screen.spots("150529")
qq <- t(as.matrix(zz, ncol = 1996)[1996:1,])
image(qq, col = "magenta3")

# RANDOM BAD SPOT NOTICED IN WHITE IMAGES                                                       ####
focus <- matrix(c(1105 + rep(c(-3:3), 7), 490 + sort(rep(c(-3:3), 7))), ncol = 2)
plot(0, type = "n", xlim = range(focus[,1]), ylim = range(focus[,2]), xlab = "", ylab = "", main = "white images")
text(focus, cex = 0.7, labels = round(pw.m[,,"white", "150828"][focus]/1000,0))

plot(0, type = "n", xlim = range(focus[,1]), ylim = range(focus[,2]), xlab = "", ylab = "", main = "grey images")
text(focus, cex = 0.7, labels = round(pw.m[,,"black", "150828"][focus]/1000,0))

plot(0, type = "n", xlim = range(focus[,1]), ylim = range(focus[,2]), xlab = "", ylab = "", main = "black images")
text(focus, cex = 0.7, labels = round(pw.m[,,"black", "150828"][focus]/1000,0))

# cluster of bright & 'no response' pixels?

# get threshold for normal behaviour in black images
bn <- qJohnson(0.999, JohnsonFit(pw.m[,,"black", "150828"]))

# identify pixels within that range in grey/white images
un.g <- which(matrix(findInterval(pw.m[,,"grey", "150828"], bn), ncol = 1996) == 0, arr.ind = T)
un.w <- which(matrix(findInterval(pw.m[,,"white", "150828"], bn), ncol = 1996) == 0, arr.ind = T)

points(un.g, pch = 16, col = adjustcolor("purple", alpha = 0.5), cex = 3)
points(un.w, pch = 16, col = adjustcolor("purple", alpha = 0.5), cex = 3)

