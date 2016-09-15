
# screen spot deregistration between black & white images

# useful to note correlation between observed & fitted values
# eg. 130701: 99.9% within bounds, 99.4% outside, 99.88% overall

library("IO.Pixels"); library("CB.Misc")

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")

####################################################################################################

# UNIVERSAL SPOT-LOCATING FUNCTION                                                              ####

# add option to check white image, linear residuals, or both
# check horizontally as well as vertically
# check positive residuals as well as negative: does this work better or worse than registration?
# clean out edge effects
# filter spots by size?

dt <- "130701"
im <- fit.w.lm(pw.m[,,,dt])      # run on linear residuals

spots <- function(im, smooth.span = 1/15, min.diam = 5, edge.trim = 10) {
    
    # offset adjustment not necessary in linear residuals - difference generally much smaller

    # apply lowess smoothing in both directions, get residuals
    col.sm <- t(apply(im, 1, function(cc) lowess(cc, f = smooth.span)$y))
    row.sm <- apply(im, 2, function(rr) lowess(rr, f = smooth.span)$y)
    
    # further trim edges & combine row & column residuals
#    col.rng <- apply(which(!is.na(col.sm), arr.ind = T), 2, range) + edge.trim * c(1,-1,1,-1)
#    row.rng <- apply(which(!is.na(row.sm), arr.ind = T), 2, range) + edge.trim * c(1,-1,1,-1)
    
    # combine row & column residuals
    res <- apply(abind(im - col.sm, im - row.sm, along = 3), 1:2, mean, na.rm = T)
    
    # truncate residual values at median
    res.high <- res.low <- res
    res.high[res.high < 0] <- NA
    res.low[res.low > 0] <- NA
    
    # morphological closing
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    cl.bright <- opening(res.high, sk)
    cl.dim <- closing(res.low, sk)
    
    cl.bright[is.infinite(cl.bright)] <- cl.dim[is.infinite(cl.dim)] <- NA
    
    # trim edges 
    cl.bright[c(1:edge.trim, ncol(cl.bright) - 0:edge.trim),] <- NA
    
    # thresholding
    px.bright <- which(cl.bright > mad(res, na.rm = T), arr.ind = T)
    px.dim <- which(cl.dim < -mad(res, na.rm = T), arr.ind = T)
    
    # filter by size
    cc <- clump(m2r(cl.bright > mad(res, na.rm = T) & !is.infinite(cl.bright)))
    
    xy.bright <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                       id = getValues(cc)[!is.na(getValues(cc))])
    df.bright <- ddply(xy.bright, .(id), summarise,
                       xm = mean(x), ym = mean(y), size = length(x))

    hist(df.bright$size, breaks = "fd")
    pixel.plot(xy.bright[xy.bright$id %in% df.bright$id[df.bright$size > 20],])
    
    pixel.plot(px.dim)
}

# original function, for reference
screen.spots <- function(im, min.diam = 5, smooth.span = 1/15, midline = 1024.5, edge.crop = 10, coords = T) {
    
    # offset correction for upper vs lower panels (if midline exists)
    if (!is.na (midline)) {
        up <- apply(im[, floor(midline) + c(1:100)], 1, median, na.rm = T)
        lp <- apply(im[, floor(midline) + c(0:-100)], 1, median, na.rm = T)
        
        im[, ceiling(midline):dim(im)[[2]]] <- im[, ceiling(midline):dim(im)[[2]]] - (up - lp)
    }
    
    # Lowess smoothing over all columns
    # (faster than Loess & easier to apply in this form)
    res <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = smooth.span), "[[", 2))
    med.res <- median(res, na.rm = T)
    
    # truncate residual values at median
    tr <- res
    tr[res > med.res] <- med.res
    
    # morphological closing
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    cl <- closing(tr, sk)
    
    # pad with median residual value
    cl[is.infinite(cl) | is.na(cl)] <- med.res
    
    # threshold at 1 SD below median
    th <- threshold(cl, method = "literal", level = - sd(res, na.rm = T))
    
    # remove any pixels that lie at edge of active area
    th[, apply(which(!is.na(res), arr.ind = T), 2, range)[,2]] <- 1 # residual artefact at edge of smoothed area
    th[c(1 + c(0:edge.crop), dim(th)[1] - c(0:edge.crop)), ] <- 1    # panel edges
    th[, c(1 + c(0:edge.crop), dim(th)[1] - c(0:edge.crop))] <- 1    # panel edges
    
    # enlarge by eroding with same structuring element
    exp <- erode(th, sk)
    
    if (coords) {return(which(exp == 0, arr.ind = T))} else {return(exp)}
}

####################################################################################################

# TRY TO FIND FAINT SPOTS IN GREY/WHITE IMAGES                                                  ####

dt <- "130701"

# plot using gradient & direct thresholding in grey, white & linear residuals of each acquisition
fpath <- "./Image-plots/spots-per-image/"

invisible(lapply(names(pw.m),
                 function(dt) {
                     bmp(paste0(fpath, dt, "-spots.bmp"), height = 480 * 2, width = 480 * 3)
                     par(mfrow = c(2,3), mar = c(2,2,3,1))
                     
                     g <- pw.m[,,"grey", dt]
                     w <- pw.m[,,"white", dt]
                     linear <- fit.w.lm(pw.m[,,,dt], res.only = F)
                 }))

####################################################################################################

# INDIVIDUAL RAW IMAGES FROM LOAN PANEL                                                          ####

# is absence of shadow in white image due to blurring?
loan <- array(dim = c(2048, 2048, 20))
loan[25:2024, 25:2024,] <- abind(lapply(list.files("./Image-data/loan/white/", pattern = "\\.tif$", full.names = T), 
                                        function(nm) {
                                            tmp <- readTIFF(nm, as.is = T)
                                            t(tmp[nrow(tmp):1, , drop = FALSE])
                                        }), along = 3)

loan.spots <- apply(loan, 3, screen.spots)

invisible(lapply(loan.spots, pixel.plot, cex = 0.2))

####################################################################################################

# SPOTS BY MORPHOLOGICAL GRADIENT                                                               ####

spots.by.gradient <- function(im, smoothing.span = 1/5, min.diam = 5, ignore.edge = 5) {
    
    smoo <- do.call("rbind", lapply(apply(im, 1, lowess, f = smoothing.span), "[[", 2))
    res <- im - smoo
    
    # distinguish between high & low residuals
    res.high <- res.low <- res
    res.high[res.high > mad(res, na.rm = T)] <- 0
    res.low[res.low < -mad(res, na.rm = T)] <- 0
    
    # find enhanced morphologica gradient of residuals
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    res.gradient <- dilate(res.high, sk) - erode(res.low, sk)
    
    # trim edges by specified width
    rng <- apply(which(!is.na(im), arr.ind = T), 2, range) + ignore.edge * c(1,-1,1,-1)
    res.gradient[c(1:rng[1,1], rng[2,1]:ncol(im)),] <- NA
    res.gradient[,c(1:rng[1,2], rng[2,2]:nrow(im))] <- NA
    
    return(which(res.gradient <= 0 & !is.na(res.gradient) & !is.infinite(res.gradient), arr.ind = T))
}

sp.white <- spots.by.gradient(pw.m[,,"white", "141009"])
sp.org <- screen.spots(pw.m[,,"white", "141009"])
# need to filter by spot size as well

# linear residuals; gradient; spots from linear residuals by gradient; values of bright & dim identified
# white residuals; gradient; spots by standard method; spots by gradient

####################################################################################################


# PLOTS OF SCREEN SPOTS + REGISTRATION IN ALL IMAGES                                            ####

# temporary support functions

smoothed.res <- function(im, midline = 1024.5) {
    
    # adjust offset across midline first
    if (!is.na(midline)) {
        up <- apply(im[, floor(midline) + c(1:100)], 1, median, 
                    na.rm = T)
        lp <- apply(im[, floor(midline) + c(0:-100)], 1, median, 
                    na.rm = T)
        im[, ceiling(midline):dim(im)[[2]]] <- im[, ceiling(midline):dim(im)[[2]]] - 
            (up - lp)
    }
    
    smoo <- do.call("rbind", lapply(apply(im, 1, lowess, f = 1/5), "[[", 2))
    im - smoo
}

morph.gradient <- function(im, min.diam = 5) {
    # distinguish between high & low residuals
    res.high <- res.low <- im
    res.high[res.high > mad(im, na.rm = T)] <- 0
    res.low[res.low < -mad(im, na.rm = T)] <- 0
    
    # find enhanced morphologica gradient of residuals
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    dilate(res.high, sk) - erode(res.low, sk)
}


fpath <- "/home/clair/Documents/Pixels/Image-plots/Linear-screen-spots/"

invisible(lapply(dimnames(pw.m)[[4]][8:21],
                 function(dt) {
                     img <- pw.m[,,,dt]
                     linear <- fit.w.lm(img)

                     sm.white <- smoothed.res(img[,,"white"]) 
                     sm.linear <- smoothed.res(linear)
                     
                     grad.linear <- morph.gradient(sm.linear)
                     grad.white <- morph.gradient(sm.white)
                     
                     px.linear <- which(grad.linear <= 0 & !is.infinite(grad.linear), arr.ind = T)
                     px.white <- which(grad.white <= 0 & !is.infinite(grad.white), arr.ind = T)
                     
                     spots.linear <- spots.white <- array(dim = dim(img[,,1]))
                     spots.linear[px.linear] <- linear[px.linear]
                     spots.white[px.white] <- img[,,"white"][px.white]
                     
                     sp.white <- screen.spots(img[,,"white"])
                     sp.linear <- screen.spots(linear)
                     
                     bmp(paste0(fpath, "linear-screen-spots-", dt, ".bmp"), height = 960, 
                         width = 480 * 3); {
                         par(mfrow = c(2, 3), mar = c(2,2,3,1))
                         
                         # plots of linear residual processing
                         pixel.image(linear, title = paste0(dt, " - linear residuals"))
                         
                         pixel.image(linear, title = "spots by m.gradient")
                         if (nrow(px.linear) > 0) draw.outlines(px.linear)
                         
                         pixel.image(linear, title = "spots by thresholding")
                         if (nrow(sp.linear) > 0) draw.outlines(sp.linear)
                         
                         # pixel.image(sm.linear, title = "Loess residuals")
                         # pixel.image(grad.linear, title = "m. gradient")
                         # pixel.image(spots.linear, break.levels = sd.levels(linear))
                         # draw.outlines(sp.linear)
                         # draw.outlines(px.white, col = "skyblue")
                         
                         # plots of white value processing
                         pixel.image(img[,,"white"], title = paste0(dt, " - white values"))
                         
                         pixel.image(img[,,"white"], title = "spots by m.gradient")
                         if (nrow(px.white) > 0) draw.outlines(px.white)
                         
                         pixel.image(img[,,"white"], title = "spots by thresholding")
                         if (nrow(sp.white) > 0) draw.outlines(sp.white)
                         #pixel.image(sm.white, title = "Loess residuals - white values")
                         #pixel.image(grad.white, title = "m. gradient")
                         #pixel.image(spots.white, break.levels = sd.levels(img[,,"white"]))
                         #draw.outlines(sp.white)
                         #draw.outlines(px.linear, col = "skyblue")
                         
                         dev.off()
                     }
                 }))


####################################################################################################

# SCREEN SPOTS IN GREY & WHITE IMAGES                                                           ####

dt <- "141118"

ss.grey <- screen.spots(pw.m[,,"grey", dt])
ss.white <- screen.spots(pw.m[,,"white", dt])

pixel.plot(ss.grey)
draw.outlines(ss.white, col = "red")

# struggling to identify screen spots at all in eg. 130705

####################################################################################################

# SCREEN SPOTS IN LINEAR RESPONSE                                                               ####

dt <- "130701"

linear <- fit.lm(pw.m[,,, dt], "g ~ b * w", res.only = T)

pixel.image(linear)
hist(linear, breaks = "fd", xlim = c(-1000,1000), ylim = c(0,100))

# should also adjust column offset for midline in final function

# get loess-smoothed residuals
smoo <- do.call("rbind", lapply(apply(linear, 1, lowess, f = 1/5), "[[", 2))

# try 1/5 as smoothing span. Edges are less important here, since we're looking for registration only
o.plot(linear[510,], xlim = c(600,900))
lines(smoo[510,], col = "blue")

res <- linear - smoo

pixel.image(res)
s.hist(res, xlim = c(-200,200), ylim = c(0,200))

res.high <- res.low <- res
res.high[res.high > mad(res, na.rm = T)] <- 0
res.low[res.low < -mad(res, na.rm = T)] <- 0

pixel.plot(which(res > 0, arr.ind = T), cex = 0.1, main = "high residuals")
points(which(res > mad(res, na.rm = T), arr.ind = T), col = "red", cex = 0.1, pch = 15)

pixel.plot(which(res < 0, arr.ind = T), cex = 0.1, main = "low residuals")
points(which(res < -mad(res, na.rm = T), arr.ind = T), col = "blue", cex = 0.1, pch = 15)

# morphological opening of higher residuals
sk <- shapeKernel(c(5, 5), type = "disc")
cl.high <- opening(res.high, sk)

cc <- clump(m2r(cl.high))

xy <- data.frame(xyFromCell(cc, which(getValues(cc) > 0 & !is.na(getValues(cc)))),
                 id = getValues(cc)[which(getValues(cc) > 0 & !is.na(getValues(cc)))])
df <- ddply(xy, .(id), summarise,
            xm = mean(x), ym = mean(y), size = length(x))

sp.high <- as.matrix(xy[xy$id %in% df$id,1:2])


# morphological closing of lower residuals
cl.low <- closing(res.low, sk)

cc <- clump(m2r(cl.low))

xy <- data.frame(xyFromCell(cc, which(getValues(cc) > 0 & !is.na(getValues(cc)))),
                 id = getValues(cc)[which(getValues(cc) > 0 & !is.na(getValues(cc)))])
df <- ddply(xy, .(id), summarise,
            xm = mean(x), ym = mean(y), size = length(x))

sp.low <- as.matrix(xy[xy$id %in% df$id,1:2])

# MORPHOLOGICAL GRADIENT                                                                        ####  

# should highlight edges between light and dark

dt <- "130701"

linear <- fit.lm(pw.m[,,, dt], "g ~ b * w", res.only = T)
smoo <- do.call("rbind", lapply(apply(linear, 1, lowess, f = 1/5), "[[", 2))
res <- linear - smoo


# distinguish between high & low residuals
res.high <- res.low <- res
res.high[res.high > mad(res, na.rm = T)] <- 0
res.low[res.low < -mad(res, na.rm = T)] <- 0

# find morphological gradient
sk <- shapeKernel(c(5, 5), type = "disc")
res.gradient <- dilate(res.high, sk) - erode(res.low, sk)
pixel.image(res.gradient)

hist(res.gradient, breaks = "fd")

# could use either 0 cutoff or asymmetric bounds. Check on further images
pixel.plot(which(res.gradient <= 0, arr.ind = T), cex = 0.1)
pixel.plot(which(res.gradient <= asymm.bounds(res.gradient)[1], arr.ind = T), cex = 0.1, col = "red")

sp <- which(res.gradient <= 0 & !is.na(res.gradient) & !is.infinite(res.gradient), arr.ind = T)

# use further dilation/erosion to expand area picked up

# extract original residual values at those points
spot.values <- array(dim = dim(res))
spot.values[sp] <- res[sp]

# cluster bright & dark patches

# match bright & dark patches to nearest neighbour

# sort nearest-neighbour distances in ascending order, stop when pairing becomes impossible

# check registration distance

# is even the fact of identifying spots like this enough evidence?

# plots
{
    bmp(paste0(fpath, "linear-residuals.bmp"), height = 2048, width = 2048)
    pixel.image(linear)
    dev.off()
}

{
    bmp(paste0(fpath, "morphological-gradient.bmp"), height = 2048, width = 2048)
    pixel.image(res.gradient)
    dev.off()
}

{
    bmp(paste0(fpath, "spots-identified.bmp"), height = 2048, width = 2048)
    pixel.image(spot.values)
    dev.off()
}

{
    sp.141009 <- screen.spots(pw.m[,,"white", "141009"])
    
    bmp(paste0(fpath, "morphological-gradient-images-aligned.bmp"), height = 2048, width = 2048)
    pixel.image(res.gradient)
    draw.outlines(sp.141009)
    dev.off()
}

#======================================================================================
# previous attempts

# original version of column smoothing algorithm
col.smoo <- apply(im, 1,
                  function(ll) {
                      # adjust for offset at midline (if it exists)
                      if (!is.na(midline)) {
                          up <- median(ll[ceiling(midline) + c(0:100)], na.rm = T)
                          lp <- median(ll[floor(midline) - c(0:100)], na.rm = T)
                          ll[ceiling(midline):length(ll)] <- ll[ceiling(midline):length(ll)] - (up - lp)
                      }
                      lowess(ll, f = smooth.span)$y
                  })


{
res.gradient1 <- dilate(res, sk) - erode(res, sk)

pixel.image(res.gradient1)

# dilation will enlarge bright regions, shrinking dark ones
res.gradient2 <- dilate(res.high, sk) - erode(res.low, sk)
pixel.image(res.gradient2)

pixel.image(linear)
pixel.image(res)

# gradient of linear image directly
res.gradient.org <- dilate(linear, sk) - erode(linear, sk)
pixel.image(res.gradient.org)
}