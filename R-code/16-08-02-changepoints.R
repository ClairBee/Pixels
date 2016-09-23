

####################################################################################################

# MARS ####

library(mda)

df <- df[!is.na(df$y),]

mm <- mars(df$x, df$y)

summary(mm)

plot(df, type = "o", pch = 20, cex = 0.4)
abline(v = mm$cuts, col = "red")
lines(mm$fitted, col = "red")

showcuts <- function(obj)
{
    tmp <- obj$cuts[obj$sel, ]
    dimnames(tmp) <- list(NULL, names(df))
    tmp
}
showcuts(mm)

plot(mm$fitted, type = "l")


####################################################################################################

# STRUCCHANGE       ####

system.time(bp <- breakpoints(y ~ x, h = 5, breaks = 4, data = df))

plot(df, type = "o", pch = 20, cex = 0.4)
abline(v = bp$breakpoints + 141, col = "red")

# finds dip exactly, but quite slow (~17s for a single half-column)


####################################################################################################

# SEGMENTED LINEAR REGRESSION                                                                   ####

library(segmented)

df <- data.frame(x = 1:1024, y = res[773,1:1024, "141009"])

# start by predicting linear model
linear <- lm(y ~ x, data = df[!is.na(df$y),])

plot(df, type = "o", pch = 20, cex = 0.4)
lines(predict(linear, newdata = data.frame(x = 1:1024)), col = "red", lwd = 2)

seg <- segmented(linear, seg.Z = ~ x, psi = NA, control = seg.control(display = F, K = 3))
plot.segmented(seg)

rp <- rpart(y ~ x, df)
plot(rp)

find.segments <- function(vv, k = 3) {
    df <- data.frame(x = 1:1024, y = vv)
    df <- df[!is.na(df$y),]
    
    # fit initial linear model
    linear <- lm(y ~ x, data = df)
    
    seg <- segmented(linear, seg.Z = ~ x, psi = NA, control = seg.control(display = F, K = k))
    if (exists(seg)) {seg$fitted.values} else {NULL}
}

system.time(zz <- apply(res[200:300, 1:1024, "141009"], 1, find.segments, k = 3))

plot(zz[,1])
####################################################################################################

# example found at https://rpubs.com/MarkusLoew/12164

# -------------------
# analyse breakpoints
# -------------------
# http://cran.r-project.org/doc/Rnews/Rnews_2008-1.pdf
library(segmented)

# have to provide estimates for breakpoints.
# after looking a the data, 
my.seg <- segmented(my.lm, 
                    seg.Z = ~ DistanceMeters, 
                    psi = list(DistanceMeters = c(4, 15)))

# When not providing estimates for the breakpoints "psi = NA" can be used.
# The number of breakpoints that will show up is not defined
#my.seg <- segmented(my.lm, 
#                    seg.Z = ~ DistanceMeters, 
#                    psi = NA)

# display the summary
summary(my.seg)

####################################################################################################
####################################################################################################

# EARLIER ATTEMPTS AT IDENTIFYING SCREEN SPOTS THROUGH MORPHOLOGICAL METHODS                    #### 

spoots <- function(bright.image, smooth.span = 1/5, min.diam = 5, midline = 1024.5, enlarge = F, auto.threshold = T, ignore.edges = 40) {
    
    # strip out padding to retain only active image region
    ar <- apply(which(!is.na(bright.image), arr.ind = T), 2, range)
    bright.im <- bright.image[ar[1,"row"]:ar[2,"row"], ar[1,"col"]:ar[2,"col"]]
    im.dims <- dim(bright.im)
    midline <- midline - min(which(!is.na(bright.image), arr.ind = T)[,"col"])
    
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    # apply lowess smoothing
    # don't use default smoothing in this case - don't want smoother to be drawn into dips, so use lower proportion
    smoo <- lowess.per.column(bright.im, midline = midline, span = smooth.span)
    
    res <- bright.im - smoo
    
    # flatten further by setting brighter pixels to mean value 
    res[res > mean(res) - sd(res)] <- mean(res) - sd(res)
    
    # dilate resulting image
    dilated <- dilate(res, sk)
    
    # erode resulting image (complete morphological closing)
    eroded <- erode(dilated, sk)
    
    # use k-means thresholding to identify spots
    # use 1-thresholded value to assign 1 to spots, 0 to background
    if (auto.threshold) {
        dim.sp <- array(1, dim = im.dims) - threshold(eroded, method = "kmeans")
    } else {
        dim.sp <- array(1, dim = im.dims) - threshold(eroded, mean(eroded) - 3*sd(eroded))
    }
    
    # if spot identification should be as large as possible, 
    if (enlarge) {dim.sp <- dilate(dim.sp, sk)}
    
    # re-pad image
    sp <- array(dim = dim(bright.image))
    sp[ar[1,"row"]:ar[2,"row"], ar[1,"col"]:ar[2,"col"]] <- dim.sp
    
    # remove any spots whose midpoint lies within the edge region
    if (ignore.edges > 0) {
        blobs <- clump(m2r(sp), dir = 4)
        
        sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                               id = getValues(blobs)[!is.na(getValues(blobs))]),
                    .(id), summarise, 
                    x.midpoint = round(mean(x),0), y.midpoint = round(mean(y),0),
                    size = length(x))
        
        sc$n.id <- sc$id
        sc$n.id[(sc$y.midpoint >= im.dims[2] - ignore.edges) | (sc$y.midpoint <= ignore.edges) | 
                    (sc$x.midpoint >= im.dims[1] - ignore.edges) | (sc$x.midpoint <= ignore.edges)] <- NA
        
        if (nrow(sc[!is.na(sc$n.id),]) > 0) {
            blobs <- subs(blobs, sc[,c("id", "n.id")])
            dim.sp <- bpx2im(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))), type = 1))
        } else {
            dim.sp <- array(0, dim = im.dims)
        }
    }
    
    if (sum(dim.sp == 1) == 0) {
        return(NULL)
    } else {
        return(which(dim.sp == 1, arr.ind = T))
    }
}

sp1 <- spoots(pw.m[,,"white", "141009"])
sp2 <- spoots(pw.m[,,"white", "131122"])

pixel.plot(spots, cex = 0.4, col = "cyan3")
points(sp1, cex = 0.4, pch = 15)

sp1.7 <- screen.spots(pw.m[,,"white", "141009"], min.diam = 7)
sp2.7 <- screen.spots(pw.m[,,"white", "131122"], min.diam = 7)

pixel.plot(spots, cex = 0.4, col = "cyan3")
points(sp1.7, cex = 0.4, pch = 15)

pixel.plot(spots2, cex = 0.4, col = "cyan3")
points(sp2.7, cex = 0.4, pch = 15)

sd1 <- clump.centres(spots, res[,,"141009"])

####################################################################################################


# CHECK DEPTH  ####

spots <- screen.spots(pw.m[,,"white", "141009"])
spots2 <- screen.spots(pw.m[,,"white", "131122"])

spot.d <- res[,,"141009"][spots]
spot2.d <- res[,,"131122"][spots2]

clump.centres <- function(px, vals) {
    
    # clump adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = c(2048, 2048))), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    xy$v <- vals[as.matrix(xy[,1:2])]
    
    df <- ddply(xy, .(id), summarise,
                xm = mean(x), ym = mean(y),
                r = ceiling(max(max(x) - min(x), max(y) - min(y)) / 2),
                q1.c = quantile(v, 0.25, na.rm = T),
                size = length(id))
    df
}

tt <- clump.centres(spots, res[,,"141009"])

tt2 <- clump.centres(spots2, res[,,"131122"])

sd(res[,,"131122"], na.rm = T)

sd(res[,,"141009"], na.rm = T)

tt$d.sd <- tt$q1.c / sd(res[,,"141009"], na.rm = T) 
tt
tt2$d.sd <- tt2$q1.c / sd(res[,,"131122"], na.rm = T) 
tt2

####################################################################################################

# MA SMOOTHING OF RESIDUALS, THEN CLOSING                                                       ####

res.smoo <- array(dim = c(2048, 2048))
n <- 11
res.smoo[,1:1024] <- t(apply(res[,1:1024,"141009"], 1, filter, rep(1/n, n)))
res.smoo[,1025:2048] <- t(apply(res[,1025:2048,"141009"], 1, filter, rep(1/n, n)))
res.smoo[res.smoo > 0] <- 0

xt <- threshold(res[,,"141009"], level = - 3 * sd(res[,,"141009"], na.rm = T), method = "literal")
xl.cl <- closing(xt, shapeKernel(c(5,5), type = "disc"))

cl <- closing(res[,,"141009"], shapeKernel(c(7,7), type = "disc"))
image(xl.cl)

o.plot(res[733,1:1024,"141009"])
lines(res.smoo[733,], col = "blue", lwd = 3)
lines(cl[733,], col = "red", lwd = 2)

image(cl)

op <- opening(res.smoo,  shapeKernel(c(5, 5), type = "disc"))
image(op)

plot(cl[773,1:1024], type = "o", pch = 20, cex = 0.4)

th <- res.smoo - cl

pixel.image(res.smoo - op)

####################################################################################################

# invert & use gaussian kernel ####

rr <- -res[,,"141009"]
pixel.image(rr)

gk <- gaussianKernel(sigma = c(1.5, 1.5))

g.cl <- opening(rr, gk)
image(g.cl)

o.plot(g.cl[733,1:1024])

####################################################################################################

# BLAH ####
ss <- screen.spots(pw.m[,,"white", "141009"])
spots <- apply(pw.m[,,"white", ], 3, screen.spots, enlarge = T, ignore.edges = 40)
spots.g <- apply(pw.m[,,"grey", ], 3, screen.spots, enlarge = T, ignore.edges = 40)

saveRDS(spots, paste0(fpath, "spots-w.rds"))
saveRDS(spots.g, paste0(fpath, "spots-g.rds"))

unlist(sapply(spots, nrow))
unlist(sapply(spots.g, nrow))

lapply(dimnames(pw.m)[[4]], 
       function(dt) {
           pixel.plot(spots[[dt]], cex = 0.4, col = "cyan3", main = paste0(dt, ": ", nrow(spots[[dt]]), ", ", nrow(spots.g[[dt]])))
           points(spots.g[[dt]], cex = 0.4, pch = 15)
       })

####################################################################################################

# LOWESS-SMOOTHING PER COLUMN                                                                   ####

smoo <- array(apply(pw.m[,,"white", ], 3, lowess.per.column, midline = 1024.5, span = 1/5), dim = c(2048, 2048, 21))
res <- pw.m[,,"white",] - smoo

rr <- res[,,"loan"]
rr[rr > median(rr, na.rm = T)] <- median(rr, na.rm = T)
pixel.image(rr)

sk <- shapeKernel(c(5, 5), type = "disc")
dil <- dilate(rr, sk)
dil[is.infinite(dil)] <- NA

pixel.image(dil)
er <- erode(dil, sk)
er[is.infinite(er)] <- NA

sd(er, na.rm = T)
pixel.plot(which(er < median(er, na.rm = T) - 3 * sd(er, na.rm = T), arr.ind  = T), cex = 0.4)

pixel.image(er, xlim = c(1000, 1200), ylim = c(500, 650))

hh <- threshold(er, method = "kmeans")

hist(pw.m[,,"white", "131122"], breaks = "fd")

df <- data.frame(sd = apply(res, 3, sd, na.rm = T),
                 mad = apply(res, 3, mad, na.rm = T),
                 sd.gain = unlist(lapply(dimnames(res)[[3]], function(dt) sd(res[,,dt][which(pw.m[,,"white", dt] > 20000)], na.rm = T))))

hist(res[,,1], breaks = "fd", ylim = c(0,30))
pixel.plot(which(res[,,1] < median(res[,,1], na.rm = T) - 6 * mad(res[,,1], na.rm = T), arr.ind = T), cex = 0.4)

pixel.plot(which(res[,,"141009"] < median(res[,,"141009"], na.rm = T) - 6 * mad(res[,,"141009"], na.rm = T), arr.ind = T), cex = 0.4)

lapply(dimnames(res)[[3]], 
       function(dt) {
           pixel.plot(which(res[,,dt] < median(res[,,dt], na.rm = T) - 4 * mad(res[,,dt], na.rm = T), arr.ind = T), 
                      cex = 0.4, main = dt, col = "cyan3")
           points(which(res[,,dt] < median(res[,,dt], na.rm = T) - 6 * mad(res[,,dt], na.rm = T), arr.ind = T), cex = 0.4, pch = 15)
       })


# median - 6 x MAD seems like a good start

pixel.image(pw.m[,,"grey", "loan"], xlim = c(1000, 1100), ylim = c(500,600))
points(ss, pch = 0, cex = 0.7)

ss2 <- screen.spots(pw.m[,,"white", "131122"], auto.threshold = F)

pixel.image(pw.m[,,"grey", "131122"])
points(ss2, pch = 0, cex = 0.3)

cl[is.infinite(cl)] <- NA
th <- threshold(cl, method = "kmeans")

# Gaussian kernel?

gk <- gaussianKernel(c(2,2), size = c(5,5))

rr <- res[,,"141009"]
pixel.image(rr)

dd <- closing(rr, sk)

image(dd)

sf <- sobelFilter(rr, axis = 0)

pixel.image(sf)

# can use sobel filter to identify background pattern!



####################################################################################################

####################################################################################################

#  MESSING ABOUT                                                                                ####


# smoothing span makes very little difference for large spots
plot(pw.m[773,,"white", "141009"], type = "l", xlim = c(1,1024))
lines(lowess(pw.m[773,1:1024,"white", "141009"], f = 1/5), col = "red", lwd = 2)
lines(lowess(pw.m[773,1:1024,"white", "141009"], f = 1/10), col = "purple", lwd = 2)
lines(lowess(pw.m[773,1:1024,"white", "141009"], f = 1/15), col = "blue", lwd = 2)

plot(pw.m[773,1:1024,"white", "141009"] - lowess(pw.m[773,1:1024,"white", "141009"], f = 1/5)$y, col = "red", type = "l")
lines(pw.m[773,1:1024,"white", "141009"] - lowess(pw.m[773,1:1024,"white", "141009"], f = 1/10)$y, col = "purple")
lines(pw.m[773,1:1024,"white", "141009"] - lowess(pw.m[773,1:1024,"white", "141009"], f = 1/15)$y, col = "blue")
lines(pw.m[773,1:1024,"white", "141009"] - lowess(pw.m[773,1:1024,"white", "141009"], f = 1/20)$y, col = "green3")

pixel.image(res[,,"131122"])
ss.131122 <- screen.spots(pw.m[,,"white", "131122"])

pixel.plot(ss.131122, cex = 0.4)

plot(pw.m[948,,"white", "131122"], type = "l", xlim = c(1,1024))
lines(lowess(pw.m[948,1:1024,"white", "131122"], f = 1/5), col = "red", lwd = 2)

plot(res[948,,"131122"], type = "l", xlim = c(1,1024))

# smooth the residuals?
pixel.image(res[,,"141009"])

res.smoo <- array(dim = c(2048, 2048))
n <- 11
res.smoo[,1:1024] <- t(apply(res[,1:1024,"141009"], 1, filter, rep(1/n, n)))
res.smoo[,1025:2048] <- t(apply(res[,1025:2048,"141009"], 1, filter, rep(1/n, n)))

pixel.image(res.smoo)

plot(res[773,,"141009"], type = "l", xlim = c(1,1024))
lines(res.smoo[773,], col = "red", lwd = 2)


####################################################################################################

# CHANGEPOINTS                                                                                  ####

library(bcp); library(cpm); library(ecp); library(segmented); library(changepoint)

cc <- bcp(res.smoo[773,])

# ecp
{
    ed <- e.divisive(md7[773,146:1019])
}

# segmented

df <- data.frame(x = c(1:1024), y = res.smoo[773,1:1024])
out.lm <- lm(y ~ x, data = df)
seg <- seg.lm.fit(res.smoo[773,])

# examples from manual
{
    set.seed(12)
    xx<-1:100
    zz<-runif(100)
    yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
    dati<-data.frame(x=xx,y=yy,z=zz)
    out.lm<-lm(y~x,data=dati)
    
    #simple example: 1 segmented variable, 1 breakpoint: you do not need to specify 
    # the starting value for psi
    o<-segmented(out.lm,seg.Z=~z)
    o<-segmented(out.lm,seg.Z=~x,psi=c(30,60),
                 control=seg.control(display=FALSE))
}
{
    library(survival)
    data(stanford2)
    
    o <- coxph(Surv(time, status) ~ age, data = stanford2)
    os <- segmented(o, ~ age, psi = 40) #estimate the breakpoint in the age effect
    summary(os) #actually it means summary.coxph(os)
    plot(os) #it does not work
    plot.segmented(os) #call explicitly plot.segmented() to plot the fitted piecewise lines
    
    plot(stanford2$time, stanford2$age)
}