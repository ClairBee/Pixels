
library("IO.Pixels"); library("CB.Misc")

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
fpath <- "./Image-plots/spots-per-image/"

####################################################################################################

# offset correction for row/panel edges                                             ####

adjust.midline.offset <- function(cc, midline = 1024.5, match.width = 100) {
    os.u <- median(cc[ceiling(midline) + c(0:match.width)], na.rm = T)
    os.l <- median(cc[floor(midline) - c(0:match.width)], na.rm = T)
    cc[ceiling(midline):length(cc)] <- cc[ceiling(midline):length(cc)] - (os.u - os.l)
    return(cc)
}

adjust.panel.offset <- function(rr, p.width = 128, match.width = 50, method = "median") {
    sp <- matrix(rr, nrow = p.width)    # 'fold' row vector into matrix of panel width
    
    # get differences across panel boundaries
    if (method == "median") {
        os <- apply(sp[p.width - (match.width:1) + 1,], 2, median, na.rm = T)[-ncol(sp)] - 
            apply(sp[1:match.width,], 2, median, na.rm = T)[-1]
    } else {
        os <- apply(sp[p.width - (match.width:1) + 1, -ncol(sp)] - sp[1:match.width, -1],
                    2, function(diff) diff[order(abs(diff))][floor(match.width * .25)])
    }
    rr + rep(c(0, cumsum(os)), each = p.width)
}

w.col.adj <- t(apply(im, 1, adjust.midline.offset))
w.row.adj <- apply(im, 2, adjust.panel.offset, method = "banana")

####################################################################################################

# lowess smoothing along rows & columns                                             ####

col.sm <- t(apply(im, 1, function(cc) lowess(cc, f = 1/15)$y))
row.sm <- apply(im, 2, function(rr) lowess(rr, f = 1/15)$y)

# average residuals in both directions: may remove need to adjust offset?
res <- apply(abind(im - col.sm, im - row.sm, along = 3), 1:2, mean, na.rm = T)

col.adj.sm <- t(apply(w.col.adj, 1, function(cc) lowess(cc, f = 1/15)$y))
row.adj.sm <- apply(w.row.adj, 2, function(rr) lowess(rr, f = 1/15)$y)

adj.res <- apply(abind(w.col.adj - col.adj.sm, w.row.adj - row.adj.sm, along = 3), 1:2, mean, na.rm = T)

####################################################################################################

# filter residuals into low vs high residuals                                       ####

res.high <- res.low <- res
res.high[res.high < 0 | is.na(res.high)] <- 0
res.low[res.low > 0 | is.na(res.low)] <- 0

adj.res.high <- adj.res.low <- adj.res
adj.res.high[adj.res.high < 0 | is.na(adj.res.high)] <- 0
adj.res.low[adj.res.low > 0 | is.na(adj.res.low)] <- 0

####################################################################################################

# morphological processing                                                          ####

min.diam <- 5
sk <- shapeKernel(c(min.diam, min.diam), type = "disc")

# morphological gradient
res.gradient <- dilate(res.high, sk) - erode(res.low, sk)

# opening of high residuals, closing of low residuals
cl.bright <- opening(res.high, sk)
cl.dim <- closing(res.low, sk)

# for offset-adjusted images
adj.res.gradient <- dilate(adj.res.high, sk) - erode(adj.res.low, sk)

adj.cl.bright <- opening(adj.res.high, sk)
adj.cl.dim <- closing(adj.res.low, sk)

####################################################################################################

# threshold to identify likely candidates                                           ####

hist(cl.bright[cl.bright > 0], breaks = "fd")
abline(v = mad(res, na.rm = T), col = "red")

hist(cl.dim[cl.dim < 0], breaks = "fd")
abline(v = -mad(res, na.rm = T), col = "red")

px.bright <- which(cl.bright > mad(res, na.rm = T), arr.ind = T)
px.dim <- which(cl.dim < -mad(res, na.rm = T), arr.ind = T)

pixel.plot(px.bright, col = "red")
pixel.plot(px.dim, col = "blue")

adj.px.bright <- which(adj.cl.bright > mad(adj.res, na.rm = T), arr.ind = T)
adj.px.dim <- which(adj.cl.dim < -mad(adj.res, na.rm = T), arr.ind = T)

####################################################################################################

# enlarge identified spots by dilation                                              ####

px.d <- erode(adj.px.dim, sk)

####################################################################################################

# general plots etc                                     ####
sp.org <- screen.spots(pw.m[,,"white", dt])
draw.outlines(sp.org, col = "blue")

####################################################################################################

# FUNCTIONS                                                                                     ####

# support functions

adjust.midline.offset <- function(cc, midline = 1024.5, match.width = 100) {
    os.u <- median(cc[ceiling(midline) + c(0:match.width)], na.rm = T)
    os.l <- median(cc[floor(midline) - c(0:match.width)], na.rm = T)
    cc[ceiling(midline):length(cc)] <- cc[ceiling(midline):length(cc)] - (os.u - os.l)
    return(cc)
}

adjust.panel.offset <- function(rr, p.width = 128, match.width = 50, method = "median") {
    sp <- matrix(rr, nrow = p.width)    # 'fold' row vector into matrix of panel width
    
    # get differences across panel boundaries
    if (method == "median") {
        os <- apply(sp[p.width - (match.width:1) + 1,], 2, median, na.rm = T)[-ncol(sp)] - 
            apply(sp[1:match.width,], 2, median, na.rm = T)[-1]
    } else {
        os <- apply(sp[p.width - (match.width:1) + 1, -ncol(sp)] - sp[1:match.width, -1],
                    2, function(diff) diff[order(abs(diff))][floor(match.width * .25)])
    }
    rr + rep(c(0, cumsum(os)), each = p.width)
}

adjust.panels <- function(im, by.row = T) {
    
    col.adj <- t(apply(im, 1, adjust.midline.offset))
    
    if (by.row) {
        row.adj <- apply(im, 2, adjust.panel.offset, method = "not.median")
        return(list(col.adj = col.adj, row.adj = row.adj))
    } else {
        return(col.adj)
    }
}

lowess.smoothing <- function(col.im, row.im, smooth.span = 1/15) {
    
    col.sm <- t(apply(col.im, 1, function(cc) lowess(cc, f = smooth.span)$y))

    if (!missing(row.im)) {
        row.sm <- apply(row.im, 2, function(rr) lowess(rr, f = smooth.span)$y)
        
        # average residuals in both directions
        res <- apply(abind(col.im - col.sm, row.im - row.sm, along = 3), 1:2, mean, na.rm = T)
    } else {
        res <- col.im - col.sm
    }
    return(res)
}

morph.gradient <- function(res, min.diam = 5) {
    
    # filter residuals
    res.high <- res.low <- res
    res.high[res.high < 0 | is.na(res.high)] <- 0
    res.low[res.low > 0 | is.na(res.low)] <- 0
    
    # define structuring element
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    # return morphological gradient
    dilate(res.high, sk) - erode(res.low, sk)
}

morph.bright.dim <- function(res, min.diam = 5, bright = F) {
    
    # filter residuals
    res.high <- res.low <- res
    res.high[res.high < 0 | is.na(res.high)] <- 0
    res.low[res.low > 0 | is.na(res.low)] <- 0
    
    # define structuring element
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    dim <- closing(res.low, sk)
    
    if (bright) {
        bright <- opening(res.high, sk)
        return(list(bright = bright, dim = dim))
    } else {
        return(dim)
    }
}

get.px <- function(im, incl.rows = T, gradient = F, incl.bright = F) {
    
    # adjust panel-edge offsets
    os.adj <- adjust.panels(im, by.row = incl.rows)
    
    # lowess smoothing (per column and, if given, per row)
    if (incl.rows) {
        sm <- lowess.smoothing(col.im = os.adj$col.adj, row.im = os.adj$row.adj)
    } else {
        sm <- lowess.smoothing(os.adj)
    }
    
    # morphological processing
    if (gradient) {
        grad <- morph.gradient(sm)
        th <- grad[grad <= 0]
    } else {
        morphed <- morph.bright.dim(sm, bright = incl.bright)

        if (incl.bright) {
            th.dim <- which(morphed$dim < -mad(res, na.rm = T), arr.ind = T)
            th.bright <- which(morphe$bright > mad(res, na.rm = T), arr.ind = T)
        } else {
            th.dim <- which(morphed < -mad(sm, na.rm = T), arr.ind = T)
        }
    }
    
    # enlarge pixels identified
    if (incl.bright) {
        
    }
    
}

####################################################################################################

# PLOT SPOTS FOUND BY EACH METHOD FOR COMPARISON                                                ####

# white image: gradient by column only, gradient by rows & columns, 
#               low values by column only, low by row & column
# linear residuals: gradient by column only, gradient by rows & columns,
#               high & low values by column only, high & low by row & column

dt <- "loan"
min.diam <- 5
sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
px.size <- 0.4

# need to decide between direct thresholding & gradient approach in linear residuals.
invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     im <- pw.m[,,"white", dt]
                     linear <- fit.w.lm(pw.m[,,,dt])
                     
                     # find spots in linear residuals
                     {
                         # offset correction & lowess residuals along columns
                         l.adj <- t(apply(linear, 1, adjust.midline.offset))
                         l.res <- l.adj - t(apply(l.adj, 1, function(cc) lowess(cc, f = 1/15)$y))
                         
                         # filter residuals into high & low
                         l.res.high <- l.res.low <- l.res
                         l.res.high[l.res.high < 0 | is.na(l.res.high)] <- 0
                         l.res.low[l.res.low > 0 | is.na(l.res.low)] <- 0
                         
                         # morphological bits
                         l.grad <- dilate(l.res.high, sk) - erode(l.res.low, sk)
                         l.bright <- opening(l.res.high, sk)
                         l.dim <- closing(l.res.low, sk)
                         
                         # thresholding
                         l.grad.bright.th <- (l.grad > asymmetric.mad(l.grad[l.grad != 0])[2]) * 1
                         l.grad.dim.th <- (l.grad < asymmetric.mad(l.grad[l.grad != 0])[1]) * 1
                         l.dim.th <- (l.dim < -mad(l.res, na.rm = T)) * 1
                         l.bright.th <- (l.bright > mad(l.res, na.rm = T)) * 1

                         # expansion
                         l.grad.bright.th <- dilate(l.grad.bright.th, sk)
                         l.grad.dim.th <- dilate(l.grad.dim.th, sk)
                         l.dim.th <- dilate(l.dim.th, sk)
                         l.bright.th <- dilate(l.bright.th, sk)
                         
                         # extraction
                         l.grad.bright.px <- which(l.grad.bright.th > 0, arr.ind = T)
                         l.grad.dim.px <- which(l.grad.dim.th > 0, arr.ind = T)
                         l.dim.px <- which(l.dim.th > 0, arr.ind = T)
                         l.bright.px <- which(l.bright.th > 0, arr.ind = T)
                      }
                    
                     # find spots in white images
                     {
                         # offset correction & lowess residuals along columns
                         w.adj <- t(apply(im, 1, adjust.midline.offset))
                         w.res <- w.adj - t(apply(w.adj, 1, function(cc) lowess(cc, f = 1/15)$y))
                         
                         # filter residuals - retain only low
                         w.res.low <- w.res
                         w.res.low[w.res.low > 0 | is.na(w.res.low)] <- 0
                         
                         # morphological bits
                         w.dim <- closing(w.res.low, sk)
                         
                         w.dim.th <- (w.dim < -mad(w.res, na.rm = T)) * 1       # thresholding
                         w.dim.th <- dilate(w.dim.th, sk)                       # expansion
                         w.dim.px <- which(w.dim.th > 0, arr.ind = T)           # extraction
                    }
                     
                     # plots
                     bmp(paste0(fpath, "spot-identification-", dt, ".bmp"),
                         height = 960 * 2, width = 960 * 4); {
                             par(mfrow = c(2, 4), mar = c(2,2,3,1))
                             
                             pixel.image(linear, title = paste0(dt, " - linear residuals"))
                             pixel.image(l.res, title = "lowess residuals")
                             pixel.plot(l.bright.px, cex = px.size, col = "red", main = "direct thresholding")
                                points(l.dim.px, cex = px.size, col = "blue", pch = 15)
                             pixel.plot(l.grad.bright.px, cex = px.size, col = "red", main = "points by gradient")
                                points(l.grad.dim.px, cex = px.size, col = "blue", pch = 15)

                                
                             pixel.image(im, title = paste0(dt, " - white values"))
                             pixel.image(w.res, title = "lowess residuals")
                             pixel.plot(w.dim.px, cex = px.size, col = "blue", main = "direct thresholding")
                             
                             dev.off()
                         }
                     

                     
                 }))

####################################################################################################

# FINAL APPROACH                                                                                ####

# use vertical offset correction in linear & white: shouldn't be an offset, but has been observed.
sspots <- function(im, smooth.span = 1/15, min.diam = 5, midline = 1024.5, match.width = 100) {

    # support function to adjust offset across midline
    adjust.midline.offset <- function(cc, midline = 1024.5, match.width = 100) {
        os.u <- median(cc[ceiling(midline) + c(0:match.width)], na.rm = T)
        os.l <- median(cc[floor(midline) - c(0:match.width)], na.rm = T)
        cc[ceiling(midline):length(cc)] <- cc[ceiling(midline):length(cc)] - (os.u - os.l)
        return(cc)
    }
    
    # get white image & residuals from linear regression
    w <- im[,,"white"]
    linear <- fit.w.lm(im)
    
    # define kernel for morphological operations
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    #----------------------------------------------------------
    
    if(!is.na(midline)) {
        l.adj <- t(apply(linear, 1, adjust.midline.offset, midline = midline, match.width = match.width))
        w.adj <- t(apply(w, 1, adjust.midline.offset, midline = midline, match.width = match.width))
    } else {
        l.adj <- linear
        w.adj <- w
    }

    # find spots in linear residuals
    {
        # lowess residuals along columns
        l.res <- l.adj - t(apply(l.adj, 1, function(cc) lowess(cc, f = smooth.span)$y))
        l.res.high <- l.res.low <- l.res
        
        # dim spots in linear residuals
        {
            l.res.low[l.res.low > 0 | is.na(l.res.low)] <- 0        # filter low residuals
            l.dim <- closing(l.res.low, sk)                         # morphological closing
            l.dim.th <- (l.dim < -mad(l.res, na.rm = T)) * 1        # threshold
            l.dim.th <- dilate(l.dim.th, sk)                        # slightly enlarge all spots
            l.dim.px <- which(l.dim.th > 0, arr.ind = T)            # get pixel coordinates
        }
        
        # bright regions in linear residuals
        {
            l.res.high[l.res.high < 0 | is.na(l.res.high)] <- 0
            l.bright <- opening(l.res.high, sk)
            l.bright.th <- (l.bright > mad(l.res, na.rm = T)) * 1
            l.bright.th <- dilate(l.bright.th, sk)
            l.bright.px <- which(l.bright.th > 0, arr.ind = T)
        }
    }
    
    # find spots in white image
    {
        # offset correction & lowess residuals along columns
        w.res <- w.adj - t(apply(w.adj, 1, function(cc) lowess(cc, f = smooth.span)$y))
        
        # filter residuals - retain only dim regions
        w.res.low <- w.res
        w.res.low[w.res.low > 0 | is.na(w.res.low)] <- 0
        
        w.dim <- closing(w.res.low, sk)                         # morphological closing
        w.dim.th <- (w.dim < -mad(w.res, na.rm = T)) * 1        # threshold
        w.dim.th <- dilate(w.dim.th, sk)                        # slightly enlarge all spots
        w.dim.px <- which(w.dim.th > 0, arr.ind = T)            # get pixel coordinates
    }

    list(nl.dim = l.dim.px, nl.bright = l.bright.pxw.dim = w.dim.px)
}

invisible(lapply(dimnames(pw.m)[[4]][5:21],
                 function(dt) {
                     if (dt == "MCT225") {
                         midline <- NA
                     } else {
                         midline <- 1024.5
                     }
                     sp <- sspots(pw.m[,,,dt], midline = midline)
                     
                     bmp(paste0(fpath, "spots-found-", dt, ".bmp"),
                         height = 960 * 1, width = 960 * 3)
                     par(mfrow = c(1,3))
                     
                     pixel.image(pw.m[,,"white", dt], title = "white image")
                         if (nrow(sp$dim) > 0) {draw.outlines(sp$dim)}
                         if (nrow(sp$bright) > 0) {draw.outlines(sp$bright, col = "red")}

                     pixel.image(fit.w.lm(pw.m[,,,dt]), title = "linear residuals")
                         if (nrow(sp$dim) > 0) {draw.outlines(sp$dim)}
                         if (nrow(sp$bright) > 0) {draw.outlines(sp$bright, col = "red")}

                     if (nrow(sp$dim) > 0) {
                         pixel.plot(sp$dim, main = paste0(dt, " - spots identified"))
                         if(nrow(sp$bright) > 0) {
                             points(sp$bright, col = "red", pch = 15, cex = 0.4)
                             }
                     } else {
                         if(nrow(sp$bright) > 0) {
                             pixel.plot(sp$bright, col = "red")
                         }
                     }
                     
                     dev.off()
                 }))


