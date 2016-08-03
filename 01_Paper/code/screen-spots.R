
# normalise images for easier comparison?

library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()

# can't run over shading-corrected image: most of the affected area has been removed
# (so, no good for assessing degree of actual damage)

####################################################################################################

# PROCEDURE

# reduce upper-vs-lower offset by subtracting offset per column
# detrend each column using Lowess smoothing (small span preferable b/c lose less at edges)
# morphological approach: closing with circular kernel
# identify depth of screen spot & filter out shallow spots
# enlarge if necessary

# check that spots not corrected by sc are covered by identified pixels?

####################################################################################################

# COMPARE SMOOTHING SPANS                                                                       ####

pixel.image(pw.m[,,"white", "141009"], xlim = c(0,1024), ylim = c(0,1024))

sp.05 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/05)
sp.15 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/15)
sp.25 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/25)
sp.35 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/35)
sp.45 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/45)


draw.outlines(sp.05, lty = 3)       # misses points at edge
draw.outlines(sp.15, lty = 2)
draw.outlines(sp.25)
draw.outlines(sp.35, lwd = 2)       # identifies everything that the others do
draw.outlines(sp.45, lwd = 2, col = "red")       # identifies everything that the others do


# need to run over all images for full comparison. Takes ~13s per image.

####################################################################################################

# RUN OVER ALL IMAGES                                                                           ####

spots.15 <- apply(pw.m[,,"white", ], 3, screen.spots)
spots.15$"MCT225" <- screen.spots(pw.m[,,"white", "MCT225"], midline = NA)

spots.35 <- apply(pw.m[,,"white", ], 3, screen.spots, smooth.span = 1/35)

fpath <- "./Image-plots/screen-spots/"
saveRDS(spots.15, "./Other-data/screen-spots-1-15.rds")
saveRDS(spots.35, "./Other-data/screen-spots-1-35.rds")

spots.15.g <- apply(pw.m[,,"grey", ], 3, screen.spots)

# summary
{
    sp.summ <- data.frame(acq = names(spots.15),
                          n.px = sapply(spots.15, nrow),
                          prop = round(100 * sapply(spots.15, nrow) / 2048^2, 3),
                          stringsAsFactors = F)
    sp.summ <- merge(sp.summ, data.frame(acq = names(mids),
                              n.spots = sapply(mids, nrow)), by = "acq", all.x = T)[,c("acq", "n.spots", "n.px", "prop")]
}

# plot all spots found
{
    invisible(lapply(dimnames(pw.m)[[4]],
                     function(dt) {
                         bmp(paste0(fpath, "ss-", dt, ".bmp"))
                         par(mar = c(2,2,3,1))
                         pixel.image(pw.m[,,"white", dt], title = paste0(dt, "; ", nrow(spots[[dt]]), "px"))
                         if (nrow(spots[[dt]]) > 0) {draw.outlines(spots[[dt]])}
                         dev.off()
                     }))
}


max(100 * sapply(spots.15, nrow) / 2048^2)
central <- sapply(spots.15, function(sp) sp[sp[,1] %in% c(500:1500) & sp[,2] %in% c(500:1500),])
max(100 * sapply(central, nrow) / 2048^2)

sapply(spots.35, nrow)

# comparison of smoothing spans
{
    invisible(lapply(dimnames(pw.m)[[4]],
                     function(dt) {
                         if (nrow(spots.15[[dt]]) > 0) {
                             pixel.plot(spots.15[[dt]], cex = 0.4, col = "cyan3")
                         } else {
                             plot(0, type = "n", xlim = c(0, 2048), ylim = c(0, 2048), xlab = "", ylab = "")
                         }
                         title(dt)
                         if (nrow(spots.35[[dt]]) > 0) points(spots.35[[dt]], pch = 15, cex = 0.4)
                         abline(h = range(which(!is.na(pw.m[1024,,"white", dt]), arr.ind = T)) + c(69, -68), col = "red")
                     }))
}

# comparison of grey vs white
{
    invisible(lapply(dimnames(pw.m)[[4]],
                     function(dt) {
                         if (nrow(spots.15[[dt]]) > 0) {
                             pixel.plot(spots.15[[dt]], cex = 0.4, col = "cyan3")
                         } else {
                             plot(0, type = "n", xlim = c(0, 2048), ylim = c(0, 2048), xlab = "", ylab = "")
                         }
                         title(dt)
                         if (nrow(spots.15.g[[dt]]) > 0) points(spots.15.g[[dt]], col = adjustcolor("gold", alpha = 0.4), pch = 15, cex = 0.4)
                         abline(h = range(which(!is.na(pw.m[1024,,"white", dt]), arr.ind = T)) + c(69, -68), col = "red")
                     }))
}

####################################################################################################

# MATCHING POINTS BETWEEN SUCCESSIVE ACQUISITIONS                                               ####

# added robustness: default to running over offsets of 100-200px in each direction?, then identify 

plot(0, type = "n", xlim = c(0,2048), ylim = c(0,2048), xlab = "", ylab = "", asp = T)
points(spots.15$"130613", pch = 15, col = "cyan3", cex = 0.4)
points(spots.15$"130701", pch = 15, col = "magenta3", cex = 0.4)

plot(0, type = "n", xlim = c(0,2048), ylim = c(0,2048), xlab = "", ylab = "", asp = T)
pcol <- c("gold", "cyan3", "green3", "blue", "red", "magenta3", "black")
invisible(lapply(c(1:2),
          function(n) {
              points(spots.15[[6 + n]], pch = 15, cex = 0.4, col = adjustcolor(pcol, alpha = 0.4)[n])
          }))

# minimise nearest-neighbour distance to get fit

# get midpoints of clusters
mids <- lapply(lapply(spots.15[sapply(spots.15, nrow) > 0], clump.centres), "[", c("xm", "ym"))

# function to match midpoints & find offset
find.offset <- function(px1, px2, start.par = c(x.adj = 0, y.adj = 0)) {
    
    im.offset <- function(par, mids1, mids2) {
        adj1 <- data.frame(mids1[,1] + par["x.adj"], mids1[,2] + par["y.adj"])
        d <- knnx.dist(mids2, adj1, k = 1)
        d <- d[d <= quantile(d, 0.75)]
        sum(d)
    }
    
    optim(start.par, im.offset, mids1 = px1, mids2 = px2)
}

# function to match pairs of points
find.offset.pairs <- function(px1, px2, start.par = c(x.adj = 0, y.adj = 0)) {
    zz <- find.offset(px1, px2)
    
    adj1 <- data.frame(px1[,1] + zz$par["x.adj"], px1[,2] + zz$par["y.adj"])
    
    i <- knnx.index(px2, adj1, k = 1)
    d <- knnx.dist(px2, adj1, k = 1)
    
    df <- data.frame(x1 = px1[,1], y1 = px1[,2],
                     x2 = px2[i,1], y2 = px2[i,2], dist = d) 
    
    df <- rbind(df,
                setNames(data.frame(NA, NA, px2[!(px2[,1] %in% df$x2 & px2[,2] %in% df$y2),], NA),
                         nm = colnames(df)))

    df <- df[order(df$d),]
    df[duplicated(df[,3:4]), 3:5] <- NA
    df
}

hh <- find.offset.pairs(mids$"141009", mids$"141118")

invisible(apply(hh, 1, function(ll) lines(ll[c(1, 3)], ll[c(2,4)], col = "red", lwd = 2)))

sapply(mids, nrow)

c1 <- find.offset.pairs(mids$"130613", mids$"130701")

c2 <- invisible(lapply(c(4:8), function(i) find.offset.pairs(mids[[i]], mids[[i+1]])))

c2[[1]]

plot(0, type = "n", xlim = c(0,2048), ylim = c(0,2048), xlab = "", ylab = "", asp = T)
pcol <- c("gold", "cyan3", "green3", "blue", "red", "magenta3", "black")
invisible(lapply(c(3:4),
                 function(n) {
                     points(spots.15[[6 + n]], pch = 15, cex = 0.4, col = adjustcolor(pcol, alpha = 0.4)[n])
                 }))

invisible(apply(c.141118, 1, function(ll) lines(ll[c(1, 3)], ll[c(2,4)], col = "red", lwd = 2)))
invisible(apply(c.150108, 1, function(ll) lines(ll[c(1, 3)], ll[c(2,4)], col = "red", lwd = 2)))


c.141118 <- find.offset.pairs(mids$"141118", mids$"141217", start.par = c(x.adj = 150, y.adj = 20))
c.150108 <- find.offset.pairs(mids$"141217", mids$"150108", start.par = c(x.adj = -300, y.adj = -20))
points(spots.15$"150108", pch = 15, cex = 0.4, col = adjustcolor("blue", alpha = 0.4))

os <- find.offset(mids$"141217", mids$"150108", start.par = c(x.adj = -300, y.adj = -20))$par

points(mids$"141217" + os, pch = 15, cex = 0.5)

####################################################################################################

# STEP-BY-STEP PLOTS                                                                            ####

im <- pw.m[,,"white", "141009"]

up <- apply(im[,1024 + c(1:100)], 1, median, na.rm = T)
lp <- apply(im[,1024 + c(0:-100)], 1, median, na.rm = T)

adj <- im
adj[,1025:2048] <- im[,1025:2048] - (up - lp)

{
    res.1.05 <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = 1/5), "[[", 2))
    res.1.15 <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = 1/15), "[[", 2))
    res.1.35 <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = 1/35), "[[", 2))
    
    pixel.image(res.1.05, xlim = x.rng, ylim = y.rng, title = "1/5")    # lose 200px at either border
    pixel.image(res.1.15, xlim = x.rng, ylim = y.rng, title = "1/15")   # lose 68px at either border
    pixel.image(res.1.35, xlim = x.rng, ylim = y.rng, title = "1/35")   # lose 8px at either border
}

res <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = 1/15), "[[", 2))

tr <- res
tr[res > median(res, na.rm = T)] <- median(res, na.rm = T)

sk <- shapeKernel(c(5, 5), type = "disc")
cl <- closing(tr, sk)
cl[is.infinite(cl) | is.na(cl)] <- median(res, na.rm = T)

th <- threshold(cl, method = "literal", level = - sd(res, na.rm = T))

exp <- erode(th, sk)

draw.outlines(which(exp == 0, arr.ind = T))

# plots
{
    x.rng <- c(100,500); y.rng <- c(900,1150); out.height <- 300
    jpeg(paste0(fpath, "ss-raw-image.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(im, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-offset-adj.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(adj, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-loess-res.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(res, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    res.levels <- sd.levels(res)
    jpeg(paste0(fpath, "ss-loess-res-trunc.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(tr, xlim = x.rng, ylim = y.rng, break.levels = res.levels)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-closing.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(cl, xlim = x.rng, ylim = y.rng, break.levels = res.levels)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-thresholded.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(th, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-enlarged.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(exp, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-overlaid.jpg"), height = out.height); {
        par(mar = c(2,2,1,1))
        pixel.image(im, xlim = x.rng, ylim = y.rng)
        draw.outlines(which(exp == 0, arr.ind = T), lwd = 2)
        dev.off()
    }
}

####################################################################################################

# RUN OVER SHADING-CORRECTED IMAGE TO ASSESS EFFECT ON RECONSTRUCTION                           ####

sc.spots <- screen.spots(shading.corrected(pw.m[,,,"141009"]))

pixel.image(shading.corrected(pw.m[,,,"141009"]))
draw.outlines(sc.spots)
draw.outlines(sp.15, lty = 2)
