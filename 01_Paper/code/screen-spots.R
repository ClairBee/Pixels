
# normalise images for easier comparison?

library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()

smoo <- array(apply(pw.m[,,"white", ], 3, lowess.per.column, midline = 1024.5, span = 1/5), dim = c(2048, 2048, 21))
res <- pw.m[,,"white",] - smoo

# can't run over shading-corrected image: most of the affected area has been removed
# (no good for assessing degree of actual damage)

# what about using MA rather than Loess?

####################################################################################################

# PROCEDURE

# reduce upper-vs-lower offset by subtracting offset per column?
# detrend each column using Lowess smoothing (small span preferable b/c lose less at edges)
# morphological approach: closing with circular kernel
# identify depth of screen spot & filter out shallow spots
# enlarge if necessary

# check that spots not corrected by sc are covered by identified pixels?

####################################################################################################

# COMPARE SMOOTHING SPANS

pixel.image(pw.m[,,"white", "141009"])

sp.15 <- screen.spots(pw.m[,,"white", "141009"], smooth.span = 1/15)
draw.outlines(sp)

####################################################################################################

# RUN OVER ALL IMAGES                                                                           ####

spots <- apply(pw.m[,,"white", ], 3, screen.spots)

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

cl <- closing(tr, sk)
cl[is.infinite(cl) | is.na(cl)] <- median(res, na.rm = T)

th <- threshold(cl, method = "literal", level = - sd(res, na.rm = T))
pixel.image(th)

exp <- erode(th, sk)

draw.outlines(which(exp == 0, arr.ind = T))

# plots
{
    x.rng <- c(100,500); y.rng <- c(900,1300)
    jpeg(paste0(fpath, "ss-raw-image.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(im, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-offset-adj.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(adj, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-loess-res.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(res.1.15, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-loess-res-trunc.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(tr, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-closing.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(cl, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-thresholded.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(th, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-enlarged.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(exp, xlim = x.rng, ylim = y.rng)
        dev.off()
    }
    jpeg(paste0(fpath, "ss-overlaid.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(im, xlim = x.rng, ylim = y.rng)
        points(which(exp == 0, arr.ind = T), pch = 0)
        dev.off()
    }
}
