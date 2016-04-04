
library("IO.Pixels")

library(mixtools)       # fitting normal mixture model to residuals
library(SuppDists)      # fitting Johnson distribution to residuals

# automatically create single array containing all pixel means
load.pixel.means()
bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)

# extract bad pixel map using alternative methods, compare to 'official' map
# 'official' map is from 16-03-14, so apply tests to same data set initially
im <- pw.m[ , , "black", "160314"]

# tmp function to quickly overlay 'official' & found points        
compare.maps <- function() {
    bp.coords <- bp[!duplicated(bp[,c(1,2)]),c(1,2,4)]
    
    if (nrow(bp.coords) > 10000) {
        print("Too many rows. Check data manually first")

    } else {
        
        cols <- c("blue", "red", "gold", "green3", "purple")
        
        plot(bpm[,1:2], xlim = c(1,1996), ylim = c(1, 2100), asp = T)
        points(bp.coords[,1:2], pch = 20, col = cols[bp.coords$type])
        
        legend("top", horiz = T, legend = levels(bp.coords$type), pch = 20, col = cols, bty = "n")
    }
    print(table(bp.coords$type))
    print(table(bp$src, bp$type))
}

###################################################################################################
#                          COMPARING & ASSESSING EFFECTIVENESS OF METHODS                      ####
###################################################################################################
# simple by-eye plots
# remove 'bad' points & check MAD/SD of remainder: should be v similar
# compare to 'official' map: behaviour of points extracted?

###################################################################################################

#                      IDENTIFY BAD PIXELS BY COMPARISON TO BAD PIXEL MAP                      ####
###################################################################################################

{
bpm.params <- xmlToList(xmlParse("./Other-data/Other-images/BadPixelMap_160314/CalibrationParameters.xml"))
    
# get 'bright offset corrected image' using flat field correction
    corr <- 60000 * (m.g - m.b) / (m.w - m.b)
    corr[is.na(corr)] <- 0      # otherwise get NA where FF == D
    corr[corr < 0] <- 0         # negative values not possible (result from different image offset)
    
    pixel.image(corr)
    # still systematic variance per panel, but magnitude much reduced
    s.hist(corr)
    
    bp <- rbind(bp, 
                data.frame(which(corr > (median(corr) * 1.5), arr.ind = T), src = "offset", type = "bright"),
                data.frame(which(corr < (median(corr) * 0.45), arr.ind = T), src = "offset", type = "dim"))
    
    compare.maps()
    
    # addition of 'noisy' category picks up large number of extra (unmatched) points. Discard.
    #bp <- rbind(bp,
    #            data.frame(which(sd.b > (median(sd.b) * 6), arr.ind = T), src = "black", type = "noisy"),
    #            data.frame(which(sd.g > (median(sd.g) * 6), arr.ind = T), src = "grey", type = "noisy"),
    #            data.frame(which(sd.w > (median(sd.w) * 6), arr.ind = T), src = "white", type = "noisy"))


    # use abs. thresholds from calibration file?
    # grey thresholds pick up too many points, but black & white seem ok
    # fills in many of the blanks, plus some additional points
    bp <- rbind(bp,
                data.frame(which(m.b > bpm.params$BlackMaxThreshold, arr.ind = T), src = "black", type = "high"),
                data.frame(which(m.w > bpm.params$WhiteMaxThreshold, arr.ind = T), src = "white", type = "high"),
                data.frame(which(m.b < bpm.params$BlackMinThreshold, arr.ind = T), src = "black", type = "low"),
                data.frame(which(m.w < bpm.params$WhiteMinThreshold, arr.ind = T), src = "white", type = "low"))

    table(bp$src, bp$type)
    
    # trying to establish pattern behind thresholds...
    {
        ecdf.b <- ecdf(m.b)
        ecdf.g <- ecdf(m.g)
        ecdf.w <- ecdf(m.w)
        
        # neither symmetrical nor consistent in quantiles
        ecdf.b(c(bpm.params$BlackMinThreshold, bpm.params$BlackMaxThreshold))    
        ecdf.w(c(bpm.params$WhiteMinThreshold, bpm.params$WhiteMaxThreshold))    
        
        as.numeric(bpm.params$BlackMaxThreshold) / median(m.b); as.numeric(bpm.params$WhiteMaxThreshold) / median(m.w)
        as.numeric(bpm.params$BlackMinThreshold) / median(m.b); as.numeric(bpm.params$WhiteMinThreshold) / median(m.w)
        
        # not symmatric about median/mean. Also not consistent in terms of distance of threshold from median/mean.
        as.numeric(bpm.params$BlackMaxThreshold) - median(m.b);  median(m.b) - as.numeric(bpm.params$BlackMinThreshold)
        as.numeric(bpm.params$WhiteMaxThreshold) - median(m.w);  median(m.w) - as.numeric(bpm.params$WhiteMinThreshold)
        
        as.numeric(bpm.params$BlackMaxThreshold) - mean(m.b);  mean(m.b) - as.numeric(bpm.params$BlackMinThreshold)
        as.numeric(bpm.params$WhiteMaxThreshold) - mean(m.w);  mean(m.w) - as.numeric(bpm.params$WhiteMinThreshold)
        
        # not linearly related to SD of mean values
        (as.numeric(bpm.params$BlackMaxThreshold) - as.numeric(bpm.params$BlackMinThreshold)) / sd(m.b)
        (as.numeric(bpm.params$WhiteMaxThreshold) - as.numeric(bpm.params$WhiteMinThreshold)) / sd(m.w)
    }
    
    # and giving up on that. Easier to try to identify new thresholds, which will pick up same/similar values
}


# per-panel plots of pixels not identified as 'bad': anything obvious?
    # U1: clear streaks of low-valued pixels in white/grey
    # U2-4: some odd-looking values in white channel (spots on screen?)

{
    tmp <- pw.m[,,"black", "160314"];  tmp[as.matrix(bpm[,1:2])] <- NA; sp.b <- subpanels(tmp)
    tmp <- pw.m[,,"grey", "160314"];  tmp[as.matrix(bpm[,1:2])] <- NA; sp.g <- subpanels(tmp)
    tmp <- pw.m[,,"white", "160314"];  tmp[as.matrix(bpm[,1:2])] <- NA; sp.w <- subpanels(tmp)
    
#    for (p in dimnames(sp.b)[[3]]) {
#        pdf(paste0("./Plots/BPM-check/160314-", p, ".pdf"))
            plot(c(sp.b[,,p]), ylim = c(0,65535), pch = 20, main = paste0("Panel ", p), 
                 col = adjustcolor("slateblue1", alpha = 0.5), xlab = "", ylab = "")
            points(c(sp.g[,,p]), ylim = c(0,65535), pch = 20, col = adjustcolor("chartreuse3", alpha = 0.5))
            points(c(sp.w[,,p]), ylim = c(0,65535), pch = 20, col = adjustcolor("gold", alpha = 0.5))
#        dev.off()
    }
}

###################################################################################################

#                               IDENTIFY BAD PIXELS BY QUANTILES                               ####
###################################################################################################

# use ECDF of each image, cut quantiles to get most extreme values
#   - haven't tested yet, but unlikely to prove effective - doesn't account for concavity/convexity


###################################################################################################

#              IDENTIFY BAD PIXELS BY LOWESS-SMOOTHING & FINDING OUTLYING RESIDUALS            ####
###################################################################################################

# use column-wise Loess smoothing (not crossing panel boundaries) to find outliers
# compare results with different smoothing spans: how robust is this?
#   - with small smoothing span, misses large groups of defective pixels (eg. at panel edges)

# initial exploratory meanderings
{
    c.smoo <- lowess.per.column(im, span = 1/15)
    c.res <- im - c.smoo
    c.mad <- mad(c.res)
    
    # find a distribution to fit. Nothing particularly close: kurtosis non-zero.
    # does poor distribution fit matter that much at the extreme tails?
    {
        hist(c.res, breaks = "fd", xlim = c(-1500,1500), prob = T)
        lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(c.res), sd = sd(c.res)), lwd = 3, col = "cornflowerblue")
        lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(c.res), sd = mad(c.res)), lwd = 3, col = "red")
        lines(c(-1500:1500), dsnorm(c(-1500:1500), mean = mean(c.res), sd = mad(c.res), xi = 1.25), lwd = 3, col = "green3")
        lines(c(-1500:1500), dlnorm(c(-1500:1500), mean = log(mean(c.res)), sd = mad(c.res)), lwd = 3, col = "purple")
        
        zz <- normalmixEM(c.res, k = 2)
        
        lines(c(-1500:1500), 
              zz$lambda[1] * dnorm(c(-1500:1500), zz$mu[1], zz$sigma[1]) + 
                  zz$lambda[2] * dnorm(c(-1500:1500), zz$mu[2], zz$sigma[2]),
              lwd = 3, col = "purple")
    }
    # Johnson distribution
    parms <- JohnsonFit(c(c.res), moment = "quant")
    plot(function(x) dJohnson(x, parms), -1500, 1500, add = T, col = "red", lwd = 2)
    
    # In meantime, try various cutoffs based on quantiles...
    c.ul <- quantile(c.res, 0.999) + (1.5 * IQR(c.res))
    c.ll <- quantile(c.res, 0.001) - (1.5 * IQR(c.res))
    
    qJohnson(c(0.001, 0.999), parms)

    col <- 99
    o.plot(im[col,])
    points(which(c.res[col,] > c.ul), im[col, which(c.res[col,] > c.ul)], col = "red")
    points(which(c.res[col,] < c.ll), im[col, which(c.res[col,] < c.ll)], col = "blue")
}

#-----------------------------------------------------------------------------------------------

# initial investigation using quantiles, smoothing with span 1/15
# misses large column-wise runs of bad pixels due to over-smoothing
{
    zz <- rbind.fill(bp.lowess.quantiles(160314, "white", lq = 0.0005, uq = 0.9995),
                     bp.lowess.quantiles(160314, "white", span = 1/7, lq = 0.0005, uq = 0.9995),
                     bp.lowess.quantiles(160314, "white", span = 1/2, lq = 0.0005, uq = 0.9995),
                     bp.lowess.quantiles(160314, "white", span = 1, lq = 0.0005, uq = 0.9995),
                     bp.lowess.johnson(160314, "white", lq = 0.0005, uq = 0.9995),
                     bp.lowess.johnson(160314, "white", span = 1/7, lq = 0.0005, uq = 0.9995),
                     bp.lowess.johnson(160314, "white", span = 1/2, lq = 0.0005, uq = 0.9995),
                     bp.lowess.johnson(160314, "white", span = 1, lq = 0.0005, uq = 0.9995))
    
    # consider first row as possible candidate, apply to all colours
    bp.lq.w <- bp.lowess.quantiles(160314, "white", lq = 0.0005, uq = 0.9995, details = T)
    bp.lq.g <- bp.lowess.quantiles(160314, "grey", lq = 0.0005, uq = 0.9995, details = T)
    bp.lq.b <- bp.lowess.quantiles(160314, "black", lq = 0.0005, uq = 0.9995, details = T)
    
    bp <- reset.bp(160314)
    
    bp <- rbind(bp,
                data.frame(bp.lq.w$low, src = "white", type = "low"),
                data.frame(bp.lq.w$high, src = "white", type = "high"),
                data.frame(bp.lq.g$low, src = "grey", type = "low"),
                data.frame(bp.lq.g$high, src = "grey", type = "high"),
                data.frame(bp.lq.b$low, src = "black", type = "low"),
                data.frame(bp.lq.b$high, src = "black", type = "high"))
    
    compare.maps()
    
    # plot bad pixels identified by 'official' map vs by this approach in each colour
    bp.pts <- data.frame(bp[!duplicated(bp[,1:2]),c(1:2, 4)], 
                         map = factor("CB", levels = c("CB", "Both", "BPM")))
    bp.matched <- merge(bp.pts, bpm, by.x = c("row", "col"), by.y = c("X", "Y"), all = T)
    
    bp.matched$map[is.na(bp.matched$type)] <- "BPM"
    bp.matched$map[!is.na(bp.matched$type) & !is.na(bp.matched$NoisyPixel)] <- "Both"
    bp.matched <- bp.matched[order(bp.matched$map),]
    
    # get sample of 'healthy' pixels to plot alongside. Must be a better way...
    # try sampling X & Y independently, then discarding any points already in list
    all <- merge(merge(x = c(1:1996), y = c(1:1996)), bp.matched[,1:3], by.x = c(1:2), by.y = c(1:2), all.x = T)
    gp <- all[is.na(all$type),]
    
    samp <- gp[sample(1:nrow(gp), nrow(bp.matched), replace = F),1:2]

    plot(bp.matched[bp.matched$map == "BPM",1:2], main = "Pixels on map not identified")
    
    # plot values of all pixels identified
    {
        plot(pw.b[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - black images, 16-03-14")
        points(pw.b[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        
        plot(pw.g[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - grey images, 16-03-14")
        points(pw.g[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        
        plot(pw.w[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - white images, 16-03-14")
        points(pw.w[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("bottomleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
    }
    
    # investigate behaviour of remainder: plot per-panel values of 'healthy' pixels
    b.healthy <- pw.b[,,"160314"]
    b.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    g.healthy <- pw.g[,,"160314"]
    g.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    w.healthy <- pw.w[,,"160314"]
    w.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    sp.b <- subpanels(b.healthy)
    sp.g <- subpanels(g.healthy)
    sp.w <- subpanels(w.healthy)
    
    plot(c(sp.b[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("slateblue1", alpha = 0.5))
    points(c(sp.g[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("chartreuse3", alpha = 0.5))
    points(c(sp.w[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("gold", alpha = 0.5))
    
    # clearly some 'bad' points left - Loess smoothing has failed to pick up corner points in white image
    o.plot(pw.w[1,993:1996,"160314"], ylim = c(0,65535))
    points(lowess(pw.w[1,993:1996,"160314"], f = 1/15), type = "l", col = "red", lwd = 2)
    o.plot(pw.w[2,993:1996,"160314"], add = T, col = adjustcolor("darkblue", alpha = 0.5))
    o.plot(pw.w[3,993:1996,"160314"], add = T, col = adjustcolor("blue", alpha = 0.5))
    o.plot(pw.w[4,993:1996,"160314"], add = T, col = adjustcolor("slateblue", alpha = 0.5))
}

#-----------------------------------------------------------------------------------------------

# now using smoothing span 1
# still fails to identify a number of bad pixels (eg. in column 1)
# 'corridor' of bad values accepted by poorly fitted Loess smoother.
# this approach is not useful on its own (perhaps if validated against row-wise smoothing?)
# however, row-wise smoothing loses a lot of information due to panel edges.
{
    # consider first row as possible candidate, apply to all colours
    bp.lq.w <- bp.lowess.quantiles(160314, "white", lq = 0.0005, uq = 0.9995, span = 1, details = T)
    bp.lq.g <- bp.lowess.quantiles(160314, "grey", lq = 0.0005, uq = 0.9995, span = 1, details = T)
    bp.lq.b <- bp.lowess.quantiles(160314, "black", lq = 0.0005, uq = 0.9995, span = 1, details = T)
    
    bp <- reset.bp(160314)
    
    bp <- rbind(bp,
                data.frame(bp.lq.w$low, src = "white", type = "low"),
                data.frame(bp.lq.w$high, src = "white", type = "high"),
                data.frame(bp.lq.g$low, src = "grey", type = "low"),
                data.frame(bp.lq.g$high, src = "grey", type = "high"),
                data.frame(bp.lq.b$low, src = "black", type = "low"),
                data.frame(bp.lq.b$high, src = "black", type = "high"))
    
    compare.maps()
    
    # plot bad pixels identified by 'official' map vs by this approach in each colour
    bp.pts <- data.frame(bp[!duplicated(bp[,1:2]),c(1:2, 4)], 
                         map = factor("CB", levels = c("CB", "Both", "BPM")))
    bp.matched <- merge(bp.pts, bpm, by.x = c("row", "col"), by.y = c("X", "Y"), all = T)
    
    bp.matched$map[is.na(bp.matched$type)] <- "BPM"
    bp.matched$map[!is.na(bp.matched$type) & !is.na(bp.matched$NoisyPixel)] <- "Both"
    bp.matched <- bp.matched[order(bp.matched$map),]
    
    # get sample of 'healthy' pixels to plot alongside. Must be a better way...
    # try sampling X & Y independently, then discarding any points already in list
    all <- merge(merge(x = c(1:1996), y = c(1:1996)), bp.matched[,1:3], by.x = c(1:2), by.y = c(1:2), all.x = T)
    gp <- all[is.na(all$type),]
    
    samp <- gp[sample(1:nrow(gp), nrow(bp.matched), replace = F),1:2]
    
    plot(bp.matched[bp.matched$map == "BPM",1:2], main = "Pixels on map not identified")
    
    # plot values of all pixels identified
    {
        plot(pw.b[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - black images, 16-03-14")
        points(pw.b[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        
        plot(pw.g[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - grey images, 16-03-14")
        points(pw.g[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        
        plot(pw.w[,,"160314"][as.matrix(samp)], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), main = "Bad pixels - white images, 16-03-14")
        points(pw.w[,,"160314"][as.matrix(bp.matched[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        legend("bottomleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
    }
    
    # investigate behaviour of remainder: plot per-panel values of 'healthy' pixels
    b.healthy <- pw.b[,,"160314"]
    b.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    g.healthy <- pw.g[,,"160314"]
    g.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    w.healthy <- pw.w[,,"160314"]
    w.healthy[as.matrix(bp.matched[,1:2])] <- NA
    
    sp.b <- subpanels(b.healthy)
    sp.g <- subpanels(g.healthy)
    sp.w <- subpanels(w.healthy)
    
    plot(c(sp.b[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("slateblue1", alpha = 0.5))
    points(c(sp.g[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("chartreuse3", alpha = 0.5))
    points(c(sp.w[,,1]), ylim = c(0,65535), pch = 20, col = adjustcolor("gold", alpha = 0.5))
    
    # clearly some 'bad' points left - Loess smoothing has failed to pick up corner points in white image
    o.plot(pw.w[1,993:1996,"160314"], ylim = c(0,65535))
    sm <- lowess(pw.w[1,993:1996,"160314"], f = 1)
    points(sm, type = "l", col = "red", lwd = 2)
    points(sm$x, sm$y + bp.lq.w[[1]]$upper, type = "l", col = "cornflowerblue", lwd = 2)
    points(sm$x, sm$y + bp.lq.w[[1]]$lower, type = "l", col = "cornflowerblue", lwd = 2)
    
    points(which(pw.w[1,993:1996,"160314"] > sm$y + bp.lq.w[[1]]$upper),
           pw.w[1,993:1996,"160314"][which(pw.w[1,993:1996,"160314"] > sm$y + bp.lq.w[[1]]$upper)],
           col = "orange", pch = 20)
    points(which(pw.w[1,993:1996,"160314"] < sm$y + bp.lq.w[[1]]$lower),
           pw.w[1,993:1996,"160314"][which(pw.w[1,993:1996,"160314"] < sm$y + bp.lq.w[[1]]$lower)],
           col = "darkblue", pch = 20)
    
    o.plot(pw.w[2,993:1996,"160314"], add = T, col = adjustcolor("darkblue", alpha = 0.5))
    o.plot(pw.w[3,993:1996,"160314"], add = T, col = adjustcolor("blue", alpha = 0.5))
    o.plot(pw.w[4,993:1996,"160314"], add = T, col = adjustcolor("slateblue", alpha = 0.5))
}


###################################################################################################

#                    FIT FULL PARAMETRIC MODELS AND USE TO IDENTIFY BAD PIXELS                 ####

im <- pw.m[,,"white","160314"]
library(SuppDists)

# circular spot (o = 2) & linear panel model
{
    circ.lm <- fit.circular.lm.poly(im, o = 3)
    circ.res <- matrix(circ.lm$residuals, ncol = 1996)
    
    panel.lm <- fit.panel.lm(circ.res)
    panel.res <- circ.res - panel.lm$fitted.values
    
    mad(panel.res)
    s.hist(panel.res, prob = T)
    lines(c(-1000:1000), dnorm(c(-1000:1000), mean = mean(panel.res), sd = sd(panel.res)), lwd = 3, col = "cornflowerblue")
    lines(c(-1000:1000), dnorm(c(-1000:1000), mean = mean(panel.res), sd = mad(panel.res)), lwd = 3, col = "orange")
    
    parms <- JohnsonFit(c(panel.res), moment = "quant")
    lines(c(-1000:1000), dJohnson(c(-1000:1000), parms), lwd = 3, col = "green3")
    
    # try to transform Johnson distribution to standard normal (for fit checking etc)
    {
        # transform Johnson data to approx.normality (easier to check fit)
        u <- (c(panel.res) - parms$xi) / parms$lambda
        f <- u + sqrt(1+u^2)    # because parms$type == "SU" - see details
        z <- parms$gamma + (parms$delta * log(f))
        
        z <- parms$gamma + (parms$delta * asinh(u))     # alternative definition - check diffs
        
        hist(z, breaks = "fd", prob = T)
        lines(c(-1000:1000)/100, dnorm(c(-1000:1000)/100, mean = 0, sd = 1), lwd = 2, col = "red")
        
        shapiro.test(sample(z, 5000))
        
        # not normal after Johnson transformation. Possibly due to lower values in LH edge panels?
        
        # check residuals per panel
        sp <- subpanels(matrix(z, ncol = 1996))
        
        hist(sp[,,1], breaks = "fd", col = "black", prob = T)
        lines(c(-1000:1000)/100, dnorm(c(-1000:1000)/100, mean = 0, sd = 1), lwd = 2, col = "red")
        
        hist(sp[,,8], breaks = "fd", col = "black", prob = T)
        lines(c(-1000:1000)/100, dnorm(c(-1000:1000)/100, mean = 0, sd = 1), lwd = 2, col = "red")
        
        # transformation hasn't accounted for all of the skewness... Need to look at distribution further
    }
    # transformation hasn't accounted for all of the skewness - need to look at distribution further
    
    

}

# automated version to allow easy comparison of multiple models
{
    zz <- rbind(bp.parametric.quantiles(160314, "white", circ.o = NA, panel.o = NA),
                bp.parametric.quantiles(160314, "white", circ.o = 1, panel.o = NA),
                bp.parametric.quantiles(160314, "white", circ.o = NA, panel.o = 1),
                bp.parametric.quantiles(160314, "white", circ.o = 1, panel.o = 1),
                bp.parametric.quantiles(160314, "white", circ.o = NA, panel.o = 2),
                bp.parametric.quantiles(160314, "white", circ.o = 2, panel.o = NA),
                bp.parametric.quantiles(160314, "white", circ.o = 1, panel.o = 2),
                bp.parametric.quantiles(160314, "white", circ.o = 2, panel.o = 1),
                bp.parametric.quantiles(160314, "white", circ.o = 2, panel.o = 2))
                
    # more of an effect after panel fitting (second operation - unsurprising?)
    # might be worth considering a model in which only panel offsets are fitted, no gradient?
    
    # 4 models have similar residual MAD (< 300), simplest being 2-spot with linear panels
    # compare results over all image sets - is model always applicable?
    bp.w <- bp.parametric.quantiles(160314, "white",  circ.o = 2, panel.o = 1, details = T)
    bp.g <- bp.parametric.quantiles(160314, "grey",  circ.o = 2, panel.o = 1, details = T)
    bp.b <- bp.parametric.quantiles(160314, "black",  circ.o = 2, panel.o = 1, details = T)
    
    comp <- rbind(bp.w$details, bp.g$details, bp.b$details)
    
    bp <- rbind(reset.bp(160314),
                data.frame(bp.w$low, src = bp.w$details$batch, type = "low"),
                data.frame(bp.w$high, src = bp.w$details$batch, type = "high"),
                data.frame(bp.g$low, src = bp.g$details$batch, type = "low"),
                data.frame(bp.g$high, src = bp.g$details$batch, type = "high"),
                data.frame(bp.b$low, src = bp.b$details$batch, type = "low"),
                data.frame(bp.b$high, src = bp.b$details$batch, type = "high"))
    bp <- bp[!duplicated(bp[,1:2]),]

    compare.maps()
    all.bp <- combine.maps(bp, bpm)
    plot.map.comparison(all.bp, 160314)
    
    # far more bad pixels picked up even using simplest parametric approach (4801 pixels vs 2842)
    # still a cluster of red pixels not identified, although many seem to display 'normal' values
    missing.bp <- all.bp[all.bp$map == "old",]
    plot(missing.bp[,1:2]); draw.panels()
    summary(pw.m[,,"white", "160314"][as.matrix(missing.bp[,1:2])])
    
    # plot remaining pixels per panel - do they look ok?
    {
        tmp <- pw.m[,,"black", "160314"];  tmp[as.matrix(all.bp[,1:2])] <- NA; sp.b <- subpanels(tmp)
        tmp <- pw.m[,,"grey", "160314"];  tmp[as.matrix(all.bp[,1:2])] <- NA; sp.g <- subpanels(tmp)
        tmp <- pw.m[,,"white", "160314"];  tmp[as.matrix(all.bp[,1:2])] <- NA; sp.w <- subpanels(tmp)
        
        #for (p in dimnames(sp.b)[[3]]) {
        #pdf(paste0("./Plots/BPM-check/160314-", p, ".pdf"))
        plot(c(sp.b[,,p]), ylim = c(0,65535), pch = 20, main = paste0("Panel ", p), 
             col = adjustcolor("slateblue1", alpha = 0.5), xlab = "", ylab = "")
        points(c(sp.g[,,p]), ylim = c(0,65535), pch = 20, col = adjustcolor("chartreuse3", alpha = 0.5))
        points(c(sp.w[,,p]), ylim = c(0,65535), pch = 20, col = adjustcolor("gold", alpha = 0.5))
        #        dev.off()
    }
}
}



###################################################################################################

#               COMPARE VALUES IN SINGLE PANEL AT A TIME: IS PARAMETRIC MODEL NEEDED?          ####
###################################################################################################

# TRY WITH & WITHOUT CIRCULAR SPOT FITTED, BUT WITHOUT PER-PANEL CORRECTION


###################################################################################################

#                      IDENTIFY BAD PIXELS IN OFFSET CORRECTION WITH PANEL ADJ                 ####
###################################################################################################




###################################################################################################

#                      IDENTIFY BAD PIXELS BY DIFFERENCING & FINDING OUTLIERS                  ####
###################################################################################################
# problem with differencing, rather than Loess-smoothing: effect is dispersed over neighbouring cells.

{# compare differences across row & column, use to identify outliers
    c.diffs <- m.w[,1:1995] - m.w[,2:1996]
    
    # pad by repeating edge rows, to allow 1996x1996 abs. diffs
    c.diffs <- c.diffs[,c(1,1:1995,1995)]
    
    # use sum of abs. neighbouring diffs to return to single point
    c.diff2 <- abs(c.diffs[,1:1996]) + abs(c.diffs[,2:1997])
    
    o.plot(m.w[4,], ylim = c(35000,65535), xlim = c(0,992))
    o.plot(c.diffs[4,] + 40000, add = T, col = adjustcolor("red", alpha = 0.3))
    o.plot(c.diff2[4,] + 40000, add = T, col = adjustcolor("blue", alpha = 0.3))
    
    o.plot(c.diff2[4,])
    mad(c.diff2[4,]); sd(c.diff2[4,])
    points(which(c.diff2[4,] > 2 * mad(c.diff2[4,])))
    
    r.diffs <- m.w[1:1995,] - m.w[2:1996,]
    o.plot(r.diffs[,4])}

