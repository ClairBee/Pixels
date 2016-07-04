
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Data-packs/acq-description/fig/"
acq <- readRDS("./02_Objects/images/pwm-160430.rds")
md <- readRDS("./02_Objects/med-diffs/md-160430.rds")

.median <- hijack(median, na.rm = TRUE)

Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "skyblue", "green3")

write("Images acquired on 16-04-30, main detector", paste0(fpath, "title.txt"))

####################################################################################################

# IMAGE PLOTS                                                                                   ####

jpeg(paste0(fpath, "img-plot-black.jpg")); {
    par(mar = c(2, 2, 1,1))
    pixel.image(acq[,,"black"])
    dev.off()
}

jpeg(paste0(fpath, "img-plot-grey.jpg")); {
    par(mar = c(2, 2, 1,1))
    pixel.image(acq[,,"grey"])
    dev.off()
}

jpeg(paste0(fpath, "img-plot-white.jpg")); {
    par(mar = c(2, 2, 1,1))
    pixel.image(acq[,,"white"])
    dev.off()
}

####################################################################################################

# HISTOGRAMS WITH THRESHOLDS                                                                    ####

# plot histograms to full height
lapply(dimnames(acq)[[3]], function(cc) {
    im <- acq[,,cc]
    jpeg(paste0(fpath, "hist-full-", cc, ".jpg"), height = 240)
    par(mar = c(4, 4, 1, 1))

    hist(im, breaks = "fd", xlab = "Pixelwise mean (grey values)", ylab = "Frequency", main = "")
    med <- median(im, na.rm = T)
    rect(med + (65535 - med)/2, 0, 65535, 60000, col = adjustcolor("red", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(0, 0, med /2, 60000, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
    rect(med /2, 0, med * 0.75, 60000, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
    rect(qJohnson(5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])),
         0, qJohnson(1-5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])), 60000, col = adjustcolor("blue", alpha = 0.3), border = NA)

    dev.off()
})

# also plot cropped histograms (freq <= 30)
lapply(dimnames(acq)[[3]], function(cc) {
    im <- acq[,,cc]
    jpeg(paste0(fpath, "hist-cropped-", cc, ".jpg"), height = 240)
    par(mar = c(4, 4, 1, 1))
    
    hist(im, breaks = "fd", ylim = c(0,30), xlab = "Pixelwise mean (grey values)", ylab = "Frequency", main = "")
    med <- median(im, na.rm = T)
    rect(med + (65535 - med)/2, 0, 65535, 60000, col = adjustcolor("red", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(0, 0, med /2, 60000, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
    rect(med /2, 0, med * 0.75, 60000, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
    rect(qJohnson(5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])),
         0, qJohnson(1-5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])), 60000, col = adjustcolor("blue", alpha = 0.3), border = NA)
    
    dev.off()
})

####################################################################################################

# SCATTERPLOTS WITH THRESHOLDS                                                                  ####    

# get thresholds
th <- apply(acq, 3, function(im) {
    med <- median(im, na.rm = T)
    c(v.dim = med * 0.5, dim = med * 0.75,
      bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
})
nr.lim <- sort(rep(qJohnson(c(5e-04, 1 - 5e-04), JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])), 2))

# smoothed scatterplot of each pair of power levels, showing thresholds
scatterplot.with.thresholds <- function(im, cols = c("black", "grey"), ...) {
    smoothScatter(im[,,cols[1]], im[,,cols[2]], nrpoints = Inf, 
                  xlab = paste0("Value in ", cols[1], " images"),
                  ylab = paste0("Value in ", cols[2], " images"), ...)
    
    # mark bright/hot pixels
    rect(th["v.bright", cols[1]], th["v.bright", cols[2]], 65534, 65534, 
         col = adjustcolor("red", alpha = 0.2), border = NA)
    rect(th["bright", cols[1]], th["bright", cols[2]], 65534, 65534, 
         col = adjustcolor("orange", alpha = 0.2), border = NA)
    abline(v = 65535, col = adjustcolor("magenta3", alpha = 0.3))
    abline(h = 65535, col = adjustcolor("magenta3", alpha = 0.3))
    
    # fit & add line through bright pixels (excl hot pixels)
    abline(line(acq[,,cols[1]][which((acq[,,cols[1]] > th["bright", cols[1]] | acq[,,cols[2]] > th["bright", cols[2]]) & acq[,,cols[2]] < 65535, arr.ind = T)],
                acq[,,cols[2]][which((acq[,,cols[1]] > th["bright", cols[1]] | acq[,,cols[2]] > th["bright", cols[2]]) & acq[,,cols[2]] < 65535, arr.ind = T)]),
           col = adjustcolor("orange", alpha = 0.6), lty = 2)
    abline(line(acq[,,cols[1]][which((acq[,,cols[1]] > th["v.bright", cols[1]] | acq[,,cols[2]] > th["v.bright", cols[2]]) & acq[,,cols[2]] < 65535, arr.ind = T)],
                acq[,,cols[2]][which((acq[,,cols[1]] > th["v.bright", cols[1]] | acq[,,cols[2]] > th["v.bright", cols[2]]) & acq[,,cols[2]] < 65535, arr.ind = T)]),
           col = adjustcolor("red", alpha = 0.6), lty = 3)
    
    # mark dim region
    rect(0,0,th["dim", cols[1]], th["dim", cols[2]], 
         col = adjustcolor("cyan3", alpha = 0.2), border = NA)
    rect(0,0,th["v.dim", cols[1]], th["v.dim", cols[2]], 
         col = adjustcolor("cornflowerblue", alpha = 0.2), border = NA)
    
    # mark non-responsive region
    nr.lim <- sort(rep(qJohnson(c(5e-04, 1 - 5e-04), JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])), 2))
    rect(nr.lim[1], nr.lim[2], nr.lim[3], nr.lim[4], col = adjustcolor("blue", alpha = 0.2), border = NA)
}

jpeg(paste0(fpath, "th-splot-black-v-grey.jpg")) ; {
    par(mar = c(4, 4, 1, 1))
    scatterplot.with.thresholds(acq, cols = c("black", "grey"))
    dev.off()
}
jpeg(paste0(fpath, "th-splot-black-v-white.jpg")) ; {
    par(mar = c(4, 4, 1, 1))
    scatterplot.with.thresholds(acq, cols = c("black", "white"))
    dev.off()
}
jpeg(paste0(fpath, "th-splot-grey-v-white.jpg")) ; {
    par(mar = c(4, 4, 1, 1))
    scatterplot.with.thresholds(acq, cols = c("grey", "white"))
    dev.off()
}

# plot abs. pixel value vs median-differenced value
{
    smoothScatter(acq[,,"black"], md[,,"black"], nrpoints = Inf)
    abline(line(acq[,,"black"][acq[,,"black"] > 10000], md[,,"black"][acq[,,"black"] > 10000]),
           col = adjustcolor("darkred", alpha = 0.5))
    abline(h = 1000, col = "red", lty = 3)
    
    smoothScatter(acq[,,"grey"], md[,,"grey"], nrpoints = Inf)
}

# shading-corrected values vs raw images
sc <- readRDS("./Other-data/Shading-corrections.rds")

jpeg(paste0(fpath, "th-splot-sc-vs-black.jpg")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(acq[3:1998,33:2028,"black"], sc[,,"160430"], nrpoints = Inf, xlim = c(0,65535),
                  ylab = "Shading-corrected", xlab = "Pixelwise mean in black images")
    
    # bright pixels
    rect(th["bright", "black"], -10000, 65535, 60000, col = adjustcolor("orange", alpha = 0.4), border = NA)
    
    dev.off()
}
jpeg(paste0(fpath, "th-splot-sc-vs-grey.jpg")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(acq[3:1998,33:2028,"grey"], sc[,,"160430"], nrpoints = Inf, xlim = c(0,65535),
                  ylab = "Shading-corrected", xlab = "Pixelwise mean in grey images")
    rect(th["bright", "grey"], -10000, 65535, 60000, col = adjustcolor("orange", alpha = 0.4), border = NA)
    rect(nr.lim[1], -10000, nr.lim[3], 60000, col = adjustcolor("blue", alpha = 0.2), border = NA)
    
    dev.off()
}
jpeg(paste0(fpath, "th-splot-sc-vs-white.jpg")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(acq[3:1998,33:2028,"white"], sc[,,"160430"], nrpoints = Inf, xlim = c(0,65535),
                  ylab = "Shading-corrected", xlab = "Pixelwise mean in white images")
    
    rect(th["bright", "white"], -10000, 65535, 60000, col = adjustcolor("orange", alpha = 0.4), border = NA)
    rect(nr.lim[1], -10000, nr.lim[3], 60000, col = adjustcolor("blue", alpha = 0.2), border = NA)
    
    dev.off()
}

####################################################################################################

# SHADING CORRECTION                                                                            ####

sc <- 60000 * (acq[,,"grey"] - acq[,,"black"]) / (acq[,,"white"] - acq[,,"black"])
sc[is.infinite(sc)] <- 0

jpeg(paste0(fpath, "shading-correction.jpg")); {
    par(mar = c(2, 2, 1,1))
    pixel.image(sc)
    dev.off()
}

jpeg(paste0(fpath, "shading-correction-histogram.jpg")); {
    par(mar = c(2, 2, 1,1))
    s.hist(sc, main = "", xlab = "Corrected grey value")
    dev.off()
}

####################################################################################################

# PARAMETRIC DESCRIPTION (QUADRATIC TREND)                                                      ####

# black image
{
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(acq[,,"black"]), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(acq[,,"black"]))
    quad.fitted[which(!is.na(acq[,,"black"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(acq[,,"black"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-black.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted, break.levels = sd.levels(acq[,,"black"]))
        draw.panels()
        dev.off()
    }
        
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-black.txt"))
}

# grey image
{
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(acq[,,"grey"]), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(acq[,,"grey"]))
    quad.fitted[which(!is.na(acq[,,"grey"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(acq[,,"grey"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-grey.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted, break.levels = sd.levels(acq[,,"grey"]))
        draw.panels()
        dev.off()
    }
    
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-grey.txt"))
}

# white image
{
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(acq[,,"white"]), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(acq[,,"white"]))
    quad.fitted[which(!is.na(acq[,,"white"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(acq[,,"white"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-white.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted, break.levels = sd.levels(acq[,,"white"]))
        draw.panels()
        dev.off()
    }
    
    pixel.image(quad.res); draw.panels()
    
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-white.txt"))
}

# shading correction
{
    sc[,,"160430"][is.na(sc[,,"160430"]) | is.infinite(sc[,,"160430"])] <- 0
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(sc[,,"160430"]), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(sc[,,"160430"]))
    quad.fitted[which(!is.na(sc[,,"160430"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(sc[,,"160430"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-sc.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted, break.levels = sd.levels(sc[,,"160430"]))
        draw.panels()
        dev.off()
    }
    
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-grey.txt"))
}

####################################################################################################

# LOCAL BRIGHTNESS                                                                              ####

# consistency of locally bright pixels between black & grey images
{
    hh <- lapply(c(1:20) * 100, 
                 function(lim) length(which(md[,,"black"] > lim & md[,,"grey"] > lim)) / 
                     length(which(md[,,"black"] > lim | md[,,"grey"] > lim)))
    
    plot(c(1:20) * 100, unlist(hh) * 100, pch = 20, ylab = "% points marked in both images", xlab = "threshold")
    abline(h = c(1:10)*10, col = adjustcolor("cyan3", alpha = 0.6))
    abline(h = c(1:10)*10 + 5, col = adjustcolor("cyan3", alpha = 0.2))
    
}