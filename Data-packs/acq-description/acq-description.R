
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Data-packs/acq-description/fig/"
#md <- readRDS("./02_Objects/med-diffs/md-160430.rds")

.median <- hijack(median, na.rm = TRUE)

Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "skyblue", "green3")

acq <- readRDS("./02_Objects/images/pwm-140128.rds")
write("Old data, after refurbishment: images acquired on 14-01-28, main detector", paste0(fpath, "title.txt"))

####################################################################################################

#### images already processed ####                                                              ####

# "./02_Objects/images/pwm-141009.rds"; "Images acquired on 14-10-09, main detector"
# "./02_Objects/images/pwm-160430.rds"; "Images acquired on 16-04-30, main detector"
# "./02_Objects/old-data/pwm-131122.rds"; "Old data: images acquired on 13-11-22, main detector"
# "./02_Objects/old-data/pwm-140128.rds"; "Old data after refurbishment: images acquired on 14-01-28, main detector"
# "./02_Objects/images/pwm-MCT225.rds"; "Images acquired on 16-07-02, Aylesbury detector MCT225"
# "./02_Objects/images/pwm-160705.rds"; "Images acquired on 16-07-05, main detector"

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

# also plot cropped histograms (freq <= 30, outside 40px removed)
lapply(dimnames(acq)[[3]], function(cc) {
    im <- acq[ , , cc]
    jpeg(paste0(fpath, "hist-cropped-", cc, ".jpg"), height = 240); {
        par(mar = c(4, 4, 1, 1))
        
        hist(im, breaks = "fd", ylim = c(0,30), col = "darkgrey", border = "darkgrey",
             xlab = "Pixelwise mean (grey values)", ylab = "Frequency", main = "")
        hist(im[41:2008, 41:2008], breaks = "fd", add = T, col = "black")
        med <- median(im, na.rm = T)
        rect(med + (65535 - med)/2, 0, 65535, 60000, col = adjustcolor("red", alpha = 0.3), border = NA)
        rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
        rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
        rect(0, 0, med /2, 60000, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
        rect(med /2, 0, med * 0.75, 60000, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
        rect(qJohnson(5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])),
             0, qJohnson(1-5e-04, JohnsonFit(acq[,,"black"][!is.na(acq[,,"black"])])), 60000, col = adjustcolor("blue", alpha = 0.3), border = NA)
        
        dev.off()
    }
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
    grey.adj <- acq[,,"grey"] - acq[,,"black"]
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(grey.adj), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(acq[,,"grey"]))
    quad.fitted[which(!is.na(acq[,,"grey"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(acq[,,"grey"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-grey.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted)
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
    white.adj <- acq[,,"white"] - acq[,,"black"]
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(white.adj), nm = c("x", "y", "gv")))
    
    quad.fitted <- quad.res <- array(dim = dim(acq[,,"white"]))
    quad.fitted[which(!is.na(acq[,,"white"]), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(acq[,,"white"]), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-white.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted)
        draw.panels()
        dev.off()
    }
    
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-white.txt"))
}

# shading correction
{
    quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
                   setNames(melt(sc), nm = c("x", "y", "gv"))[!is.na(sc),])
    
    quad.fitted <- quad.res <- array(dim = dim(sc))
    quad.fitted[which(!is.na(sc), arr.ind = T)] <- quad.lm$fitted.values
    quad.res[which(!is.na(sc), arr.ind = T)] <- quad.lm$residuals
    
    jpeg(paste0(fpath, "quad-trend-fitted-sc.jpeg")); {
        par(mar = c(2, 2, 1, 1))
        pixel.image(quad.fitted, break.levels = sd.levels(sc))
        draw.panels()
        dev.off()
    }
    
    write(paste(c("$", apply(cbind(c("", c("", "+")[(coef(quad.lm)[-1] > 0) + 1]),
                                   format(round(coef(quad.lm), 4), scientific = F),
                                   c("", "x", "y", "x^2", "y^2", "xy")),
                             1, paste, collapse = ""), "$"), collapse = ""),
          paste0(fpath, "quad-coef-sc.txt"))
}

####################################################################################################

# PREDICT WHITE/SC VALUES FROM GREY & BLACK IMAGES                                              ####

# white values predicted from black & grey
{
    sc[is.na(sc)] <- 0
    # data frame of all variables for active region of image
    df <- setNames(data.frame(melt(acq[,,"black"]), 
                              melt(acq[,,"grey"]),
                              melt(acq[,,"white"]),
                              melt(sc))[,c("X1", "X2", "value", "value.1", "value.2", "value.3")],
                   nm = c("x", "y", "b", "g", "w", "sc"))
    df <- df[!is.na(df$b),]
    
    # fit linear model to healthy px only, check line through fitted/actual
    w.lm <- lm(w ~ b * g, data = df)                         # 236.2960      0.9949
    
    df$w.fv <- w.lm$fitted.values
    df$w.res <- w.lm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(w.lm)$adj.r.squared, 3), "; ",
                 "RMSE ", round(summary(w.lm)$sigma, 2)),
          paste0(fpath, "fitted-wv-all.txt"))
    
    pdf(paste0(fpath, "fitted-wv-all-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$w, df$w.fv, xlim = c(0,65535), ylim = c(0,65535),
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Observed white value", ylab = "Fitted white value")
        abline(line(df$w, df$w.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
        dev.off()
    }
}

# shading-corrected values predicted from black & grey
{
    # fit linear model, check line through fitted/actual
    sc.lm <- lm(sc ~ b * g, data = df)                         # 236.2960      0.9949
    
    df$sc.fv <- sc.lm$fitted.values
    df$sc.res <- sc.lm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(sc.lm)$adj.r.squared, 3), "; ",
                 "RMSE ", round(summary(sc.lm)$sigma, 2)),
          paste0(fpath, "fitted-sc-all.txt"))
    
    pdf(paste0(fpath, "fitted-sc-all-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$sc, df$sc.fv,
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Observed value", ylab = "Fitted value")
        abline(line(df$sc, df$sc.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
        dev.off()
    }
}
