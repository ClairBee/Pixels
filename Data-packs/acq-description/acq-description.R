
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Data-packs/acq-description/fig/"
acq <- import.acq("./Image-data/160430")
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