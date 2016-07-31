
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()
models <- read.csv("./Other-data/Gaussian-spots.csv", row.names = 1)

####################################################################################################

# SUMMARY STATISTICS                                                                            ####

df <- data.frame("acq" = sort(rep(dimnames(pw.m)[[4]], 3)),
                 "im" = rep(dimnames(pw.m)[[3]], dim(pw.m)[[4]]),
                 "mean" = c(apply(pw.m, 3:4, mean, na.rm = T)),
                 "median" = c(apply(pw.m, 3:4, median, na.rm = T)),
                 "sd" = c(apply(pw.m, 3:4, sd, na.rm = T)),
                 "mad" = c(apply(pw.m, 3:4, mad, na.rm = T)))

write.csv(df, paste0(fpath, "image-summary.csv"), quote = F, row.names = F)

df <- read.csv(paste0(fpath, "image-summary.csv"), as.is = T)

# convert to .csv format for export
qq <- do.call("rbind", lapply(unique(df$acq), function(dt) unlist(df[df$acq == dt, 3:5])))

write.csv(qq, paste0(fpath, "summary-statistics.csv"), quote = F, row.names = F)

####################################################################################################

# DESCRIPTIVE                                                                                   ####

# plot colour scale
{
    jpeg(paste0(fpath, "image-scale.jpg"), height = 100); {
        par(mar = c(4, 1, 1, 1))
        image.scale(-9:9, c(-9,9), col = sd.colours(), breaks = c(-9,-6,-3,-2,-1,-0.5, 0, 0.5, 1,2,3,6,9))
        title(xlab = expression(paste("Mean value + ", x * sigma), collapse = ""))
        abline(v = 0, lty = 2)
        dev.off()
    }

}

# pixelwise mean images
{
    # 'healthy' detector: 14-10-09 (first images after refurbishment)
    jpeg(paste0(fpath, "pwm-black-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"black", "141009"])
        draw.panels(lty = 3)
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-grey-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "141009"])
        draw.panels(lty = 3)
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-white-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "141009"])
        draw.panels(lty = 3)
        dev.off()
    }
    
    # detector with dark lines: 
    {
        
    }
}

hist.scale <- function(data, sc.offset = -100, xlim = c(min(data, na.rm = T), max(data, na.rm = T)), scale.colours = sd.colours(), scale = sd.levels(data), pch = 15, ...) {
    cl <- cut(xlim[1]:xlim[2], scale)
    points(xlim[1]:xlim[2], rep(sc.offset, length(xlim[1]:xlim[2])), pch = pch, col = scale.colours[cl], ...)
}

# pixelwise mean histograms
{
    jpeg(paste0(fpath, "pwm-black-141009-hist.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"black", "141009"], breaks = "fd", xlim = c(0,65535), xlab = "", ylab = "", main = "")
        #hist.scale(pw.m[,,"black", "141009"], sc.offset = -800)
        dev.off()
    }
    
    jpeg(paste0(fpath, "pwm-black-141009-hist-cropped.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"black", "141009"], breaks = "fd", xlim = c(0,65535), ylim = c(0,30), xlab = "Grey Value observed", ylab = "Frequency", main = "")
        dev.off()
    }
    
    jpeg(paste0(fpath, "pwm-grey-141009-hist.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"grey", "141009"], breaks = "fd", xlim = c(0,65535), xlab = "Grey Value observed", ylab = "Frequency", main = "")
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-grey-141009-hist-cropped.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"grey", "141009"], breaks = "fd", xlim = c(0,65535), ylim = c(0,30), xlab = "Grey Value observed", ylab = "Frequency", main = "")
        dev.off()
    }
    
    jpeg(paste0(fpath, "pwm-white-141009-hist.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"white", "141009"], breaks = "fd", xlim = c(0,65535), xlab = "Grey Value observed", ylab = "Frequency", main = "")
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-white-141009-hist-cropped.jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        hist(pw.m[,,"white", "141009"], breaks = "fd", xlim = c(0,65535), ylim = c(0,30), xlab = "Grey Value observed", ylab = "Frequency", main = "")
        dev.off()
    }
}



####################################################################################################

# PARAMETRIC MODELLING                                                                          ####

# Gaussian spot with constrained centre: 16-07-05
{
    # constrained
    gs.160705 <- gaussian.spot.ls(pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"],
                                  c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                                  x0.target = c(768, 1280), y0.target = c(768, 1280))
    gs.160705.u  <- gaussian.spot.ls(pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"],
                                     c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
    jpeg(paste0(fpath, "Gaussian-spot-constraints-160705.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"])
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.160705$par), add = T, drawlabels = F, lty = 1, lwd = 2)
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.160705.u$par), add = T, drawlabels = F, lty = 2)
        rect(768, 768, 1280, 1280, lty = 3)
        points(gs.160705$par["x0"], gs.160705$par["y0"], pch = 15, cex = 1.5)
        points(gs.160705.u$par["x0"], gs.160705.u$par["y0"], pch = 17, cex = 1.5)
        dev.off()
    }
    
    # create array containing residuals from final fitted spot
    res.160705.u <- pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"] - gaussian.spot.mat(gs.160705.u$par)
    hr.dens.160705 <- r2m(focal(m2r(bpx2im(data.frame(which(abs(res.160705.u) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))),
                                matrix(1/51^2, ncol = 51, nrow = 51)))
    jpeg(paste0(fpath, "spot-residuals-160705.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.160705.u)
        plot(which(abs(res.160705.u) > 1204, arr.ind = T), pch = 15, cex = 0.6)
        points(which(hr.dens.160705 > 0.1, arr.ind = T), pch = 15, cex = 0.6, col = adjustcolor("red", alpha = 0.4))
        dev.off()
    }
    length(which(abs(res.160705.u) > 1204))     # 3370
}

# Gaussian spot with constraint unnecessary: 16-04-30
{
    gs.160430 <- gaussian.spot.ls(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"],
                                  c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                                  x0.target = c(768, 1280), y0.target = c(768, 1280))
    gs.160430.u  <- gaussian.spot.ls(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"],
                                     c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
    jpeg(paste0(fpath, "Gaussian-spot-constraints-160430.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.160430$par), add = T, drawlabels = F, lty = 1, lwd = 2)
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.160430.u$par), add = T, drawlabels = F, lty = 2)
        rect(768, 768, 1280, 1280, lty = 3)
        points(gs.160430$par["x0"], gs.160430$par["y0"], pch = 15, cex = 1.5)
        points(gs.160430.u$par["x0"], gs.160430.u$par["y0"], pch = 17, cex = 1.5)
        dev.off()
    }
    
    # create array containing residuals from final fitted spot
    res.160430 <- pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"] - gaussian.spot.mat(gs.160430$par)
    
    # get local density
    hr.dens.160430 <- r2m(focal(m2r(bpx2im(data.frame(which(abs(res.160430) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))),
                         matrix(1/51^2, ncol = 51, nrow = 51)))
    
    jpeg(paste0(fpath, "spot-residuals-160430.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.160430.u)
        plot(which(abs(res.160430) > 1204, arr.ind = T), pch = 15, cex = 0.6)
        points(which(hr.dens.160430 > 0.1, arr.ind = T), pch = 15, cex = 0.6, col = adjustcolor("cyan3", alpha = 0.4))
        points()
        dev.off()
    }
    length(which(abs(res.160430) > 1204))     # 2932
}

# Gaussian spot with constraint required: MCT225
{
    gs.MCT225 <- gaussian.spot.ls(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"],
                                  c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                                  x0.target = c(768, 1280), y0.target = c(768, 1280))
    gs.MCT225.u  <- gaussian.spot.ls(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"],
                                     c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
    jpeg(paste0(fpath, "Gaussian-spot-constraints-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"])
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.MCT225$par), add = T, drawlabels = F, lty = 1, lwd = 2)
        contour(1:2048, 1:2048, gaussian.spot.mat(gs.MCT225.u$par), add = T, drawlabels = F, lty = 2)
        rect(768, 768, 1280, 1280, lty = 3)
        points(gs.MCT225$par["x0"], gs.MCT225$par["y0"], pch = 15, cex = 1.5)
        points(gs.MCT225.u$par["x0"], gs.MCT225.u$par["y0"], pch = 17, cex = 1.5)
        dev.off()
    }
    
    # create array containing residuals from final fitted spot
    res.MCT225 <- pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] - gaussian.spot.mat(gs.MCT225$par)
    hr.dens.MCT225 <- r2m(focal(m2r(bpx2im(data.frame(which(abs(res.MCT225) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))),
                                matrix(1/51^2, ncol = 51, nrow = 51)))
    jpeg(paste0(fpath, "spot-residuals-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.MCT225)
        plot(which(abs(res.MCT225) > 1204, arr.ind = T), pch = 15, cex = 0.6)
        points(which(hr.dens.MCT225 > 0.1, arr.ind = T), pch = ".", col = adjustcolor("red", alpha = 0.1))
        dev.off()
    }
    length(which(abs(res.MCT225) > 1204))     # 122681
}

####################################################################################################

# SCREEN SPOTS                                                                                  ####

clump.centres <- function(px) {
    
    # clump adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = c(2048, 2048))), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    df <- ddply(xy, .(id), summarise,
                xm = mean(x), ym = mean(y),
                r = ceiling(max(max(x) - min(x), max(y) - min(y)) / 2))
    df
}

# get centres of screen spots
spots <- screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40)
ss <- clump.centres(spots)

# white image with screen spots manually marked
jpeg(paste0(fpath, "pwm-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(pw.m[,,"white", "141009"])
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# shading-corrected image with screen spots manually marked
jpeg(paste0(fpath, "sc-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(shading.corrected(pw.m[,,,"141009"]))
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# position of screen spots in successive acquisitions overplotted
jpeg(paste0(fpath, "spots-overplotted.jpg")); {
    par(mar = c(2,2,1,1))
    plot(screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40),
         col = "cyan3", pch = ".", xlab = "", ylab = "")
    points(screen.spots(pw.m[,,"white", "141118"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("gold", alpha = 0.4), pch = ".")
    points(screen.spots(pw.m[,,"white", "141217"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("magenta3", alpha = 0.4), pch = ".")
    legend("topright", pch = 15, col = c("cyan3", "gold", "magenta3"), bty = "n",
           legend = sapply(dimnames(pw.m)[[4]][1:3], fancy.date))
    dev.off()
}

####################################################################################################
