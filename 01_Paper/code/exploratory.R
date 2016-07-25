
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

load.pixel.means.2048()

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
    jpeg(paste0(fpath, "spot-residuals-160705.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.160705.u)
        plot(which(abs(res.160705.u) > 1000, arr.ind = T), pch = 15, cex = 0.6)
        dev.off()
    }
    length(which(abs(res.160705.u) > 1000))     # 3370
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
    jpeg(paste0(fpath, "spot-residuals-160430.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.160430.u)
        plot(which(abs(res.160430.u) > 1000, arr.ind = T), pch = 15, cex = 0.6)
        dev.off()
    }
    length(which(abs(res.160430) > 1000))     # 2932
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
    jpeg(paste0(fpath, "spot-residuals-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        #pixel.image(res.MCT225)
        plot(which(abs(res.MCT225) > 1000, arr.ind = T), pch = 15, cex = 0.6)
        dev.off()
    }
    length(which(abs(res.MCT225) > 1000))     # 122681
}

####################################################################################################
# TEMPORARY BLOCK - SETTING THRESHOLDING LEVELS                                                 ####

.smoothScatter(which(abs(res.MCT225.u) > 2 * sd(res.MCT225.u, na.rm = T), arr.ind = T), main = "MCT225")
.smoothScatter(which(abs(res.160430.u) > 2 * sd(res.160430.u, na.rm = T), arr.ind = T), main = "160430")
.smoothScatter(which(abs(res.160705.u) > 2 * sd(res.160705.u, na.rm = T), arr.ind = T), main = "160705")

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
