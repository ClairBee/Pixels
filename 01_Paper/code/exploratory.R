
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()
pw.sd <- load.pixel.sds()


####################################################################################################

# SUMMARY STATISTICS                                                                            ####

# table of summary statistics for all images
df <- data.frame("acq" = sort(rep(dimnames(pw.m)[[4]], 3)),
                 "im" = rep(dimnames(pw.m)[[3]], dim(pw.m)[[4]]),
                 "mean" = c(apply(pw.m, 3:4, mean, na.rm = T)),
                 "median" = c(apply(pw.m, 3:4, median, na.rm = T)),
                 "sd" = c(apply(pw.m, 3:4, sd, na.rm = T)),
                 "mad" = c(apply(pw.m, 3:4, mad, na.rm = T)),
                 "skew" = c(apply(pw.m, 3:4, function(im) skewness(c(im), na.rm = T))))
                 

write.csv(df, paste0(fpath, "image-summary.csv"), quote = F, row.names = F)

df <- read.csv(paste0(fpath, "image-summary.csv"), as.is = T)

# reformat for export into paper
qq <- do.call("rbind", lapply(unique(df$acq), function(dt) unlist(df[df$acq == dt, 3:5])))

write.csv(qq, paste0(fpath, "summary-statistics.csv"), quote = F, row.names = F)

####################################################################################################

# DESCRIPTIVE                                                                                   ####

# create standalone plot of colour scale
{
    jpeg(paste0(fpath, "image-scale.jpg"), width = 640, height = 100); {
        par(mar = c(4, 1, 0, 1))
        image.scale(-9:9, c(-9,9), col = sd.colours(), breaks = c(-9,-6,-3,-2,-1,-0.5, 0, 0.5, 1,2,3,6,9))
        title(xlab = expression(paste("Mean value + ", x * sigma), collapse = ""))
        abline(v = c(-9:9), lty = 2)
        abline(v = 0)
        dev.off()
    }

}

# pixelwise mean images - black
{
    jpeg(paste0(fpath, "pwm-black-131122.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"black", "131122"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-black-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"black", "141009"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-black-loan.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"black", "loan"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-black-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"black", "MCT225"])
        draw.panels(p = list(x = panel.edges()$x, y = c(1, 2049)))
        dev.off()
    }

}

# pixelwise mean histograms - black
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

# pixelwise mean images - grey
{
    jpeg(paste0(fpath, "pwm-grey-131122.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "131122"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-grey-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "141009"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-grey-loan.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "loan"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-grey-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", "MCT225"])
        draw.panels(p = list(x = panel.edges()$x, y = c(1, 2049)))
        dev.off()
    }
}

# pixelwise mean images - white
{
    jpeg(paste0(fpath, "pwm-white-131122.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "131122"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-white-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "141009"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-white-loan.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "loan"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "pwm-white-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "MCT225"])
        draw.panels(p = list(x = panel.edges()$x, y = c(1, 2049)))
        dev.off()
    }
}

# offset mean images - white
{
    jpeg(paste0(fpath, "offset-white-131122.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "131122"] - pw.m[,,"black", "131122"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "offset-white-141009.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "141009"] - pw.m[,,"black", "141009"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "offset-white-loan.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "loan"] - pw.m[,,"black", "loan"])
        draw.panels()
        dev.off()
    }
    jpeg(paste0(fpath, "offset-white-MCT225.jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "MCT225"] - pw.m[,,"black", "MCT225"])
        draw.panels()
        dev.off()
    }
}

# pixelwise SD
{
    pixel.image(pw.sd[,,"grey","141009"])
    draw.panels(lty = 3)
    .smoothScatter(pw.sd[,,"grey", "141009"], pw.m[,,"grey", "141009"])
}

####################################################################################################

# PARAMETRIC MODELLING                                                                          ####

models <- read.csv("./Other-data/Gaussian-spots.csv", row.names = 1)

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
                r = ceiling(max(max(x) - min(x), max(y) - min(y)) / 2),
                size = length(id))
    df
}

spots <- apply(pw.m[,,"white", ], 3, screen.spots, enlarge = T, ignore.edges = 40)

# get centres of screen spots
spots <- screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40)
ss <- clump.centres(spots)

# "131002" "131122" "140128" "140129" "141009" "141118" "141217" "150108" "150113" "150126" "150529" "150730" "150828" "151015" "160314" "160430" "160705" "loan"   "MCT225"
spots[["131002"]] <- screen.spots(pw.m[,,"white", "131002"], enlarge = T, ignore.edges = 40)


# white image with screen spots manually marked
jpeg(paste0(fpath, "pwm-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(pw.m[,,"white", "141009"], xlim = c(0,1024), ylim = c(0,1024))
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# shading-corrected image with screen spots manually marked
jpeg(paste0(fpath, "sc-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(shading.corrected(pw.m[,,,"141009"]), xlim = c(0,1024), ylim = c(0,1024))
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# position of screen spots in successive acquisitions overplotted
jpeg(paste0(fpath, "spots-overplotted.jpg")); {
    par(mar = c(2,2,1,1))
    plot(screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40),
         col = "cyan3", pch = ".", xlab = "", ylab = "", xlim = c(0,1024), ylim = c(0,1024))
    points(screen.spots(pw.m[,,"white", "141118"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("gold", alpha = 0.4), pch = ".")
    points(screen.spots(pw.m[,,"white", "141217"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("magenta3", alpha = 0.4), pch = ".")
    legend("bottomleft", pch = 15, col = c("cyan3", "gold", "magenta3"), bty = "n",
           legend = sapply(dimnames(pw.m)[[4]][1:3], fancy.date))
    dev.off()
}

####################################################################################################

####################################################################################################

# PARTIAL-COLUMN DEFECTS                                                                        ####

ll <- find.lines(pw.m[,,"grey", "131122"], dim.lines = T)
ll.px <- summarise.lines(ll)

n.plot <- function(cc, im, add = F, ...) {
    
    if (add) {
        lines(im[cc, ] - im[cc - 1, ], ...)
    } else {
        plot(im[cc, ] - im[cc - 1, ], type = "l", xlab = "", ylab = "", main = paste0(cc, "-", cc-1), ...)
        abline(h = 0, col = "red")
    }
}

# dark column (1302 1526 1893 1398 1998  610 1327   88 1994)
n.plot(1301, pw.m[,,"grey", "131122"], xlim = c(1025, 2048), col = "sky")

plot(pw.m[1301,,"grey", "131122"], xlim = c(1025, 2048), type = "l", col = "skyblue")
lines(pw.m[1303,,"grey", "131122"], col = "cyan3")
lines(pw.m[1302,,"grey", "131122"])


# dim & bright (232)
n.plot(232, pw.m[,,"grey", "131122"], xlim = c(1025, 2048))

plot(pw.m[231,,"grey", "131122"], type = "l", col = "skyblue", xlim = c(1025, 2048))
lines(pw.m[233,,"grey", "131122"], col = "cyan3")
lines(pw.m[232,,"grey", "131122"])






####################################################################################################

# FULL-COLUMN DEFECTS                                                                           ####

jpeg(paste0(fpath, "Column-defect-407.jpg"), height = 240); {
    par(mar = c(2,2,1,1))
    plot(pw.m[407,,"white", "loan"], type = "l", xlim = c(0,1024), ylim = c(4500, 5500), col = "gold", xlab = "", ylab = "")
    lines(pw.m[407,,"grey", "loan"], col = "green3")
    lines(pw.m[407,,"black", "loan"])
    abline(v = 1024.5, col = "red")
    
    legend("bottomleft", col = c("black", "green", "gold", "red"), legend = c("Black", "Grey", "White", "Panel midline"), lty = 1, bty = "n")
    dev.off()
}
apply(pw.m[407,200:1024,, "loan"], 2, sd)
apply(pw.m[408,200:1024,, "loan"], 2, sd)

jpeg(paste0(fpath, "Column-defect-1766.jpg"), height = 240); {
    par(mar = c(2,2,1,1))
    plot(pw.m[1302,,"white", "130613"], type = "l", xlim = c(1025, 2048), ylim = c(4500, 5500), col = "gold", xlab = "", ylab = "")
    lines(pw.m[1302,,"grey", "130613"], col = "green3")
    lines(pw.m[1302,,"black", "130613"])
    abline(v = 1024.5, col = "red")
    
    legend("bottomleft", col = c("black", "green", "gold", "red"), legend = c("Black", "Grey", "White", "Panel midline"), lty = 1, bty = "n")
    dev.off()
}

apply(pw.m[1766,200:1024,, "131122"], 2, sd, na.rm = T)
apply(pw.m[1767,200:1024,, "131122"], 2, sd, na.rm = T)

apply(pw.m[1302,1025:1848,, "130613"], 2, sd, na.rm = T)
apply(pw.m[1303,1025:1848,, "130613"], 2, sd, na.rm = T)

o.plot(pw.m[1100:1500,1600,"white","130613"])

o.plot(pw.m[1302,,"white","130613"])

####################################################################################################


####################################################################################################

####################################################################################################

# READOUT DARK (PER SUBPANEL)                                                                   ####

fpath <- "./Image-plots/subpanels/"
library(spatial)

surface.trend <- function(im, order = 2) {
    gv <- setNames(melt(im), nm = c("x", "y", "z"))
    s.ls <- surf.ls(order, gv[!is.na(gv$z),])
    trmat(s.ls, 1, 2048, 1, 2048, 2047)$z
}

subpanel.lm <- function(im) {
    
    # fit linear model to each subpanel
    apply(array(im, dim = c(128, 16, 1024, 2)), c(2, 4),
          function(s) lm(z ~ x + y, data = setNames(melt(s), nm = c("x", "y", "z"))))
}

predict.subpanels <- function(sp.fitted) {
    
    sp.template <- setNames(melt(array(dim = c(128, 1024))), nm = c("x", "y", "z"))
    
    sp.predicted <- lapply(sp.fitted, predict, newdata = sp.template)
    array(abind(abind(lapply(sp.predicted[1:16], array, dim = c(128, 1024)), along = 1.5),
                abind(lapply(sp.predicted[17:32], array, dim = c(128, 1024)), along = 1.5), along = 4),
          dim = c(2048, 2048))
}

subpanel.gradients <- function(sp.fitted) {
    
    # create matrix of fitted coefficients
    models <- do.call("rbind", lapply(sp.fitted, "[[", "coefficients"))
    
    setNames(data.frame(models, models[,2] * 128 + models[,3]*1024),
             nm = c("intercept", "x", "y", "g"))
}

#----------------------------------------------------------------------------------------
# fit quadratic surface to all images in 'main seqence'
qs <- array(apply(pw.m[,,"black", 7:19], 3, surface.trend), dim = c(2048, 2048, 13))
adj <- pw.m[,,"black", 7:19] - qs

saveRDS(adj, paste0(fpath, "quadratic-surface-removed.rds"))

sp.models <- list()
sp.models$"141009" <- subpanel.lm(adj[,,"141009"])
sp.models$"141118" <- subpanel.lm(adj[,,"141118"])
sp.models$"141217" <- subpanel.lm(adj[,,"141217"])
sp.models$"150108" <- subpanel.lm(adj[,,"150108"])
sp.models$"150113" <- subpanel.lm(adj[,,"150113"])
sp.models$"150126" <- subpanel.lm(adj[,,"150126"])
sp.models$"150529" <- subpanel.lm(adj[,,"150529"])
sp.models$"150730" <- subpanel.lm(adj[,,"150730"])
sp.models$"150828" <- subpanel.lm(adj[,,"150828"])
sp.models$"151015" <- subpanel.lm(adj[,,"151015"])
sp.models$"160314" <- subpanel.lm(adj[,,"160314"])
sp.models$"160430" <- subpanel.lm(adj[,,"160430"])
sp.models$"160705" <- subpanel.lm(adj[,,"160705"])

saveRDS(sp.models, paste0(fpath, "sp-fitted.rds"))

g <- abind(lapply(lapply(sp.models, subpanel.gradients), as.matrix), along = 3)

pdf(paste0(fpath, "Sp-gradient-plots.pdf")); {
    par(mfrow = c(4,2), mar = c(2,2,3,1))
    
    matplot(g[1:16,"g",], type = "l", main = "Gradient (lower)", xlab = "", ylab = "")
    matplot(g[17:32,"g",], type = "l", main = "Gradient (upper)", xlab = "", ylab = "")
    
    matplot(g[1:16,"intercept",], type = "l", main = "Intercept (lower)", xlab = "", ylab = "")
    matplot(g[17:32,"intercept",], type = "l", main = "Intercept (lower)", xlab = "", ylab = "")
    
    matplot(g[1:16,"x",], type = "l", main = "X-gradient (lower)", xlab = "", ylab = "")
    matplot(g[17:32,"x",], type = "l", main = "X-gradient (lower)", xlab = "", ylab = "")
    
    matplot(g[1:16,"y",], type = "l", main = "Y-gradient (lower)", xlab = "", ylab = "")
    matplot(g[17:32,"y",], type = "l", main = "Y-gradient (lower)", xlab = "", ylab = "")
    
    dev.off()
}


####################################################################################################

# READOUT DARK (WHOLE IMAGE)                                                                    ####

fpath <- "./Image-plots/Diffs-between-acquisitions/"

diffs.old <- pw.m[,,"black", c("130701", "131002", "131122")] - pw.m[,,"black", c("130613", "130701", "131002")]
diffs <- pw.m[,,"white", 8:19] - pw.m[,,"white", 7:18]


# plots of changes in old images
invisible(lapply(dimnames(diffs.old)[[3]], 
                function(dd) {
                    bmp(paste0(fpath, "diff-black-", dd, ".bmp"))
                    par(mar = c(2,2,1,1))
                    pixel.image(diffs.old[,,dd])
                    draw.panels(lty = 3)
                    dev.off()
                }))

# plots of changes in new images
invisible(lapply(dimnames(diffs)[[3]], 
                 function(dd) {
                     bmp(paste0(fpath, "diff-white-", dd, ".bmp"))
                     par(mar = c(2,2,1,1))
                     pixel.image(diffs[,,dd])
                     draw.panels(lty = 3)
                     dev.off()
                 }))

# change in dark lines - they're actually LESS changeable than the rest of the panel.
plot(pw.m[1314,,"black", "131122"], xlim = c(0,1024), type = "l", xlab = "", ylab = "",
       ylim = range(pretty(pw.m[1314,1:1024,"black",1:4])))
lines(pw.m[1314,,"black", "131002"], col = "blue")
lines(pw.m[1314,,"black", "130701"], col = "forestgreen")
lines(pw.m[1314,,"black", "130613"], col = "green3")

# change in corner value
plot(pw.m[2024,,"black", "131122"], xlim = c(1025, 2048), type = "l", xlab = "", ylab = "",
     ylim = range(pretty(pw.m[2024,1025:2048,"black",1:4])))
lines(pw.m[2024,,"black", "131002"], col = "blue")
lines(pw.m[2024,,"black", "130701"], col = "forestgreen")
lines(pw.m[2024,,"black", "130613"], col = "green3")

####################################################################################################

# EDGE EFFECTS                                                                                  ####

# divide image into edge region & border of same size from within panel, compare

box.strip <- function(im, start, width) {
    xo <- c(start, ncol(im) - start)
    xi <- xo + c(width, - width)
    
    yo <- c(start, nrow(im) - start)
    yi <- yo + c(width, - width)
    
    c(im[xo[1]:xi[1], yo[1]:yo[2]], 
      im[xo[2]:xi[2], yo[1]:yo[2]],
      im[xi[1]:xi[2], yo[1]:yi[1]],
      im[xi[1]:xi[2], yo[2]:yi[2]])
}

st <- box.strip(pw.m[,,"black", "160430"], start = 1, width = 60)
st.inner <- box.strip(pw.m[,,"black", "160430"], start = 61, width = 64)
st.inmost <- box.strip(pw.m[,,"black", "160430"], start = 126, width = 69)
st.central <- box.strip(pw.m[,,"black", "160430"], start = 196, width = 75)


hist(st, breaks = "fd", xlim = c(0,10000))
hist(st.inner, breaks = "fd", border = adjustcolor("cyan3", alpha = 0.3), add = T)
hist(st.inmost, breaks = "fd", border = adjustcolor("gold", alpha = 0.3), add = T)
hist(st.central, breaks = "fd", border = adjustcolor("orange", alpha = 0.3), add = T)

l <- list(st, st.inner, st.inmost, st.central)
round(unlist(lapply(l, mean, na.rm = T)), 1)
round(unlist(lapply(l, sd, na.rm = T)), 1)

ss <- lapply(c(1:10)*10, function(ww) sd(box.strip(pw.m[,,"black", "160430"], start = 1, width = ww), na.rm = T))
ss.g <- lapply(c(1:10)*10, function(ww) sd(box.strip(pw.m[,,"grey", "160430"], start = 1, width = ww), na.rm = T))
ss.w <- lapply(c(1:10)*10, function(ww) sd(box.strip(pw.m[,,"white", "160430"], start = 1, width = ww), na.rm = T))

plot(c(1:10)*10, ss.w, pch = 20, xlab = "border width", ylab = "SD")

ss.b <- apply(pw.m[,,"black",], 3, 
              function(im) {
                  unlist(lapply(c(1:10)*10, function(ww) sd(box.strip(im, start = 1, width = ww), na.rm = T) / sd(im, na.rm = T)))
              })

matplot(c(1:10)*10, ss.b, type = "l", xlab = "Border width", ylab = "Standard deviation")

ss.g <- apply(pw.m[,,"grey",], 3, 
              function(im) {
                  unlist(lapply(c(1:10)*10, function(ww) sd(box.strip(im, start = 1, width = ww), na.rm = T) / sd(im, na.rm = T)))
              })
matplot(c(1:10)*10, ss.g, type = "l", xlab = "Border width", ylab = "Regional vs global SD")
abline(h = 1, col = "yellow")

ss.w <- apply(pw.m[,,"white",], 3, 
              function(im) {
                  unlist(lapply(c(1:10)*10, function(ww) sd(box.strip(im, start = 1, width = ww), na.rm = T) / sd(im, na.rm = T)))
              })
matplot(c(1:10)*10, ss.w, type = "l", xlab = "Border width", ylab = "Regional vs global SD")
abline(h = 1, col = "yellow")