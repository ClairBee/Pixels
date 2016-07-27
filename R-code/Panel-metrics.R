
# IDEAS FOR PANEL METRICS

library("IO.Pixels"); library("CB.Misc")

models <- read.csv("./Other-data/Gaussian-spots.csv", row.names = 1)
load.pixel.means.2048()
fpath <- "./Notes/Panel-metrics/fig/" 
.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

# actual bad pixels vs bad pixel map
# proportion of screen covered by bad pixels (compare to official thresholds?)
# circularity of spot response
# consistency of subpanels (gradient/offset)
# screen spots
# power required to obtain given mean value
# SD of observed values
# pixelwise SD of observed values

####################################################################################################

# should there be a single x-gradient instead of per-subpanel?

####################################################################################################

# FIT CENTRE-CONSTRAINED GAUSSIAN SPOT                                                          ####

# fit Gaussian spot to all acquisitions (around 1 min per acquisition)
# centre constrained to lie in central region

os <- array(apply(pw.m, 4, function(acq) acq[,,"grey"] - acq[,,"black"]), dim = dim(pw.m[,,1,]))
os.gs <- apply(os, 3,
               function(im) gaussian.spot.ls(im,
                                             c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                                             x0.target = c(768, 1280), y0.target = c(768, 1280)))

models <- data.frame(do.call(rbind, lapply(os.gs, "[[", 1)))
row.names(models) <- dimnames(pw.m)[[4]]

# manually adjust model for 16-07-05: spot was not centred
models["160705",1:5] <- gs.160705.u$par

write.csv(models, "./Other-data/Gaussian-spots.csv", row.names = T)


####################################################################################################

# RESIDUALS                                                                                     ####

fv <- array(apply(models, 1, gaussian.spot.mat, arr.dim = c(2048, 2048)),
            dim = c(2048, 2048, nrow(models)),
            dimnames = list(NULL, NULL, rownames(models)))

res <- abind(sapply(rownames(models), function(dt) pw.m[,,"grey", dt] - pw.m[,,"black", dt] - fv[,, dt],
                    simplify = F), along = 3)

models$rmse <- apply(res, 3, sd, na.rm = T)   # SD of residuals
write.csv(round(models, 3), paste0(fpath, "spot-models.csv"))

pixel.image(fv[,,"loan"])
pixel.image(res[,,"MCT225"])
px <- which(abs(res[,,"MCT225"]) > 2 * sd(res[,,"MCT225"], na.rm = T), arr.ind = T)

models$gt.2rmse <- apply(res, 3, function(r) length(which(abs(r) > 2 * sd(r, na.rm = T))))
models$gt.1000 <- apply(res, 3, function(r) length(which(abs(r) > 1000)))
models$gt.1500 <- apply(res, 3, function(r) length(which(abs(r) > 1500)))
models$gt.1204 <- apply(res, 3, function(r) length(which(abs(r) > 1204)))

plot(models$rmse, models$gt.1204)
mean(models$rmse[7:18])

####################################################################################################

# SUMMARY PLOTS OF ALL IMAGES                                                                   ####

for (dt in row.names(models)) {
    jpeg(paste0(fpath, "offset-image-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt])
        contour(1:2048, 1:2048, fv[,,dt], add = T, drawlabels = F)
        dev.off()
    }
    jpeg(paste0(fpath, "model-transects-h-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot(pw.m[,1024,"grey", dt] - pw.m[,1024,"black", dt], type = "l", xlab = "", ylab = "")
        lines(fv[,1024,dt], col = "red", lwd = 2)
        dev.off()
    }
    jpeg(paste0(fpath, "model-transects-v-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot(pw.m[1024,,"grey", dt] - pw.m[1024,,"black", dt], type = "l", xlab = "", ylab = "")
        lines(fv[1024, , dt], col = "red", lwd = 2)
        dev.off()
    }
    jpeg(paste0(fpath, "residuals-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(res[,,dt])
        contour(1:2048, 1:2048, fv[,,dt], add = T, drawlabels = F)
        dev.off()
    }
    jpeg(paste0(fpath, "large-residuals-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        plot(which(abs(res[,,dt]) > 1000, arr.ind = T),
             xlim = c(0,2048), ylim = c(0,2048), col = "cyan3", pch = 15, cex = 0.6)
        points(which(abs(res[,,dt]) > 1204, arr.ind = T),
             xlim = c(0,2048), ylim = c(0,2048), col = "black", pch = 15, cex = 0.6)
        contour(1:2048, 1:2048, fv[,,dt], add = T, drawlabels = F)
        dev.off()
    }
}

# refit 16-07-05 without x0, y0 constraint (spot genuinely looks off-centre)
{
    gs.u <- gaussian.spot.ls(pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"],
                             c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
    gs.160705 <- gaussian.spot.mat(gs.u$par)
    pixel.image(pw.m[,,"grey", "160705"] - pw.m[,,"black", "160705"])
    contour(1:2048, 1:2048, gs.160705, add = T, drawlabels = F)
    
    gs.u$par
    models["160705",1:5] <- gs.u$par
}

####################################################################################################

# DENSE REGIONS OF DEFECTS                                                                      ####

qq <- bpx2im(data.frame(which(abs(res[,,"MCT225"]) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))
qq <- r2m(focal(m2r(qq), matrix(1/51^2, ncol = 51, nrow = 51)))

pixel.image(qq)
hist(qq, breaks = "fd")

hres.density <- apply(res, 3, 
                      function(im) {
                          tmp <- bpx2im(data.frame(which(abs(im) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))
                          r2m(focal(m2r(tmp), matrix(1/51^2, ncol = 51, nrow = 51), na.rm = T))
})
hr <- array(hres.density, dim = c(2048, 2048, 21))
dimnames(hr) <- dimnames(pw.m[,,"black",])
saveRDS(hr, "./Other-data/spot-residual-density.rds")


pixel.plot(which(hr[,,"loan"] > 0.1, arr.ind = T))

hr.loan <- bpx2im(data.frame(which(abs(res[,,"loan"]) > 1204, arr.ind = T), type = 1), im.dim = c(2048, 2048))
hr.loan <- r2m(focal(m2r(hr.loan), matrix(1/51^2, ncol = 51, nrow = 51)))
pixel.plot(which(hr.loan > 0.1, arr.ind = T))

pixel.plot(which(abs(res[,,"loan"]) > 1204, arr.ind = T))
pixel.image(hr.loan)
points(which(abs(res[,,"loan"]) > 1204, arr.ind = T), pch = ".")
hr.loan[1520, 69]

lapply(dimnames(hr)[[3]], 
       function(dt) {
           pdf(paste0(fpath, "hd-res-", dt, ".pdf"))
           par(mar = c(2,2,1,1))
           pixel.plot(which(hr[,,dt] > 0.1, arr.ind = T))
           dev.off()
       })

####################################################################################################

bpm <- as.data.frame(t(as.matrix(as.data.frame(xmlToList(xmlParse("./Other-data/Old-data/140129/BadPixelMap.bpm.xml"))$BadPixels))), stringsAsFactors = F)
bpm[,1:3] <- sapply(bpm[,1:3], as.numeric)

bpm$y <- 1999 - bpm$Y + 20
bpm$x <- bpm$X + 25

# measured from top to bottom, with no offset

pixel.plot(bpm[,c("x", "y")], col = "red", cex = 0.4, xlim = c(700,800), ylim = c(1000, 1100))
points(which(pw.m[,,"grey","140129"] - pw.m[,,"black", "140129"] < 10000, arr.ind = T), cex = 0.4, pch = 15)
draw.panels(col = "gold")

o.plot(pw.m[769,,"grey", "140129"], xlim = c(0, 1024))
abline(v = 1024.5, col = "red")

####################################################################################################

# SUBPANEL SNR                                                                                  ####

# may be useful as tool for assessing impact of removal of bad pixels.
sp <- array(pw.m[,,, "160430"], dim = c(128, 16, 1024, 2, 3),
            dimnames = list(NULL, NULL, NULL, c("L", "U"), c("black", "grey", "white")))
sp.means <- apply(sp, c(1,2,4,5), mean, na.rm = T)
sp.sds <- apply(sp, c(1,2,4,5), sd, na.rm = T)

snr.per.col <- sp.means / sp.sds

# SNR of shading-corrected image
{
    sc.sp <- array(shading.corrected(pw.m[,,,"160430"]), dim = c(128, 16, 1024, 2),
                   dimnames = list(NULL, NULL, NULL, c("L", "U")))
    sc.means <- apply(sc.sp, c(1,2,4), mean, na.rm = T)
    sc.sds <- apply(sc.sp, c(1,2,4), sd, na.rm = T)
    snr.sc <- sc.means / sc.sds
    
    # per subpanel
    snr.sc.sp <- apply(sc.means, c(2:3), mean, na.rm = T) / apply(sc.sds, c(2:3), mean, na.rm = T)
    matplot(snr.sc.sp, type = "l", col = c("blue", "green3"), lty = 1)
    
    # plot SNR across whole panel
    {
        snr.sc.c <- array(snr.sc, dim = c(2048, 2), dimnames = list(NULL, c("L", "U")))
        plot(snr.sc.c[,"U"], type = "l", ylab = "", xlab = "", ylim = range(pretty(snr.sc.c)))
        lines(snr.sc.c[,"L"], col = "cyan3")
        abline(v = 128 * (1:15), col = "red", lty = 1)
    }
    # plot individual subpanels
    ul <- "U"; p = 4
    {
        plot(snr.sc[,p,ul], type = "l", ylim = range(pretty(snr.sc[,p,ul])), xlab = "", ylab = "",
             main = paste0(ul, p, "; gradient ", round(coef(line(snr.sc[,p,ul]))[2], 2), "; subpanel SNR ", round(snr.sc.sp[p,ul], 2)))
        # trends in SNR across panel?
        abline(line(snr.sc[,p,ul]), col = "red", lty = 3)
    }
}


# plot across whole detector
{
    snr.per.col.c <- array(snr.per.col, dim = c(2048, 2, 3), dimnames = dimnames(snr.per.col[1,,,]))
    plot(snr.per.col.c[,"U", "black"], type = "l", ylab = "", xlab = "", ylim = range(pretty(snr.per.col.c[,"U",])))
    lines(snr.per.col.c[,"U", "grey"], col = "skyblue")
    lines(snr.per.col.c[,"U", "white"], col = "red")
    abline(v = 128 * (1:15), col = "gold", lty = 2)
    
    plot(snr.per.col.c[,"L", "black"], type = "l", ylab = "", xlab = "", ylim = range(pretty(snr.per.col.c[,"L",])))
    lines(snr.per.col.c[,"L", "grey"], col = "skyblue")
    lines(snr.per.col.c[,"L", "white"], col = "red")
    abline(v = 128 * (1:15), col = "gold", lty = 2)
    
    # strange shape in 2nd & 3rd panels in from edge is probably to do with penumbra/spot edge.
}
# some exploratory plots
{
    ul <- "U"; p = 16
    plot(snr.per.col[,p,ul,"black"], type = "l", ylim = range(pretty(snr.per.col[,p,ul,])), xlab = "", ylab = "", main = paste0(ul, p))
    lines(snr.per.col[,p,ul,"grey"], col = "green3")
    lines(snr.per.col[,p,ul,"white"], col = "orange")
    
    # trends in SNR across panel?
    abline(line(snr.per.col[,p,ul,"black"]), col = "red", lty = 3)
    abline(line(snr.per.col[,p,ul,"grey"]), col = "red", lty = 3)
    abline(line(snr.per.col[,p,ul,"white"]), col = "red", lty = 3)
    print(round(apply(snr.per.col[,p,ul,], 2, function(cc) coef(line(cc))[2]), 2))
}

# SNR per subpanel
{
    snr.per.sp <- apply(sp.means, c(2:4), mean, na.rm = T) / apply(sp.sds, c(2:4), mean, na.rm = T)
    matplot(snr.per.sp[,"L",], type = "l", col = c("black", "green3", "orange"), lty = 1, ylab = "", ylim = range(pretty(snr.per.sp)))
    matplot(snr.per.sp[,"U",], type = "l", col = c("black", "green3", "orange"), lty = 2, add = T)
}

# what about local SNR?
{
    n <- 51; k <- matrix(1, nrow = n, ncol = n)
    local.m <- r2m(focal(m2r(pw.m[,,"black", "160430"]), k, mean, na.rm = T))
    local.sd <- r2m(focal(m2r(pw.m[,,"black", "160430"]), k, sd, na.rm = T))
    
    local.snr <- local.m / local.sd
    pixel.image(local.snr)
}
# too influenced by v local factors to be a particularly useful measure of anything, so no.