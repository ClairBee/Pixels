
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


