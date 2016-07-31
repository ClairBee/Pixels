
# parametric modelling of expected spot shape

library("IO.Pixels"); library("CB.Misc")
library(spatial)

fpath <- "./Notes/Parametric-spot/fig/"

acq <- readRDS("./02_Objects/images/pwm-loan.rds")
os.g <- acq[,,"grey"] - acq[,,"black"]
gv <- setNames(melt(os.g), nm = c("x", "y", "z"))

pw.m <- abind(sapply(c("131122", "140128", "MCT225", "160430", "loan"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)

# fit circular spot to each panel & report deviations

####################################################################################################

# CLEAN COPY OF CODE                                                                            ####

hh <- gaussian.spot.ls(os.g, c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                       x0.target = c(768, 1280), y0.target = c(768, 1280))

hh.fv <- gaussian.spot.mat(hh$par)
hh.res <- os.g - hh.fv
pixel.image(hh.res)

# lapply

# gaussian elliptical (rho = 0), no constraints on parameters
{
    # return vector of fitted values
    gaussian.ellipse.ls <- function(obs, param) {
            A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
            sig.x <- param["sig.x"]; sig.y <- param["sig.y"]
            
            est <- A * exp(- 0.5 * ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2))
            sum((est - obs$z)^2, na.rm = T)
    }
    
    gaussian.ellipse.mat <- function(param, obs) {
        
        est <- param["A"] * exp(-0.5 * ((((obs$x - param["x0"])/param["sig.x"])^2) + ((obs$y - param["y0"])/param["sig.y"])^2))
        
        array(est, dim = c(2048, 2048))
    }
    
    qq <- apply(pw.m, 4, 
                function(acq) {
                    gv <- setNames(melt(acq[,,"grey"] - acq[,,"black"]), nm = c("x", "y", "z"))
                    
                    # optimise with no constraints
                    optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                          gaussian.ellipse.ls, obs = gv)
                    })
    
    gaussian.ellipse.fv <- sapply(dimnames(pw.m)[[4]],
                                  function(dt) {
                                      gaussian.ellipse.mat(qq[,dt],
                                                           setNames(melt(pw.m[,,"grey", dt] - pw.m[,,"black", dt]),
                                                                    nm = c("x", "y", "z")))
                                  }, simplify = F)
    
    plot(pw.m[1024,,"grey", "loan"] - pw.m[1024,,"black", "loan"], type = "l")
    lines(gaussian.ellipse.fv$"loan"[1024,], col = "cyan3")
    
    # produce jpeg plots of all sets
    lapply(dimnames(pw.m)[[4]], function(dt) {
        jpeg(paste0(fpath, "surface-ge-unconst-", dt, ".jpg")); {
            par(mar = c(2,2,1,1))
            pixel.image(gaussian.ellipse.fv[[dt]][,])
            dev.off()
        }
        jpeg(paste0(fpath, "surface-res-ge-unconst-", dt, ".jpg")); {
            par(mar = c(2,2,1,1))
            pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt] - gaussian.ellipse.fv[[dt]][,])
            dev.off()
        }
        jpeg(paste0(fpath, "surface-trans-v-ge-unconst-", dt, ".jpg"), height = 240); {
            par(mar = c(2,2,1,1))
            plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[1024,], type = "l")
            lines(gaussian.ellipse.fv[[dt]][1024,], col = "orange", lwd = 2)
            dev.off()
        }
        jpeg(paste0(fpath, "surface-trans-h-ge-unconst-", dt, ".jpg"), height = 240); {
            par(mar = c(2,2,1,1))
            plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[,1024], type = "l", lwd = 2)
            lines(gaussian.ellipse.fv[[dt]][,1024], col = "orange", lwd = 2)
            dev.off()
        }
    })
}

# Gaussian elliptical (rho = 0), midpoints must be in central region
# tighten constraint: must be in central quarter (768:1280)
{
    qq.const <- apply(pw.m, 4, 
                function(acq) {
                    gv <- setNames(melt(acq[,,"grey"] - acq[,,"black"]), nm = c("x", "y", "z"))
                    
                    # optimise with constraints
                    # in particular, spot centre must be in central third of panel
                    optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                          gaussian.ellipse.ls, obs = gv, method = "L-BFGS-B",
                          lower = c(-Inf, 768, 768, 0, 0), upper = c(Inf, 1280, 1280, Inf, Inf))
                })
    
    gaussian.ellipse.const.fv <- sapply(dimnames(pw.m)[[4]],
                                  function(dt) {
                                      gaussian.ellipse.mat(qq.const[[dt]]$par,
                                                           setNames(melt(pw.m[,,"grey", dt] - pw.m[,,"black", dt]),
                                                                    nm = c("x", "y", "z")))
                                  }, simplify = F)
    
    # produce jpeg plots of all sets
    lapply(dimnames(pw.m)[[4]], function(dt) {
        jpeg(paste0(fpath, "surface-ge-const-", dt, ".jpg")); {
            par(mar = c(2,2,1,1))
            pixel.image(gaussian.ellipse.const.fv[[dt]][,])
            dev.off()
        }
        jpeg(paste0(fpath, "surface-res-ge-const-", dt, ".jpg")); {
            par(mar = c(2,2,1,1))
            pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt] - gaussian.ellipse.const.fv[[dt]][,])
            dev.off()
        }
        jpeg(paste0(fpath, "surface-trans-v-ge-const-", dt, ".jpg"), height = 240); {
            par(mar = c(2,2,1,1))
            plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[1024,], type = "l")
            lines(gaussian.ellipse.const.fv[[dt]][1024,], col = "orange", lwd = 2)
            dev.off()
        }
        jpeg(paste0(fpath, "surface-trans-h-ge-const-", dt, ".jpg"), height = 240); {
            par(mar = c(2,2,1,1))
            plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[,1024], type = "l", lwd = 2)
            lines(gaussian.ellipse.const.fv[[dt]][,1024], col = "orange", lwd = 2)
            dev.off()
        }
    })
}

####################################################################################################

# SURFACES FITTED TO ALL GREY IMAGES                                                            ####

# helper function to return surface fitted by least squares
surface.trend <- function(im, order = 2) {
    gv <- setNames(melt(im), nm = c("x", "y", "z"))
    s.ls <- surf.ls(order, gv[!is.na(gv$z),])
    trmat(s.ls, 1, 2048, 1, 2048, 2047)$z
}

st.2 <- array(apply(pw.m, 4, function(acq) surface.trend(acq[,,"grey"] - acq[,,"black"], order = 2)),
              dim = c(2048, 2048, dim(pw.m)[[4]]), dimnames = list(NULL, NULL, dimnames(pw.m)[[4]]))


lapply(dimnames(pw.m)[[4]], function(dt) {
    jpeg(paste0(fpath, "offset-image-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey",dt] - pw.m[,,"black",dt])
        dev.off()
    }
    jpeg(paste0(fpath, "surface-o2-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(st.2[,,dt])
        dev.off()
    }
    jpeg(paste0(fpath, "surface-res-o2-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt] - st.2[,,dt])
        dev.off()
    }
    jpeg(paste0(fpath, "surface-trans-v-o2-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[1024,], type = "l")
        lines(st.2[1024,,dt], col = "orange", lwd = 2)
        dev.off()
    }
    jpeg(paste0(fpath, "surface-trans-h-o2-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[,1024], type = "l", lwd = 2)
        lines(st.2[,1024,dt], col = "orange", lwd = 2)
        dev.off()
    }
})

st.2.res <- pw.m[,,"grey",] - pw.m[,,"black",] - st.2

pixel.image(st.2.res[,,"131122"])
hist(st.2.res[,,"131122"], breaks = "fd")

px.131122 <- which
plot(which(st.2.res[,,"131122"] > 1000, arr.ind = T), pch = 15, cex = 0.5, col = "red", xlim = c(0,2048), ylim = c(0,2048))
points(which(st.2.res[,,"131122"] < -1000, arr.ind = T), pch = 15, cex = 0.5, col = "blue")

sd(st.2.res[,,"131122"], na.rm = T)

plot(which(st.2.res[,,"160430"] > 1000, arr.ind = T), main = "160430", pch = 15, cex = 0.5, col = "red", xlim = c(0,2048), ylim = c(0,2048))
points(which(st.2.res[,,"160430"] < -1000, arr.ind = T), pch = 15, cex = 0.5, col = "blue")

# can now identify high-residual areas & identify dense defect regions

plot(which(st.2.res[,,"loan"] > 1000, arr.ind = T), main = "160430", pch = 15, cex = 0.5, col = "red", xlim = c(0,2048), ylim = c(0,2048))
points(which(st.2.res[,,"loan"] < -1000, arr.ind = T), pch = 15, cex = 0.5, col = "blue")

plot(which(st.2.res[,,"160430"] > 1000, arr.ind = T), main = "160430", pch = 15, cex = 0.5, col = "red", xlim = c(900,1000), ylim = c(900,1000))
points(which(st.2.res[,,"160430"] < -1000, arr.ind = T), pch = 15, cex = 0.5, col = "blue")

pw.m[,,"grey", "160430"][which(st.2.res[,,"160430"] < -1000, arr.ind = T)]

#----------------------------------------------------------------------------------------------------
# bivariate Gaussian, no interaction (~6m to run for 5 images)
xy.gaussian.trend <- function(im) {

}

gt.ind <- apply(pw.m, 4, function(acq) xy.gaussian.trend(acq[,,"grey"] - acq[,,"black"]))
gt.ind.mat <- function(param, obs) {
    
    A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]
    
    est <- A * exp(- 0.5 * ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2))   
    
    array(est, dim = c(2048, 2048))
}
lapply(dimnames(pw.m)[[4]], function(dt) {
    jpeg(paste0(fpath, "surface-gt-ind-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(gt.ind[,,dt])
        dev.off()
    }
    jpeg(paste0(fpath, "surface-res-gt-ind-", dt, ".jpg")); {
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt] - gt.ind[,,dt])
        dev.off()
    }
    jpeg(paste0(fpath, "surface-trans-v-gt-ind-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[1024,], type = "l")
        lines(gt.ind[1024,,dt], col = "orange", lwd = 2)
        dev.off()
    }
    jpeg(paste0(fpath, "surface-trans-h-gt-ind-", dt, ".jpg"), height = 240); {
        par(mar = c(2,2,1,1))
        plot((pw.m[,,"grey", dt] - pw.m[,,"black", dt])[,1024], type = "l", lwd = 2)
        lines(gt.ind[,1024,dt], col = "orange", lwd = 2)
        dev.off()
    }
})

zz <- lapply(models, gt.ind.mat, obs = gv)

lines(zz[[5]][1024,], col = "cyan3")

# split out workhorse function from wrappers

####################################################################################################

# QUADRATIC TREND WITH CENTRE PARAMETER                                                         ####

qt.const <- function(par, x, y, z) {
    cc <- par["cc"]; a1 <- par["a1"]; a2 <- par["a2"]; b1 <- par["b1"]; b2 <- par["b2"];
    x0 <- par["x0"]; y0 <- par["y0"]
    est <- cc + a1 * (x-x0) + a2 * (x - x0)^2 + b1 * (y-y0) + b2 * (y-y0)^2
    
    sum((est - z)^2, na.rm = T)
}

gv.tst <- setNames(melt(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"]), nm = c("x", "y", "z"))
zz <- optim(c(cc = 15000, x0 = 1024.5, y0 = 1024.5, a1 = 1, a2 = 1, b1 = 1, b2 = 1),
            qt.const, x = gv.tst$x, y = gv.tst$y, z = gv.tst$z, method = "L-BFGS-B",
            lower = c(-Inf, 683, 683, -Inf, -Inf, -Inf, -Inf), upper = c(Inf, 1365, 1365, Inf, Inf, Inf, Inf))

####################################################################################################

# ELLIPTICAL SPOT                                                                               ####

# mystery hybrid function...
# elliptical, no eccentricity (independent x & y), exponential curve; manual function (ls, optim)
{
    ellipse.exp.ls <- function(obs, param) {
        A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
        sig.x <- param["sig.x"]; sig.y <- param["sig.y"]
        
        est <- A * exp(- 0.5 * ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2))
        sum((est - obs$z)^2, na.rm = T)
    }
    
    zz <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700),
                ellipse.exp.ls, obs = gv, )
    gv$fv.g2d.ls <- zz$par["A"] / (2 * pi * zz$par["sig.x"] * zz$par["sig.y"]) * exp(-((((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + ((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2))
    fv.g2d.ls <- array(gv$fv.g2d.ls, dim = c(2048, 2048))
    res.g2d.ls <- os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048))
    
    # plot fitted model
    {
        pixel.image(array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Fitted values")
        pixel.image(os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Residuals")
        
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1:1024], ylim = c(0,1000), breaks = "fd", main = "Residuals after circular spot")
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1025:2048], breaks = "fd", 
             add = T, border = adjustcolor("red", alpha = 0.3))
        legend("topright", col = c("black", "red"), pch = 15, legend = c("Lower", "Upper"), bty = "n")
        
        image(1:2048, 1:2048, os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), asp = T,
              breaks = c(-20000, -500, 500, 20000), col = c("blue", NA, "red"), xlim = c(1900,2048), ylim = c(0,100))
        
        plot(os.g[1024,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1024,], col = "blue")
        
        plot(os.g[1025,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1025,], col = "green3")
        
        plot(os.g[, 1024], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[, 1024], col = "green3")
        
        plot(os.g[,1025], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[,1025], col = "green3")
        
        plot(os.g[512,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[512,], col = "green3")
    }
    
    # now fit subpanels to model
    {
        panel.os <- lm(fv.g2d.ls ~ tt, data = gv[!is.na(gv$z),])
        pp <- array(dim = c(2048, 2048))
        pp[25:2024, 25:2024] <- panel.os$fitted.values - coef(panel.os)[1]
        pixel.image(pp)
        pixel.image(res.g2d.ls - pp)
        pp <- predict(panel.os, newdata = gv$tt)
    }
    
    # function to apply directly to image
    fit.ellipse <- function(im, A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700) {
        
        surf.dat <- setNames(melt(im), nm = c("x", "y", "z"))
        
        optim(c(A = A, x0 = x0, y0 = y0, sig.x = sig.x, sig.y = sig.y), ellipse.exp.ls, obs = surf.dat)
    }
    
    zz <- fit.ellipse(im = pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])
    fv.g2d.ls <- array(zz$par["A"] * exp(-((((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + ((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2)),
                       dim = c(2048, 2048))
    
    pixel.image(fv.g2d.ls)
    pixel.image(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"] - fv.g2d.ls)
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1024,], type = "l")
    lines(fv.g2d.ls[1024,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1025,], type = "l")
    lines(fv.g2d.ls[1025,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1024], type = "l")
    lines(fv.g2d.ls[,1024], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1025], type = "l")
    lines(fv.g2d.ls[,1025], col = "blue")
}

# elliptical, scaled quadratic curve
{
    ellipse.ls <- function(obs, param) {
        A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
        sig.x <- param["sig.x"]; sig.y <- param["sig.y"]
        
        est <- A + ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2)
        sum((est - obs$z)^2, na.rm = T)
    }
    
    zz <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700),
                ellipse.ls, obs = gv)
    fv.ellipse <- array(zz$par["A"] + (((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + (((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2),
                    dim = c(2048, 2048))
    pixel.image(fv.ellipse)
        res.g2d.ls <- os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048))
    
    # plot fitted model
    {
        pixel.image(array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Fitted values")
        pixel.image(os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Residuals")
        
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1:1024], ylim = c(0,1000), breaks = "fd", main = "Residuals after circular spot")
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1025:2048], breaks = "fd", 
             add = T, border = adjustcolor("red", alpha = 0.3))
        legend("topright", col = c("black", "red"), pch = 15, legend = c("Lower", "Upper"), bty = "n")
        
        image(1:2048, 1:2048, os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), asp = T,
              breaks = c(-20000, -500, 500, 20000), col = c("blue", NA, "red"), xlim = c(1900,2048), ylim = c(0,100))
        
        plot(os.g[1024,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1024,], col = "blue")
        
        plot(os.g[1025,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1025,], col = "green3")
        
        plot(os.g[, 1024], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[, 1024], col = "green3")
        
        plot(os.g[,1025], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[,1025], col = "green3")
        
        plot(os.g[512,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[512,], col = "green3")
    }
    
    # now fit subpanels to model
    {
        panel.os <- lm(fv.g2d.ls ~ tt, data = gv[!is.na(gv$z),])
        pp <- array(dim = c(2048, 2048))
        pp[25:2024, 25:2024] <- panel.os$fitted.values - coef(panel.os)[1]
        pixel.image(pp)
        pixel.image(res.g2d.ls - pp)
        pp <- predict(panel.os, newdata = gv$tt)
    }
    
    # function to apply directly to image
    fit.ellipse <- function(im, A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700) {
        
        surf.dat <- setNames(melt(im), nm = c("x", "y", "z"))
        
        optim(c(A = A, x0 = x0, y0 = y0, sig.x = sig.x, sig.y = sig.y), ellipse.exp.ls, obs = surf.dat)
    }
    
    zz <- fit.ellipse(im = pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])
    fv.g2d.ls <- array(zz$par["A"] * exp(-((((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + ((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2)),
                       dim = c(2048, 2048))
    
    pixel.image(fv.g2d.ls)
    pixel.image(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"] - fv.g2d.ls)
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1024,], type = "l")
    lines(fv.g2d.ls[1024,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1025,], type = "l")
    lines(fv.g2d.ls[1025,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1024], type = "l")
    lines(fv.g2d.ls[,1024], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1025], type = "l")
    lines(fv.g2d.ls[,1025], col = "blue")
}

####################################################################################################

# 2D GAUSSIAN SPOT                                                                              ####
# allows for eccentricity through change in covariance matrix
bvn <- function(x, y, x0, y0, sig.x, sig.y, rho) {
    const <- (2 * pi * sig.x * sig.y * sqrt(1 - rho^2))^-1
    rr <- (2 * (1 - rho^2))^-1
    
    z.x <- (x - x0) / sig.x
    z.y <- (y - y0) / sig.y
    
    const * exp(- rr * (z.x^2 + z.y^2 - (2 * rho * z.x * z.y)))  
}

gauss.2d.ls <- function(param, x, y, z) {
    A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- A * bvn(x, y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - z)^2, na.rm = T)
}


# try MLE using built-in formula
gauss.2d.ml <- function(obs, param) {
    A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- A + bvn(x = obs$x, y = obs$y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - obs$z)^2, na.rm = T)
    
    LL <- function(mu, sigma) {
        +     R = dnorm(x, mu, sigma)
        +     #
            +     -sum(log(R))
        + }
}


    
g2d <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5,
               sig.x = sd(gv$z, na.rm = T), sig.y = sd(gv$z, na.rm = T), rho = 0),
             gauss.2d.ls, x = gv$x, y = gv$y, z = gv$z, method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 2048, 2048, Inf, Inf, 1))
g2d$par
}

hh <- bvn(gv$x, gv$y, g2d$par["x0"], g2d$par["y0"], g2d$par["sig.x"], g2d$par["sig.y"], g2d$par["rho"])
pixel.image(array(hh, dim = c(2048, 2048)))

####################################################################################################

# MWE FOR STACK EXCHANGE                                                                        ####

# bivariate normal density function
bvn <- function(x, y, x0, y0, sig.x, sig.y, rho) {
    const <- 1 / (2 * pi * sig.x * sig.y * sqrt(1 - rho^2))
    rr <- 1 / (2 * (1 - rho^2))^-1
    
    z.x <- (x - x0) / sig.x
    z.y <- (y - y0) / sig.y
    
    const * exp(- rr * (z.x^2 + z.y^2 - (2 * rho * z.x * z.y)))  
}

# create some fake data
x <- seq(0, 500, 25); y <- seq(0, 500, 25)
z <- outer(x, y, bvn, x0 = 250, y0 = 250, sig.x = 100, sig.y = 100, rho = 0)

fdat <- setNames(melt(z), nm = c("x", "y", "z"))

gauss.2d.ls <- function(param, x, y, z) {
    x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- bvn(x, y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - z)^2, na.rm = T)
}

g2d <- optim(c(x0 = 240, y0 = 240, sig.x = 80, sig.y = 80, rho = 0),
             gauss.2d.ls, x = fdat$x, y = fdat$y, z = fdat$z, method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 250, 250, Inf, Inf, 1))
g2d$par

#           x0           y0        sig.x        sig.y          rho 
# 2.400000e+02 2.400000e+02 8.000000e+01 8.000000e+01 1.985233e-20


####################################################################################################

# QUADRATIC TREND SURFACE BY LS                                                                 ####

library(spatial)
s.ls <- surf.ls(2, gv[!is.na(gv$z), c("x", "y", "z")])
trsurf <- trmat(s.ls, 1, 2048, 1, 2048, 2047)
pixel.image(trsurf$z)

points(which(trsurf$z == max(trsurf$z), arr.ind = T), pch = 3)      # 982, 1105

plot(os.g[982,], type = "l")
lines(trsurf$z[982,], col = "blue")

plot(os.g[,1105], type = "l")
lines(trsurf$z[,1105], col = "blue")



# does my bivariate normal function fit successfully here?
{
    g2d.ls <- optim(c(A = 15000, x0 = 25, y0 = 25,
                      sig.x = sd(trsurf$z), sig.y = sd(trsurf$z), rho = 0),
                    gauss.2d.ls, obs = setNames(melt(trsurf$z), nm = c("x", "y", "z")), method = "L-BFGS-B",
                    lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 25, 25, Inf, Inf, 1))
    g2d.ls$par
}

# nope.

hh <- bvn(trsurf$x, trsurf$y, g2d.ls$par["x0"], g2d.ls$par["y0"], g2d.ls$par["sig.x"], g2d.ls$par["sig.y"], g2d.ls$par["rho"])
pixel.image(array(hh, dim = c(2048, 2048)))

# residual plot
pixel.image(os.g - trsurf$z)
o.plot(os.g[1024,])
lines(trsurf$z[1024,], col = "green3")

####################################################################################################

# 2D EXAMPLES                                                                                   ####

# Rosenbrock surface
{
    f <- function(x1,y1) (1-x1)^2 + 100*(y1 - x1^2)^2
    x <- seq(-2,2,by=.15)
    y <- seq(-1,3,by=.15)
    z <- outer(x,y,f)
    persp(x,y,z,phi=45,theta=-45,col="yellow",shade=.00000001,ticktype="detailed")
    
    # When using optim for multidimensional optimization, the input in your function definition must be a single vector
    f <- function(x) (1-x[1])^2 + 100*(x[2]-x[1]^2)^2
    
    # starting values must be a vector now
    optim( c(0,0), f )$par      # [1] 0.9999564 0.9999085
    
}

# Himmelblau's function
{
    f <- function(x1,y1) (x1^2 + y1 - 11)^2 + (x1 + y1^2 - 7)^2
    x <- seq(-4.5,4.5,by=.2)
    y <- seq(-4.5,4.5,by=.2)
    z <- outer(x,y,f)
    persp(x,y,z,phi=-45,theta=45,col="yellow",shade=.65 ,ticktype="detailed")
    
    f <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
    optim(c(-4,-4), f)$par
    optim(c(2,-2), f)$par
    optim(c(2,2), f)$par
    optim(c(-4,4),f)$par
}

# solving a maximisation problem
{
    f <- function(x,y) 2*x*(y**2)+2*(x**2)*y+x*y
    
    # set x & y values, get z
    x <- seq(-0.5,0.5, len=200)
    y <- seq(-0.5,0.5, len=200)
    z <- outer(x,y,f)
    image(x,y,z)
    
    # partial derivative functions
    fx <- function(x,y,h=0.001) (f(x+h,y)-f(x,y))/h
    fy <- function(x,y,h=0.001) (f(x,y+h)-f(x,y))/h
    
    zfx <- outer(x,y,fx)
    zfy <- outer(x,y,fy)
    
    contour(x,y,zfx,level=0)
    contour(x,y,zfy,level=0, add=T, col="red")
    
    # create new function for optim
    fbb<-function(x) f(x[1],x[2])
    optim(c(0.5,0.5),fbb,control=list(fnscale=-1))
    
    fxb <- function(x) fx(x[1],x[2])
    fyb <- function(x) fy(x[1],x[2])
    
    sumssq <- function(x) fxb(x)**2+fyb(x)**2       # sum of squares
    
    optim(c(0.1,0.1),sumssq)
    optim(c(-0.4,0.1),sumssq)
    
    optim(c(0.4,0.4),fbb,control=list(fnscale=-1))
}
####################################################################################################

# UPPER VS LOWER PANEL OFFSET                                                                   ####

# create factor to fit offset
gv$upper <- gv$y > 1024.5

ul.os <- lm(z ~ upper, data = gv[!is.na(gv$z),])
coef(ul.os)["upperTRUE"]

gv$p.adj <- gv$z - (gv$upper * coef(ul.os)["upperTRUE"])

pixel.image(array(gv$p.adj, dim = c(2048, 2048)))
o.plot(array(gv$z, dim = c(2048, 2048))[1024,])
o.plot(array(gv$z, dim = c(2048, 2048))[512,])

# or, maybe manually will be less mad - aim to smooth centre line by offsetting upper panel
# use median to remove effect of dead lines on either side
d <- apply(os.g[,1025 + c(0:19)], 1, median) - apply(os.g[,1024 -  + c(0:19)], 1, median)
median(d, na.rm = T)

gv$p.adj <- gv$z - (gv$upper * median(d, na.rm = T))

pixel.image(array(gv$p.adj, dim = c(2048, 2048)))
o.plot(array(gv$p.adj, dim = c(2048, 2048))[1024,])
o.plot(array(gv$p.adj, dim = c(2048, 2048))[512,])

gv2 <- setNames(gv[,c("x", "y", "p.adj")], nm = c("x", "y", "z"))
zz.adj <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700),
            g2d.ls, obs = gv2)
zz.adj$par
zz$par
# plot fitted model
{
    gv$fv.g2d.ls <- zz.adj$par["A"] * exp(-((((gv$x - zz.adj$par["x0"]) / zz.adj$par["sig.x"])^2) + ((gv$y - zz.adj$par["y0"]) / zz.adj$par["sig.y"])^2))
    pixel.image(array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Fitted values")
    pixel.image(os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Residuals")
    
    hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1:1024], breaks = "fd", ylim = c(0,30))
    hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1025:2048], breaks = "fd", 
         add = T, border = adjustcolor("red", alpha = 0.3))
    
    plot(os.g[1024,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1024,], col = "green3")
    
    plot(os.g[1025,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1025,], col = "green3")
    
    plot(os.g[, 1024], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[, 1024], col = "green3")
    
    plot(os.g[,1025], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[,1025], col = "green3")
    
    plot(os.g[512,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[512,], col = "green3")
}

####################################################################################################

# PER-PANEL OFFSET                                                                              ####

# linear regression (simpler approach)
{
    d <- apply(os.g[,1025 + c(0:19)], 1, mean) - apply(os.g[,1024 -  + c(0:19)], 1, mean)
    o.plot(d)
    
    pixel.image(os.g[,1025:2048])
    
    sp <- array(dim = c(128, 16, 1024, 2), dimnames = list(NULL, NULL, NULL, c("L", "U")))
    
    for (ul in dimnames()) {
        for (p in 1:16) {
            sp[,p,,ul] <- paste0(ul, formatC(p, width = 2, flag = "0"))
        }
    }
    tt <- array(sp, dim = c(2048, 2048))
    
    gv$tt <- c(tt)
    
    panel.os <- lm(z ~ tt, data = gv[!is.na(gv$z),])
    
    pp <- array(predict(panel.os), dim = c(2000,2000))
    os.g.pa <- os.g
    os.g.pa[25:2024, 25:2024] <- os.g[25:2024, 25:2024] - pp
    pixel.image(os.g.pa)
}

####################################################################################################

# NLM / NLS / OPTIM: BUILD FUNCTION                                                             ####

ff <- function(param, obs) {
    
    est <- (obs$x - param["x0"])^2 + (obs$y - param["y0"])^2
    sum((est - obs$z)^2, na.rm = T)
}

ff.optim <- optim(c(x0 = 1000, y0 = 2000), ff, obs = gv)
{
    ff.nls <- nls(z ~ (x - x0)^2 + (y - y0)^2, data = gv, start = c(x0 = 1000, y0 = 2000))
    # repeated failure to converge
    
    ff.nlm <- nlm(ff, c(x0 = 1000, y0 = 2000), obs = gv)
    # failed to converge: gradient too close to 0.
}

#---------------------------------------------------------------------------------
# back to optim. Let's try to build up the complexity...
# 2d gaussian, no interaction
ff <- function(param, obs) {
    
    est <- param["A"] * exp(-0.5 * ((((obs$x - param["x0"])/param["sig.x"])^2) + ((obs$y - param["y0"])/param["sig.y"])^2))
    sum((est - obs$z)^2, na.rm = T)
}
ff.mat <- function(param, obs) {
    
    est <- param["A"] * exp(-0.5 * ((((obs$x - param["x0"])/param["sig.x"])^2) + ((obs$y - param["y0"])/param["sig.y"])^2))
    
    array(est, dim = c(2048, 2048))
}

ff.optim <- optim(c(A = 15000, x0 = 1000, y0 = 2000, sig.x = 500, sig.y = 500), ff, obs = gv)
ff.optim$par
pixel.image(ff.mat(ff.optim$par, gv))

plot(os.g[1024,], type = "l")
lines(ff.mat(ff.optim$par, gv)[1024,], col = "orange")

plot(os.g[,1024], type = "l")
lines(ff.mat(ff.optim$par, gv)[,1024], col = "orange")

#---------------------------------------------------------------------------------
# circular spot, centre can move
sp.lm <- function(param, obs) {
    obs$dd <- sqrt((obs$x - param["x0"])^2 + (obs$y - param["y0"])^2)
    
    lm(z ~ poly(dd, 2), obs)
}

sp <- sp.lm(c(x0 = 1024, y0 = 1024), obs = gv)
pixel.image(array(sp$fitted.values, dim = c(2000,2000)))

ff <- function(param, obs) {
    dd <- sqrt((obs$x - param["x0"])^2 + (obs$y - param["y0"])^2)
    
    est <- param["a"] + param["b"] * dd + param["c"] * dd^2
    
    sum((est - obs$z)^2, na.rm = T)
}

ff.optim <- optim(c(a = 15000, b = -2000000, c = -40000, x0 = 1024, y0 = 1024), ff, obs = gv)
ff.optim$par

ff.mat <- function(param, obs) {
    
    dd <- sqrt((obs$x - param["x0"])^2 + (obs$y - param["y0"])^2)
    est <- param["a"] + param["b"] * dd + param["c"] * dd^2    
    
    array(est, dim = c(2048, 2048))
}
pixel.image(ff.mat(ff.optim$par, gv))
pixel.image(os.g)

rr <- nlm(ff, c(cc = 15000, x0 = 1024, y0 = 1024, a1 = 50, a2 = 50, b1 = 50, b2 = 50), obs = gv)

qq <- spot.lm(os.g)

####################################################################################################

# BIVARIATE GAUSSIAN MODEL FITTED TO BLACK IMAGES                                               ####

pixel.image(pw.m[,,"black", "160705"])

# suppose no particular reason to assume circularity, so stick with quadratic spot (more flexible)
library(spatial)

gv <- setNames(melt(pw.m[,,"black", "160705"]), nm = c("x", "y", "z"))
s.ls <- surf.ls(2, gv[!is.na(gv$z), c("x", "y", "z")])
trsurf <- trmat(s.ls, 1, 2048, 1, 2048, 2047)
pixel.image(trsurf$z)

qs.res <- pw.m[,,"black", "160705"] - trsurf$z
pixel.image(qs.res)
rect(60,60,1988,1988, lty = 3)

matplot(pw.m[1:128, c(c(0:15) * 128 + 1, 2048), "black", "160705"], type = "l", ylim = c(5000,8000))

gv.res <- setNames(melt(qs.res), nm = c("x", "y", "z"))
click()

sp <- array(qs.res, dim = c(128, 16, 1024, 2))
pixel.image(sp[,1,,2])

####################################################################################################

# MEAN IMAGES                                                                                   ####

wmg <- apply(pw.m[,,,7:19], 3, apply, c(1:2), mean, na.rm = T)
wmg <- array(wmg, dim = c(2048, 2048, 3), dimnames = list(NULL, NULL, c("black", "grey", "white")))

pixel.image(wmg[,,"white"])
