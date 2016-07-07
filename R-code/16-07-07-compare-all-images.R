
# PARAMETRIC MODELLING OF ALL IMAGES

library("IO.Pixels"); library("CB.Misc")

load.pixel.means()

gv.cols <- colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red", "magenta3"), space = "Lab")

####################################################################################################

# QUADRATIC TREND                                                                               ####

# fit to flat-field adjusted grey/white images (subtract black image first)
pwm.ff <- pw.m
pwm.ff[,,"white",] <- pwm.ff[,,"white",] - pwm.ff[,,"black",]
pwm.ff[,,"grey",] <- pwm.ff[,,"grey",] - pwm.ff[,,"black",]

# can't fit automatically (too big) - use function to output fitted values & coeffs separately
qt.lm <- function(im) {
    fit.lm <- lm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y),
                 setNames(melt(im), nm = c("x", "y", "gv")))
    list(fitted = array(fit.lm$fitted.values, dim = dim(im)),
         df = setNames(data.frame(t(c(coef(fit.lm),
                                    r2 = summary(fit.lm)$adj.r.squared,
                                    rmse = summary(fit.lm)$sigma))),
                       nm = c("Intercept", "x", "y", "x^2", "y^2", "xy", "r^2", "rmse")))
}

# fit all models
models <- apply(pwm.ff, 4, apply, 3, qt.lm)


# extract fitted values as array
fv <- array(dim = c(2048, 2048, 3, 12), dimnames = dimnames(pw.m))
fv[3:1998, 33:2028,,] <- abind(lapply(models, 
                   function(acq) abind(lapply(acq,
                                              function(model) model$fitted),
                                       along = 3)),
            along = 4)

saveRDS(fv, "./Image-plots/Quadratic-trends/qt.rds")
fv <- readRDS("./Image-plots/Quadratic-trends/qt.rds")


# extract coefficients etc as dataframe, for reference
df <- do.call("rbind", lapply(models, 
             function(model) do.call("rbind", lapply(model, "[[", "df"))))

# plot all quadratic trends
range(fv[3:1998,33:2028,"black",])   # 4700.49 6575.86
for (dt in dimnames(fv)[[4]]) {
    jpeg(paste0("./Image-plots/Quadratic-trends/quadratic-trends-ff-black-", dt, ".jpeg"))
            par(mar = c(2,2,1,1))
            image(1:2048, 1:2048, fv[,,"black",dt], breaks = c(46:66) * 100, col = gv.cols(20),
                  xlab = "", ylab = "")
            draw.panels(lty = 2)
        dev.off()
}

range(fv[3:1998,33:2028,"grey",])    # 9140.507 18022.526
for (dt in dimnames(fv)[[4]]) {
    jpeg(paste0("./Image-plots/Quadratic-trends/quadratic-trends-ff-grey-", dt, ".jpeg"))
    par(mar = c(2,2,1,1))
    image(1:2048, 1:2048, fv[,,"grey",dt], breaks = c((30:60) * 300) + 100, col = gv.cols(30),
          xlab = "", ylab = "")
    draw.panels(lty = 2)
    dev.off()
}

range(fv[3:1998,33:2028,"white",])   # 28271.12 46778.15
for (dt in dimnames(fv)[[4]]) {
    jpeg(paste0("./Image-plots/Quadratic-trends/quadratic-trends-ff-white-", dt, ".jpeg"))
    par(mar = c(2,2,1,1))
    image(1:2048, 1:2048, fv[,,"white",dt], breaks = c(28:47) * 1000, col = gv.cols(19),
          xlab = "", ylab = "")
    draw.panels(lty = 2)
    dev.off()
}

# create pdf image showing scale
pdf("./Image-plots/Quadratic-trends/qtscale-white.pdf", height = 1.5, width = 7); {
    par(mar = c(2, 1, 2, 1))
    image.scale(fv[,,"white", 1], zlim = c(28, 47) * 1000, 
                col = gv.cols(19), breaks = c(28:47) * 1000)
    title("White")
    dev.off()
}
pdf("./Image-plots/Quadratic-trends/qtscale-grey.pdf", height = 1.5, width = 7); {
    par(mar = c(2, 1, 2, 1))
    image.scale(fv[,,"grey", 1], zlim = c(30, 60) * 100, 
                col = gv.cols(30), breaks = c((30:60) * 300) + 100)
    title("Grey")
    dev.off()
}
pdf("./Image-plots/Quadratic-trends/qtscale-black.pdf", height = 1.5, width = 7); {
    par(mar = c(2, 1, 2, 1))
    image.scale(fv[,,"black", 1], zlim = c(46:66) * 100, 
                col = gv.cols(20), breaks = c(46:66) * 100)
    title("Black")
    dev.off()
}

####################################################################################################

# PLOT TRANSECT ALONG CENTRE LINE FOR MORE INTUITIVE COMPARISON                                 ####

# transect along central column
pdf("./Image-plots/Quadratic-trends/qt-col-hist-black.pdf"); {
    par(mar = c(4, 4, 1, 1))
    matplot(fv[1024,,"black",], type = "l", col = gv.cols(12), lty = 1, ylab = "Grey value", xlab = "Distance from bottom edge")
    legend("top", dimnames(fv)[[4]],  col = gv.cols(12), lty = 1, ncol = 2, bty = "n")
    dev.off()
}

# transect along central column
pdf("./Image-plots/Quadratic-trends/qt-col-hist-grey.pdf"); {
    par(mar = c(4, 4, 1, 1))
    matplot(fv[1024,,"grey",], type = "l", col = gv.cols(12), lty = 1, ylab = "Grey value", xlab = "Distance from bottom edge")
    legend("bottom", dimnames(fv)[[4]],  col = gv.cols(12), lty = 1, ncol = 2, bty = "n")
    dev.off()
}

pdf("./Image-plots/Quadratic-trends/qt-col-hist-white.pdf"); {
    par(mar = c(4, 4, 1, 1))
matplot(fv[1024,,"white",], type = "l", col = gv.cols(12), lty = 1, ylab = "Grey value", xlab = "Distance from bottom edge")
legend("bottom", dimnames(fv)[[4]],  col = gv.cols(12), lty = 1, ncol = 2, bty = "n")
dev.off()
}
