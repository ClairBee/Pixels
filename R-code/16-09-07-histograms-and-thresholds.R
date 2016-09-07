
library("IO.Pixels"); library("CB.Misc")

# LINE DETECTION ALGORITHM UNSATISFACTORY

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")

####################################################################################################

# JOHNSON DISTRIBUTIONS OF ALL BLACK & GREY IMAGES                                              ####

fpath <- "./Image-plots/histograms/"

normal.Johnson.QQ <- function (data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99), title = "Normal Q-Q plot", ...) {
    
    data <- data[!is.na(data)]
    m <- mean(data)
    s <- mad(data)
    jf <- JohnsonFit(data, moment = "quant")
    
    plot(qnorm(quantiles, m, s), quantile(data, quantiles), 
             pch = 20, asp = T, ylab = "Observed quantile", xlab = "Fitted quantile", 
             main = title, col = adjustcolor("seagreen", alpha = 0.3), ...)
    abline(h = quantile(data, grid.quantiles), col = "lightseagreen", lty = 2)
    abline(v = qnorm(grid.quantiles, m, s), col = "lightseagreen", lty = 2)
    
    points(qJohnson(quantiles, jf), quantile(data, quantiles), pch = 20)
    abline(v = qJohnson(grid.quantiles, jf), col = "skyblue", lty = 2)
    abline(v = qnorm(grid.quantiles, m, s), col = "skyblue", lty = 2)
    
    abline(0, 1, col = "red", lty = 2)

    legend("topleft", cex = 0.7, title = "Gridlines: ", legend = paste0("Q", grid.quantiles * 100), lty = 2, col = "skyblue", box.col = "white")
    legend("bottomright", cex = 0.7, title = "Points", col = c("seagreen", "black"), pch = 20, box.col = "white",
           legend = c("Normal (MAD)", "Johnson"))
}

# histograms with fit, Q-Q plot & identified defects: black images
invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "hist-black-", dt, ".pdf"))
                     par(mfrow = c(2,2), mar = c(3,3,3,1))
                     
                     JF <- JohnsonFit(pw.m[,,"black", dt][!is.na(pw.m[,,"black", dt])])
                     
                     hist(pw.m[,,"black", dt], breaks = "fd", prob = T, xlim = c(0,10000), main = paste0(dt, " - black"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), sd(pw.m[,,"black", dt], na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), mad(pw.m[,,"black", dt], na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = median(pw.m[,,"black", dt], na.rm = T) + c(-6,-5,5,6) * mad(pw.m[,,"black", dt], na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topright", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     hist(pw.m[,,"black", dt], breaks = "fd", prob = T, ylim = c(0,0.0005), xlim = c(0,10000), main = paste0(dt, " - black (cropped)"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), sd(pw.m[,,"black", dt], na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), mad(pw.m[,,"black", dt], na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = median(pw.m[,,"black", dt], na.rm = T) + c(-6,-5,5,6) * mad(pw.m[,,"black", dt], na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     legend("topright", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     normal.Johnson.QQ(pw.m[,,"black", dt], title = paste0("Q-Q plot - ", dt, " black"))

                     med <- median(pw.m[,,"black", dt], na.rm = T)
                     sp <- mad(pw.m[,,"black", dt], na.rm = T)
                     
                     pixel.plot(which(pw.m[,,"black", dt] > med + 6 * sp | pw.m[,,"black", dt] < med - 6 * sp, arr.ind = T), col = "gold", cex = 0.3, main = paste0("Defects - ", dt))
                     points(which(pw.m[,,"black", dt] > qJohnson(pnorm(5, 0, 1), JF) | pw.m[,,"black", dt] < qJohnson(pnorm(-5, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3, col = "cyan3")
                     points(which(pw.m[,,"black", dt] > qJohnson(pnorm(6, 0, 1), JF) | pw.m[,,"black", dt] < qJohnson(pnorm(-6, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3)
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))


# histograms with fit, Q-Q plot & identified defects: grey images
invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "hist-grey-", dt, ".pdf"))
                     par(mfrow = c(2,2), mar = c(3,3,3,1))
                     
                     f.im <- im <- pw.m[,,"grey", dt]
                     
                     zz <- density(im, n = 65536, na.rm = T)
                     mu <- zz$x[which.max(zz$y)]
                     
                     # filter image to remove dark pixels
                     f.im <- f.im[which(f.im > 10000)]
                     JF <- JohnsonFit(f.im)
                     
                     hist(im, breaks = "fd", prob = T, xlim = c(0, 26000), main = paste0(dt, " - grey"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), sd(f.im, na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), mad(f.im, na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = mu + c(-6,-5,5,6) * mad(im, na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topleft", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     hist(im, breaks = "fd", prob = T, ylim = c(0,0.0005), xlim = c(0, 26000), main = paste0(dt, " - grey (cropped)"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), sd(f.im, na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), mad(f.im, na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = mu + c(-6,-5,5,6) * mad(im, na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topleft", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     normal.Johnson.QQ(f.im, title = paste0("Q-Q plot - ", dt, " grey"))
                     
                     med <- median(f.im, na.rm = T)
                     sp <- mad(f.im, na.rm = T)
                     
                     pixel.plot(which(im > mu + 6 * sp | im < mu - 6 * sp, arr.ind = T), col = "gold", cex = 0.3, main = paste0("Defects - ", dt))
                     points(which(im > qJohnson(pnorm(5, 0, 1), JF) | im < qJohnson(pnorm(-5, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3, col = "cyan3")
                     points(which(im > qJohnson(pnorm(6, 0, 1), JF) | im < qJohnson(pnorm(-6, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3)
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))

####################################################################################################

# HALF-NORMAL DISTRIBUTION FOR HISTOGRAMS                                                      ####

library(VGAM)

dat <- pw.m[,,"black", "160430"]
hist(dat, breaks = c(0:65535), xlim = c(4000,6000), prob = T)
zz <- density(dat, n = 65536, na.rm = T)

mu <- zz$x[which.max(zz$y)]
abline(v = mu, col = "red")

# break at point of maximum density
dat.l <- abs(dat[which(dat <= mu)] - mu)
dat.u <- dat[which(dat > mu)] - mu

sig.l <- sqrt(mean(dat.l^2))
sig.u <- sqrt(mean(dat.u^2))

hist(dat.l, breaks = "fd", xlim = c(0,5000), prob = T)
lines(0:5000, dhalfnorm(0:5000, sd2theta(sig.l)), col = "red", lwd = 2)

hist(dat.u, breaks = "fd", xlim = c(0,5000), prob = T)
lines(0:5000, dhalfnorm(0:5000, sd2theta(sig.u)), col = "red", lwd = 2)

# images with strong upward drift are poorly fitted by half-normal. Hey ho.

####################################################################################################

# INDEPENDENT SD FOR DEFECT IDENTIFICATION                                                      ####

# treat data as two approximately half-normal distributions
# cut at modal density
# find MAD of upper & lower parts

# is MAD of half-normal same as that of mirrored normal?

dat <- pw.m[,,"black", "160430"]

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "") {
    
    zz <- density(dat, n = 65536, na.rm = T)
    mu <- zz$x[which.max(zz$y)]

    # break at point of maximum density
    dat.l <- dat[which(dat <= mu)] - mu
    dat.u <- dat[which(dat > mu)] - mu
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "")
    abline(v = mu, col = "red", lwd = 1)
    abline(v = mu - mad(c(dat.l, -dat.l)) * c(5,6), lty = c(2,3), col = "red")
    abline(v = mu + mad(c(dat.u, -dat.u)) * c(5,6), lty = c(2,3), col = "red")
}

asymm.bounds <- function(dat, n = 6) {
    
    zz <- density(dat, n = 65536, na.rm = T)
    
    # break at point of maximum density
    mu <- zz$x[which.max(zz$y)]
    
    # calculate MAD on either side of breakpoint, centred at breakpoint
    c(mu - n * mad(dat[which(dat <= mu)], center = mu),
      mu + n * mad(dat[which(dat > mu)], center = mu))
}

hist.with.boundaries(pw.m[,,"black", "130701"], xlim = c(0, 10000))

dt <- "130701"

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "asymm-normal-", dt, ".pdf"))
                     par(mfrow = c(2,2), mar = c(3,3,3,1))
                     
                     im.b <- pw.m[,,"black", dt]
                     im.g <- pw.m[,,"grey", dt]

                     hist.with.boundaries(im.b, title = paste0("Black histogram - ", dt))
                     hist.with.boundaries(im.g, title = paste0("Grey histogram - ", dt), xlim = c(0,25000))
                     # add bounds with dark lines removed, for comparison
                     abline(v = asymm.bounds(im.g[which(im.g > 10000)]), col = "gold")
                     
                     pixel.plot(which(array(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2),
                                            dim = dim(im.b)), arr.ind = T),
                                main = paste0("Defects - ", dt, " black"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2),
                                            dim = dim(im.g)), arr.ind = T),
                                main = paste0("Defects - ", dt, " grey"))
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))
                     


####################################################################################################

# BEHAVIOUR OF ROW DEFECTS IN LOAN PANEL                                                        ####

o.plot(pw.m[, 1025, "black", "loan"] - pw.m[, 1026, "black", "loan"])
abline(h = 0, col = "red")

o.plot(pw.m[, 77, "black", "loan"] - pw.m[, 78, "black", "loan"])
abline(h = 0, col = "red")

####################################################################################################