
library("IO.Pixels"); library("CB.Misc")

# LINE DETECTION ALGORITHM UNSATISFACTORY

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7")

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

####################################################################################################

# FUNCTIONS                                                                                     ####

# find MAD thresholds based on two half-normal densities, conjoined at modal density
asymm.bounds <- function(dat, n = 6) {
    
    zz <- density(dat, n = 65536, na.rm = T)
    
    # break at point of maximum density
    mu <- zz$x[which.max(zz$y)]
    
    # calculate MAD on either side of breakpoint, centred at breakpoint
    c(mu - n * mad(dat[which(dat <= mu)], center = mu),
      mu + n * mad(dat[which(dat > mu)], center = mu))
}

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "", ...) {
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "", ...)
    abline(v = asymm.bounds(dat, n = 6), lty = 2, col = "red")
    abline(v = asymm.bounds(dat, n = 5), lty = 3, col = "red")
    
    abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JohnsonFit(dat[!is.na(dat)])), col = "cyan3", lty = c(2,3,3,2))
}

####################################################################################################

# THRESHOLDING BY ASYMMETRIC MAD                                                                ####

fpath <- "./Image-plots/thresholds/"


px <- list()

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "thresholds-", dt, ".pdf"), width = 14)
                     par(mfrow = c(2, 4), mar = c(3,3,3,1))
                     
                     im.b <- pw.m[,,"black", dt]; im.g <- pw.m[,,"grey", dt]
                     res.b <- md7[,,"black", dt]; res.g <- md7[,,"grey", dt]
                     
                     hist.with.boundaries(im.b, title = paste0(dt, " - black values"))
                     hist.with.boundaries(im.g, title = paste0(dt, " - grey values"), xlim = c(0,25000))
                     
                     hist.with.boundaries(res.b, title = paste0(dt, " - black residuals"), xlim = c(-1000,1000))
                     hist.with.boundaries(res.g, title = paste0(dt, " - grey residuals"), xlim = c(-1000,1000))
                     
                     n.val.b <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2)))
                     n.val.g <- length(which(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2)))
                     n.res.b <- length(which(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2)))
                     n.res.g <- length(which(findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     n.all.b <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2) | findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2)))
                     n.all.g <- length(which(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2) | findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     n.all <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2) | findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2) |
                                               findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2) | findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     px[[dt]] <<- data.frame(n.val.b, n.val.g, n.res.b, n.res.g, n.all.b, n.all.g, n.all)
                     
                     pixel.plot(which(array(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2),
                                            dim = dim(im.b)), arr.ind = T), 
                                main = paste0(dt, " - black values (", n.val.b, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2),
                                            dim = dim(im.g)), arr.ind = T),
                                main = paste0(dt, " - grey values  (", n.val.g, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2),
                                            dim = dim(res.b)), arr.ind = T), 
                                main = paste0(dt, " - black residuals (", n.res.b, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2),
                                            dim = dim(res.g)), arr.ind = T), 
                                main = paste0(dt, " - grey residuals (", n.res.g, ")"))
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))
                     

df <- rbind.fill(px)
rownames(df) <- dimnames(pw.m)[[4]]

write.csv(df, paste0(fpath, "px-identified.csv"), quote = F)

####################################################################################################

# CLASSIFICATION BY ASYMMETRIC MAD                                                              ####

dt <- "160314"

# get all abnormal pixels, by mean value & median-smoothed residual
bpx <- which(array(findInterval(pw.m[,,"black", dt], asymm.bounds(pw.m[,,"black", dt])) %in% c(0,2) |
                       findInterval(pw.m[,,"grey", dt], asymm.bounds(pw.m[,,"grey", dt])) %in% c(0,2) |
                       findInterval(md7[,,"black", dt], asymm.bounds(md7[,,"black", dt])) %in% c(0,2) |
                       findInterval(md7[,,"grey", dt], asymm.bounds(md7[,,"grey", dt])) %in% c(0,2), dim = dim(pw.m[,,"black", dt])),
             arr.ind = T)

os <- data.frame(x = bpx[,1],
                 y = bpx[,2],
                 b = pw.m[,,"black", dt][bpx],
                 g = pw.m[,,"grey", dt][bpx],
                 res.b = md7[,,"black", dt][bpx],
                 res.g = md7[,,"grey", dt][bpx])


abline(h = asymm.bounds(md7[,,"grey", dt]), col = "green3", lty = 2)
abline(v = asymm.bounds(md7[,,"black", dt]), col = "green3", lty = 2)


# cut data midway between bounds & max values
class.boundaries <- function(dat) {
    
    rng <- c(floor(min(dat, na.rm = T)), ceiling(max(dat, na.rm = T)) + 1)
    inner.bounds <- asymm.bounds(dat)
    
    vb <- inner.bounds[2] + ((rng[2] - inner.bounds[2]) * c(0.25, 0.5, 0.75))

    sort(c(rng, inner.bounds, vb))
}

.smoothScatter(os$res.b, os$res.g, nrpoints = Inf, xlab = "Black residuals", ylab = "Grey residuals", main = paste0(dt))
abline(0,1, col = "red")

abline(v = class.boundaries(md7[,,"black", dt]), col = "skyblue", lty = 2)
abline(h = class.boundaries(md7[,,"grey", dt]), col = "skyblue", lty = 2)

table(findInterval(md7[,,"black", dt], class.boundaries(md7[,,"black", dt])), useNA = "ifany")

os$class.b <- findInterval(os$b, class.boundaries(os$b))
os$class.g <- findInterval(os$g, class.boundaries(os$g))

table("black" = os$class.b, "grey" = os$class.g)


cor(os[!is.na(os$res.b), c("res.b", "g")])
cor(os[!is.na(os$res.g), c("res.g", "b")])





cor(os[!is.na(os$res.b) & os$g < 65535 & os$b > 0, c("res.b", "res.g")])

smoothScatter(os$res.b, os$res.g, nrpoints = Inf)
abline(0,1, col = "red")

smoothScatter(os$b, os$res.b, nrpoints = Inf, xlab = "Values", ylab = "Residuals", main = paste0(dt, " - black"))
abline(a = -5000, b = 1, col = "red")

smoothScatter(os$g, os$res.g, nrpoints = Inf, xlab = "Values", ylab = "Residuals", main = paste0(dt, " - grey"))
abline(a = -17500, b = 1, col = "red")

plot(os, lower.panel = panel.cor)

####################################################################################################

# CHECK LINEARITY ISSUES VS OFFSET ERRORS                                                       ####

# decide whether to include as separate category
# are most nonlinear pixels already identified using residual/extreme-value approach?

wlm <- fit.w.lm(pw.m[,,,dt])

.smoothScatter(wlm$df$g.x, wlm$df$g.y)

hist.with.boundaries(wlm$df$res, xlim = c(-2000,2000))

nlpx <- wlm$df[findInterval(wlm$df$res, asymm.bounds(wlm$df$res)) %in% c(0,2),]
             
.smoothScatter(nlpx$res.b, nlpx$res.g)
abline(0,1)

pixel.plot(nlpx, col = "cyan3")
pixel.plot(bpx, col = "orange")

px <- merge(data.frame(nlpx, type = "nl"), data.frame(os, type = "os"), by = c("x", "y"), all = T)
px$b <- apply(px[,c("b.x", "b.y")], 1, mean, na.rm = T)
px$g <- apply(px[,c("g.x", "g.y")], 1, mean, na.rm = T)

px$res.b <- md7[,,"black", dt][as.matrix(px[,c("x", "y")])]
px$res.g <- md7[,,"grey", dt][as.matrix(px[,c("x", "y")])]

.smoothScatter(px[is.na(px$type.y) & !is.na(px$type.x), c("res.b", "res.g")])
abline(v = class.boundaries(md7[,,"black", dt]), lty = 2)
abline(h = class.boundaries(md7[,,"grey", dt]), lty = 2)

# check number of nonlinear pixels NOT identified as offset errors in all images



####################################################################################################

# BEHAVIOUR OF ROW DEFECTS IN LOAN PANEL                                                        ####

o.plot(pw.m[, 1025, "black", "loan"] - pw.m[, 1026, "black", "loan"])
abline(h = 0, col = "red")

o.plot(pw.m[, 77, "black", "loan"] - pw.m[, 78, "black", "loan"])
abline(h = 0, col = "red")

####################################################################################################

####################################################################################################

# DISCONTINUED                                                                                  ####

# JOHNSON DISTRIBUTIONS OF ALL BLACK & GREY IMAGES

{fpath <- "./Image-plots/histograms/"

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
}
####################################################################################################

# HALF-NORMAL DISTRIBUTION FOR HISTOGRAMS 
{
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
lines(0:5000, dhalfnorm(0:5000, sd2theta(sig.u)), col = "red", lwd = 2)}

# images with strong upward drift are poorly fitted by half-normal. Hey ho.

####################################################################################################
