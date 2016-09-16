
library("IO.Pixels"); library("CB.Misc")

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")

# should thresholding of linear residuals be done on MAD or on absolute threshold?
fpath <- "./Image-plots/linear-res/"

####################################################################################################

# GET ALL LINEAR MODELS & RESIDUALS                                                             ####

# fit all linear models
linear <- invisible(apply(pw.m, 4, fit.w.lm, res.only = T))

saveRDS(linear, paste0(fpath, "models.rds"))

# get SD/MAD of all residuals
m <- lapply(linear, function(mod) asymmetric.mad(mod$df$res, n = 1))
s <- lapply(linear, function(mod) sd(mod$df$res, na.rm = T))

cbind(do.call("rbind", m), unlist(s))

summ <- data.frame(mad = unlist(lapply(linear, function(mod) mad(mod$df$res, na.rm = T))),
                   r2 = unlist(lapply(linear, "[", "r2")))

loan.res <- array(linear$"loan"$df$res, dim = c(2048, 2048))
lh <- hist(loan.res, breaks = "fd", ylim = c(0,10000), xlim = c(-1000,1000))
hist(loan.res[1500:2048, 1:200], breaks = lh$breaks, add = T, border = "red")
abline(v = c(-500,500), col = "orange")

pixel.plot(which(loan.res > 500, arr.ind = T), col = "red")
points(which(loan.res < -500, arr.ind = T), col = "blue", cex = 0.4, pch = 15)

# get correlation of grey & fitted values for all images

qq <- invisible(lapply(linear[c("130613", "141009", "loan", "MCT225")], 
                       function(ff) ff$df$res))
mad(unlist(qq), na.rm = T)

####################################################################################################

# LOAN PANEL - NEW APPROACH TO THRESHOLDING?                                                    ####

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "", JF = F, ...) {
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "", ...)
    abline(v = asymmetric.mad(dat, n = 6), lty = 2, col = "red")
    
}

dt <- "loan"
dt <- "140129"

linear <- fit.w.lm(pw.m[,,,dt])

# actual MAD of each half of distribution
mu <- modal.density(linear)
a.mad <- abs(mu - asymmetric.mad(linear, n = 1))

#hist(linear, breaks = "fd", xlim = c(-2000,2000), main = paste0("linear residuals - ", dt, " panel"), xlab = "", ylab = "")
#abline(v = modal.density(linear), col = "orange", lwd = 2)
#abline(v = c(-6,6) * a.mad[2], col = "red")
#abline(v = c(-6,6) * a.mad[1], col = "orange")

#pixel.plot(which(linear < -6 * a.mad[2], arr.ind = T), col = "blue")



jpeg(paste0("./Notes/01_asymm-thresholds/",dt,"-linear-residuals.jpg"))
pixel.image(linear)
dev.off()

jpeg(paste0("./Notes/01_asymm-thresholds/",dt,"-linear-residuals-hist.jpg"))
hist.with.boundaries(linear, xlim = c(-2000,2000))
#abline(v = mu, col = "orange", lwd = 2)
abline(v = mu + c(-6,6) * a.mad[2], col = "cyan3")
dev.off()

jpeg(paste0("./Notes/01_asymm-thresholds/",dt,"-linear-threshold-shift.jpg"))
pixel.plot(which(linear < mu - 6 * a.mad[2], arr.ind = T))
points(which(linear < mu - 6 * a.mad[1], arr.ind = T), col = "red", pch = 15)
dev.off()

dt <- "loan"
hh <- hist(linear, breaks = "fd", xlim = c(-2000,2000))
hist(linear[which(linear < mu - 6 * a.mad[2], arr.ind = T)], add = T, breaks = hh$breaks, col = "red", border = "red")

sc <- shading.corrected(pw.m[,,,"loan"])
hh <- hist(sc, breaks = "fd", xlim = c(15000,25000), ylim = c(0,5000))
hist(sc[which(linear < mu - 6 * a.mad[2], arr.ind = T)], add = T, breaks = hh$breaks, col = "red", border = "red")
