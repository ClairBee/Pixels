
# parametric modelling of white images

library("IO.Pixels")

pw.w <- readRDS("./Other-data/Pixelwise-means-white.rds")
pw.g <- readRDS("./Other-data/Pixelwise-means-grey.rds")
pw.b <- readRDS("./Other-data/Pixelwise-means-black.rds")

pw.m <- pw.w[,,"141009"]

pixel.image(pw.m)

# may need to be carried out iteratively - EM-style algorithm?
###################################################################################################
# FITTING PROCEDURE (short version)

circ.lm <- fit.circular.lm.poly(pw.m, o = 3)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)

panel.lm <- fit.panel.lm(circ.res)
panel.res <- circ.res - panel.lm$fitted.values

qq <- lm(value ~ poly(X2, 2), data = melt(panel.res[990:1000, 1:992]))
o.plot(panel.res[992,1:992])
points(predict(qq, data.frame(X2 = c(1:992))), col = "blue", lwd = 2,  type = "l")

qq3.lower <- lm(value ~ poly(X2, 3), data = melt(panel.res[, 1:992]))
qq3.upper <- lm(value ~ poly(X2, 3), data = melt(panel.res[, 993:1996]))

o.plot(panel.res[992,1:992])
points(predict(qq3, data.frame(X2 = c(1:992))), col = "blue", lwd = 2,  type = "l")

curve.lower <- matrix(rep(predict(qq3.lower, data.frame(X2 = c(1:992))), 1996), nrow = 1996, byrow = T)
curve.upper <- matrix(rep(predict(qq3.upper, data.frame(X2 = c(993:1996))), 1996), nrow = 1996, byrow = T)

res.3 <- matrix(ncol = 1996, nrow = 1996)
res.3[,1:992] <- panel.res[,1:992] - curve.lower
res.3[,993:1996] <- panel.res[,993:1996] - curve.upper

pixel.image(res.3[,1:992])

###################################################################################################
# seems to be a slow-moving wave passing through the image, leaving 45-deg lines. Is this flicker noise?
pixel.image(panel.res)

# use detailed lowess smoothing to get underlying pattern
sm <- lowess.per.column(panel.res, span = 1/50)

cols <- c("blue", "purple" ,"red", "orange", "gold", "green", "cyan3", "skyblue")

plot(sm[1200, 993:1950], type = "l")
for (i in 1:9) {
    points(sm[1200 + i, 993:1950], type = "l", col = cols[i])
}
 
par(mfrow = c(3,3))
for (i in 0:8) {
    pacf(sm[1200 + i, 993:1950], lag = 50)
}
par(mfrow = c(1,1))

pacf(sm[1200, 993:1950], lag = 50, ci.type = "ma")
pacf(sm[1700, 993:1950], lag = 50, ci.type = "white")

plot(pacf(sm[1160, 993:1950], lag = 50, plot = F)$acf, pch = 20)
abline(h = 0)

for (i in 1:9) {
    points(pacf(sm[1160 + i, 993:1950], lag = 50, plot = F)$acf, pch = 20, col = cols[i])
}
# some evidence of a wave pattern, but v faint

spec.ar(sm[1160, 993:1950], order = 20, ylim = c(0.001, 1e+05))
for (i in 1:9) {
    spec.ar(sm[1160 + i, 993:1950], order = 20, add = T, col = cols[i])
}


###################################################################################################
# FIT CIRCULAR SPOT PATTERN, THEN PER-PANEL GRADIENT

circ.2 <- fit.circular.lm.poly(pw.m, o = 2)
circ.3 <- fit.circular.lm.poly(pw.m, o = 3)

plot(c(0:1433), predict(circ.2, data.frame(z = (c(0:1433)))), type = "l",
     xlab = "Distance from centre of panel", ylab = "Fitted value",
     main = "Profile of fitted circular spot")
points(c(0:1433), predict(circ.3, data.frame(z = (c(0:1433)))), type = "l", col = "green3")
     
# look at behaviour over single transect
z.dist <- merge(x = c(1:1996), y = c(1:1996))
z.dist$z <- sqrt((z.dist$x -1023.5)^2 + (z.dist$y - 992.5)^2)
zz <- as.matrix(cast(z.dist, x ~ y))

o.plot(pw.m[992, 993:1996])
points(c(1:1004), predict(circ.2, data.frame(z = zz[992, 993:1996])), type = "l", col = "blue")
points(c(1:1004), predict(circ.3, data.frame(z = zz[992, 993:1996])), type = "l", col = "green3")

circ.res.2 <- matrix(circ.2$residuals, ncol = 1996)
circ.res.3 <- matrix(circ.3$residuals, ncol = 1996)
o.plot(circ.res.2[992, 993:1996], main = "residuals after circular spot removed")

mean((circ.res.2)); mean((circ.res.3))
mad(circ.res.2); mad(circ.res.3)
sd(circ.res.2); sd(circ.res.3)

hist(circ.res.2, breaks = "fd", xlim = c(-500,500))
hist(circ.res.3, breaks = "fd", xlim = c(-500,500))

# fit panelwise gradients, then get MAD & compare result to direct loess-smoothing.
panel.lm.2 <- fit.panel.lm(circ.res.2)
panel.lm.3 <- fit.panel.lm(circ.res.3)

points(panel.lm.2$fitted.values[992, 993:1996], type = "l", col = "blue", lwd = 2)
points(panel.lm.3$fitted.values[992, 993:1996], type = "l", col = "green3", lwd = 2)

panel.res.2 <- circ.res.2 - panel.lm.2$fitted.values
panel.res.3 <- circ.res.3 - panel.lm.3$fitted.values

o.plot(panel.res.2[992, 993:1996])
abline(v = c(256, 512, 768), col = "red")
points(panel.res.3[992, 993:1996], col = adjustcolor("blue", alpha = 0.5), type = "l")

mean((panel.res.2)); mean((panel.res.3))
mad(panel.res.2); mad(panel.res.3)
sd(panel.res.2); sd(panel.res.3)

# compare to direct loess smoothing
smoo <- lowess.per.column(pw.m)
smoo.res <- pw.m - smoo
mean(smoo.res); mad(smoo.res); sd(smoo.res)

points(smoo.res[992, 993:1996], col = adjustcolor("green3", alpha = 0.5), type = "l")

# how consistent are per-column residuals vs lowess?
col.res.clm.u <- apply(panel.res.2[,993:1996], 1, mean) 
col.res.clm.l <- apply(panel.res.2[,1:992], 1, mean) 

o.plot(col.res.clm.l, ylim = c(-500,500))
o.plot(col.res.clm.u, add = T, col = adjustcolor("blue", alpha = 0.5))

hist(col.res.clm.u, breaks = "fd", col = adjustcolor("blue", alpha = 0.5), xlim = c(-500,500))
hist(col.res.clm.l, breaks = "fd", col = adjustcolor("green3", alpha = 0.5), add = T)

col.res.low.u <- apply(smoo.res[,993:1996], 1, mean) 
col.res.low.l <- apply(smoo.res[,1:992], 1, mean) 

o.plot(col.res.low.l, ylim = c(-500,500))
o.plot(col.res.low.u, add = T, col = adjustcolor("blue", alpha = 0.5))

o.plot(apply(panel.res.2[,993:1996], 1, mad), ylim = c(0,500), main = "MAD of parametric residuals")
o.plot(apply(panel.res.2[,1:992], 1, mad), add = T, col = adjustcolor("blue", alpha = 0.5))

o.plot(apply(smoo.res[,1:992], 1, mad), add = T, col = adjustcolor("green3", alpha = 0.5))
o.plot(apply(smoo.res[,993:1996], 1, mad), add = T, col = adjustcolor("red3", alpha = 0.5))


hist(col.res.low.l, breaks = "fd", col = adjustcolor("blue", alpha = 0.5))
hist(col.res.low.u, breaks = "fd", col = adjustcolor("green3", alpha = 0.5), add = T)
# lowess smoothing may be over-smoothing? Takes out what may be bad pixels in TL panel.

###################################################################################################
# get 'sickness score' per pixel based on MAD
sickness.2lm <- panel.res.2 / mad(panel.res.2)
sickness.smoo <- smoo.res / mad(smoo.res)

hist(sickness.2lm, breaks = "fd", prob = T, ylim = c(0,0.0005), main = "Sickness scores from parametric model")
lines(c(-400:400), dnorm(c(-400:400), mean = mean(sickness.2lm), sd = sd(sickness.2lm)), lwd = 3, col = "cornflowerblue")

hist(sickness.smoo, breaks = "fd", prob = T, ylim = c(0,0.0005), main = "Sickness scores from direct lowess smoothing")
lines(c(-400:400), dnorm(c(-400:400), mean = mean(sickness.2lm), sd = sd(sickness.2lm)), lwd = 3, col = "cornflowerblue")

###################################################################################################
# should compare results with different degrees of lowess smoothing

smoo.1 <- lowess.per.column(pw.m, f = 1)

o.plot(pw.m[992, 1:996])
points(lowess(pw.m[992, 1:996], f = 1), type = "l", col = "red", lwd = 2)       # min. smoothing possible
points(lowess(pw.m[992, 1:996], f = 0.75), type = "l", col = "purple", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/2), type = "l", col = "blue", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/3), type = "l", col = "seagreen", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/4), type = "l", col = "green3", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/5), type = "l", col = "green", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/6), type = "l", col = "greenyellow", lwd = 2)     # closest fit?
points(lowess(pw.m[992, 1:996], f = 1/8), type = "l", col = "yellow", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/10), type = "l", col = "gold", lwd = 2)
points(lowess(pw.m[992, 1:996], f = 1/12), type = "l", col = "orange", lwd = 2)

spans <- c(1, 3/4, 1/2, 1/3, 1/4, 1/5, 1/6, 1/8, 1/10, 1/12, 1/15)
span.labels <- c("1", "3/4", "1/2", "1/3", "1/4", "1/5", "1/6", "1/8", "1/10", "1/12", "1/15")


res <- data.frame(mean = numeric(), abs.mean = numeric(), sd = numeric(), mad = numeric())

for (i in 1:length(spans)) {
    res[i,] <- c(mean(abs(pw.m[992, 1:996] - lowess(pw.m[992, 1:996], f = spans[i])$y)),
                 mean((pw.m[992, 1:996] - lowess(pw.m[992, 1:996], f = spans[i])$y)),
                 sd(pw.m[992, 1:996] - lowess(pw.m[992, 1:996], f = spans[i])$y),
                 mad(pw.m[992, 1:996] - lowess(pw.m[992, 1:996], f = spans[i])$y))
}
res
scale.res <- scale(res, center = T, scale = F)

o.plot(scale.res[,1], ylim = c(-40, 140), xaxt = "none")
axis(1, at = c(1:11), labels = span.labels)
o.plot(scale.res[,2], add = T, col = "red")
o.plot(scale.res[,3], add = T, col = "blue")
o.plot(scale.res[,4], add = T, col = "green3")
legend("topright", col = c("black", "red", "blue", "green3"), pch = 20, bty = "n",
       legend = c("Mean residual", "Mean abs. residual", "Residual SD", "Residual MAD"))

###################################################################################################

# possible alternative: use loess-smoothing in both directions and somehow combine the scores to get 'atypical' points

# plots created by plot.lm:
#   - residuals vs fitted
#   - normal Q-Q
#   - scale-location
#   - residuals

circ.fitted <- matrix(circ.lm$fitted.values, ncol = 1996)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)

pixel.image(circ.res)
hist(circ.res, breaks = "fd", xlim = c(-5000,5000))

qq <- fit.panel.lm(circ.res)
pixel.image(qq$fitted.values)

panel.res <- circ.res - qq$fitted.values
pixel.image(panel.res)

s.hist(panel.res, main = "")
o.plot(upper.res[992,])
abline(h = 0, col = "red", lwd = 2)
abline(h = mad(upper.res[992,]) * c(-1,1), col = "blue", lwd = 2)
abline(h = mad(upper.res[992,]) * 2 * c(-1,1), col = "purple", lwd = 2)

qq <- lowess.per.column(panel.res)

pixel.image(qq)
s.hist(qq)

###################################################################################################

# start by fitting per-panel gradients
# ignoring possible half-panel breaks for now


# then fit circular pattern

# check residuals


# MSE as measure of fit