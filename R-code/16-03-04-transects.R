
library("IO.Pixels")

load.images(150828, "black")

#=================================================================================
# plot transect across mean/sd value. Use loess smoothing to establish outliers?
pw.m <- pixelwise.mean(b.150828)
pw.sd <- pixelwise.sd(b.150828)

mean.sd <- sd(pw.m)

# get smoothed mean & rough confidence interval
# better to use proper loess confidence interval, if possible
smoothed.mean <- lowess(pw.m[1, c(1:992)])
lower <- list(x = smoothed.mean$x, y = smoothed.mean$y - mean.sd)
upper <- list(x = smoothed.mean$x, y = smoothed.mean$y + mean.sd)

plot(pw.m[1, c(1:992)], type = "o", pch = 20)
lines(lowess(pw.m[1, c(1:992)], f = 1/3), col = "red", lwd = 3)
lines(lower, col = "red3", lwd = 2)
lines(upper, col = "red3", lwd = 2)

bright <- which(pw.m[1,c(1:992)] > upper$y)
dim <- which(pw.m[1,c(1:992)] < lower$y)

#=================================================================================
minidat <- b.150828[c(1150:1278), c(1868:1996),]

pw.m.b.150828 <- pixelwise.mean(b.150828)

pixel.image(minidat[,,1], break.levels = sd.levels(pw.m.b.150828))
pw.m <- pixelwise.mean(minidat)

plot(pw.m[1, ], type = "o", pch = 20)

# 'jitter' clearly visible



