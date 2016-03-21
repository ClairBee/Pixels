
library("IO.Pixels")

# original plan was as outlined below. However, more artefacts are visible in the data

###################################################################################################

# trying to identify bad pixels by deviations from dominant amplitude along columns in dark images
# get amplitude by Loess smoothing/differencing
# use NMAD/SD to define bounds on acceptable values
# flag as bad pixels any that fall outside of that range

# compare to 'official' bad pixel map
# compare to thresholded definitions

###################################################################################################

load.images(160314, "black")
pw.m <- pixelwise.mean(b.160314)

###################################################################################################
# tried using simple differencing - will flag up 'return' point as an error as well as actual error.
{
diffs.lower <- pw.m[, 1:991] - pw.m[, 2:992]
diffs.upper <- pw.m[, 993:1995] - pw.m[,994:1996]

o.plot(diffs.lower[993,], ylim = c(-200, 200))

cw.sd.diff.l <- apply(abs(diffs.lower), 1, sd)
cw.sd.diff.u <- apply(abs(diffs.upper), 1, sd)

cw.m.diff.l <- apply(abs(diffs.lower), 1, mean)
cw.m.diff.u <- apply(abs(diffs.upper), 1, mean)

cw.mad.diff.l <- apply(abs(diffs.lower), 1, mad)
cw.mad.diff.u <- apply(abs(diffs.upper), 1, mad)

# look at single column to start with
abs.sd <- sd(abs(diffs.lower[993,]))                # 28.6771
abs.mean <- mean(abs(diffs.lower[993,]))            # 89.45303
abs.mad <- mad(abs(diffs.lower[993,]))              # 19.34793
abs.outlier.u <- quantile(abs(diffs.lower[993,]), 0.75) + 1.5 * IQR(abs(diffs.lower[993,]))
abs.outlier.l <- quantile(abs(diffs.lower[993,]), 0.25) - 1.5 * IQR(abs(diffs.lower[993,]))


o.plot(abs(diffs.lower[993,]), ylim = c(0, 200))

abline(h = abs.mean, col = "red", lwd = 2)       #   0.3074168
abline(h = abs.mean + (abs.sd * c(-1,1)), col = "blue", lwd = 2)
abline(h = abs.mean + (abs.sd * 1.96 * c(-1,1)), col = "green3", lwd = 2)

abline(h = abs.outlier.u, col = "purple", lwd = 2)
abline(h = abs.outlier.l, col = "purple", lwd = 2)


# could get row-wise SD to give bounds?
bp <- which(abs(diffs.lower[993,]) > (abs.mean + (abs.sd * 1.96)), arr.ind = T)

o.plot(diffs.lower[993,])
points(cbind(bp, diffs.lower[993, c(bp)], ncol = 2), col = "red")

o.plot(pw.m[993,1:992])
points(cbind(bp, pw.m[993, c(bp)], ncol = 2), col = "red")
}

###################################################################################################
# now with Loess smoothing to get centre line
c <- 992; r <- c(1, 992)
transect <- pw.m[c, r[1]:r[2]]

#
{

o.plot(transect)

# jump at halfway point, maybe another at quarters? Need to look into this.
abline(v = 992-(1024 * c(0.25, 0.5, 0.75)), col = "red")

smoo <- lowess(transect, f = 1/15)$y
points(smoo, type = "l", col = "red", lwd = 2)

res <- transect - smoo
abs.res <- abs(res)

abs.mean <- mean(abs.res)       # 45.14857
abs.sd <- sd(abs.res)           # 21.10757
abs.outlier.u <- quantile(abs.res, 0.75) + (1.5 * IQR(abs.res))
abs.outlier.l <- quantile(abs.res, 0.25) - (1.5 * IQR(abs.res))

bp <- which(abs.res > (abs.mean + 1.96 * abs.sd))

plot(transect, type = "l", col = adjustcolor("grey", alpha = 0.5))
points(transect, pch = 20, cex = 0.7)

points(smoo + abs.mean, type = "l", col = "blue", lwd = 2)
points(smoo - abs.mean, type = "l", col = "blue", lwd = 2)

points(smoo - abs.mean - (1.96 * abs.sd), type = "l", col = "purple", lwd = 2)
points(smoo + abs.mean + (1.96 * abs.sd), type = "l", col = "purple", lwd = 2)

points(cbind(bp, transect[bp], ncol = 2), col = "red")

# how do those pixels behave in white/grey panels?
#load.images(160314, "grey"); pw.g <- pixelwise.mean(g.160314)
#load.images(160314, "white"); pw.w <- pixelwise.mean(w.160314)

o.plot(pw.g[c,1:992])
points(cbind(bp, pw.g[c,1:992][bp], ncol = 2), col = "red")

o.plot(pw.w[c,1:992])
points(cbind(bp, pw.w[c,1:992][bp], ncol = 2), col = "red")
}

# try using different measures against Loess smoothing to identify residuals
{
    
}
###################################################################################################
# data is highly asymmetric - divide points according to those above & those below centre
{
transect <- pw.m[992,1:992]

smoo <- lowess(transect, f = 1/15)$y
points(smoo, type = "l", col = "red", lwd = 2)

res <- transect - smoo

o.plot(res)
points(cbind(c(1:992)[res >= 0], res[res >= 0]), pch = 20, col = "purple")
points(cbind(c(1:992)[res < 0], res[res < 0]), pch = 20, col = "blue")
abline(h = 0, col = "red")

res.high <- res[res >= 0]
res.low <- res[res < 0]

mean(res.high)
sd(res.high)

mean(res.low)
sd(res.low)
}
# abs. value of point is affected by preceding difference,
# and by accuracy of Loess fit
# maybe try using diffs, normalised by residuals?

###################################################################################################
# rather than neighbouring points, try difference of next-but-one neighbours
{
transect <- pw.m[992,1:992]
diff1 <- transect[1:991] - transect[2:992]
diff2 <- transect[1:990] - transect[3:992]

plot(abs(diff2), pch = 20, cex = 0.7)

m <- mean(abs(diff2))
sd <- sd(abs(diff2))

abline(h = m, col = "red", lwd = 2)
abline(h = m + (1.96 * c(-1, 1) * sd), col = "blue", lwd = 2)

bp <- which(abs(diff2) > m + (1.96 * sd), arr.ind = T)

points(cbind(c(1:992)[bp], abs(diff2)[bp]), col = "red")
points(x = 552, y = abs(diff2)[552], col = "blue")

o.plot(transect)
points(x = 552, y = transect[552], col = "red")

ratio <- diff2 / transect[1:990]
o.plot(ratio)

# try MA
smoo <- filter(transect, filter = c(1/3, 1/3, 1/3), sides = 2)
o.plot(smoo)

smoo.diff <- smoo[1:991] - smoo[2:992]
o.plot(smoo.diff)
o.plot(abs(smoo.diff))

m <- mean(abs(smoo.diff), na.rm = T)
sd <- sd(abs(smoo.diff), na.rm = T)
abline(h = m, col = "red", lwd = 2)
abline(h = m + (1.96 * c(-1,1) * sd), col = "blue", lwd = 2)
}

###################################################################################################
# row-wise SD/MAD of oscillations
