
library("IO.Pixels")

load.images(150828, "black")
pw.m.b.150828 <- pixelwise.mean(b.150828)
pw.sd.b.150828 <- pixelwise.sd(b.150828)
lvls <- sd.levels(pw.m.b.150828)

x <- c(1023, 1150); y <- c(1,992)
panel <- pw.m.b.150828[x[1]:x[2], y[1]:y[2]]
transect <- panel[4,]

#########################################################################################
# Looking at patterns of oscillations across panels

# look at single-column transect: smooth data using Lowess
# need to look at impact of different smoothing rates.
smoo <- lowess(transect, f = 1/15)$y

plot(transect, type = "o", pch = 20, cex = 0.7)
points(smoo, type = "o", pch = 20, cex = 0.7, col = adjustcolor("red", alpha = 0.5))

res <- transect - smoo
plot(res, type = "o", pch = 20, cex = 0.7)
mean(res)

res.amp <- movingFun(abs(res), 5, median, 'around')
plot(res.amp, type = "o", pch = 20, cex = 0.7)

#########################################################################################

# apply smoother to columns of whole panel
# check effect of changes to smoothing filter
smoothed <- do.call(rbind, lapply(apply(panel, 1, lowess, f = 1/15), "[[", 2))

# get residuals
panel.res <- panel - smoothed       
hist(panel.res, breaks = "fd", xlim = c(-1,1) * 200)

pixel.image(smoothed, break.levels = iqr.levels(panel))
pixel.image(panel)

# check that result is same as for transect: smoothing applied in correct direction?
o.plot(panel[4,], col = "black")
o.plot(smoothed[4,], col = adjustcolor("red", alpha = 0.5), add = T)

o.plot(panel.res[4,]); abline(h = 0, col = "red3", lwd = 2)

#########################################################################################
# check ?random subsamples of points along row; how do mean & SD change?
m <- matrix(panel.res[4,], ncol = 8, nrow = 124)

mean(sample(panel.res[4,], length(panel.res[4,]), replace = T))

res.means <- apply(panel.res, 1, mean)
res.medians <- apply(panel.res, 1, median)
res.sds <- apply(panel.res, 1, sd)

o.plot(res.means)
o.plot(res.medians)
o.plot(res.sds)

# all columnwise residual means are positive b/c in black channel, greater range for high values
# columnwise residual medians also mostly positive. Not sure why this is.
# No readily apparent pattern across columns
#-----------------------------------------------------------------------------
# 
cor(panel.res[1,], panel[1,])

