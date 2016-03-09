
library("IO.Pixels")

o.plot <- function(data, new = F, ...) {
    if (new) {
        plot(data, type = "o", pch = 20, cex = 0.7, ...)
    } else {
        points(data, type = "o", pch = 20, cex = 0.7, ...)
    }
}

load.images(150828, "black")
pw.m.b.150828 <- pixelwise.mean(b.150828)
pw.sd.b.150828 <- pixelwise.sd(b.150828)
lvls <- sd.levels(pw.m.b.150828)

x <- (1023, 1150); y <- c(1,992)
panel <- pw.m.b.150828[x[1]:x[2], y[1]:y[2]]
transect <- panel[4,]


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


# apply to whole panel
smoothed <- do.call(rbind, lapply(apply(panel, 1, lowess, f = 1/15), "[[", 2))

pixel.image(smoothed, break.levels = sd.levels(panel))
pixel.image(panel)

# check that result is same as for transect: smoothing applied in correct direction?
o.plot(panel[4,], col = "black", new = T)
o.plot(smoothed[4,], col = adjustcolor("red", alpha = 0.5))

panel.res <- abs(panel - smoothed)
# ok up to this point

z <- apply(panel.res, 2, movingFun, n = 5, fun = median, type = 'around')
zz <- matrix(NA, ncol = ncol(panel.res), nrow = nrow(panel.res))
for (i in 1:ncol(panel.res)) {
    
}
z[4,] == res.amp

o.plot(z[4,], new = T)
