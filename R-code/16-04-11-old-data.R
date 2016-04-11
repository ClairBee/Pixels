
# TESTING PARAMETRIC MODEL OVER OLD DATA TO LOOK FOR HOT LINES

library("IO.Pixels")

# panel offsets are different - need to adjust accordingly

########################################################################################################

# load old data, crop by 4px on RHS to give standard 1996x1996 image
b.130613 <- t(readTIFF("./Other-data/Old-data/130613-bad-pixel-map/BadPixelMapBlack.tif", as.is = T)[1996:1,,drop = F])[1:1996,]
g.130613 <- t(readTIFF("./Other-data/Old-data/130613-bad-pixel-map/BadPixelMapGrey.tif", as.is = T)[1996:1,,drop = F])[1:1996,]
w.130613 <- t(readTIFF("./Other-data/Old-data/130613-bad-pixel-map/BadPixelMapWhite.tif", as.is = T)[1996:1,,drop = F])[1:1996,]

pixel.image(b.130613)
pixel.image(g.130613)
pixel.image(w.130613)

########################################################################################################

# SQUIRCULAR MODEL                                                                                  ####

spot <- he.spot.lm(w.130613, n = 2, o = 2, robust = T)
spot.res <- matrix(spot$residuals, ncol = 1996)

pixel.image(matrix(spot$fitted.values, ncol = 1996))

panel <- panel.lm(spot.res, "x + y", robust = T, upper.pad = 24, left.pad = 24)
res <- spot.res - panel$fitted.values

pixel.image(res)
o.plot(res[992,])

# not sure how that's happened. But it's very odd.

pixel.image(w.130613)
# need to fit separate model to top & bottom of image?

