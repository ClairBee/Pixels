
library("IO.Pixels"); library("CB.Misc")

load.pixel.means.2048()

# rather than as a model for thresholding, consider using circular spot to to describe shape of data

# alternative: allow polynomial in x and y separately, rather than combined into hypoteneuse
# could also fit a ring to the data? Or maybe a second circle after quadratic trend?

####################################################################################################

# PARAMETRIC SPOT-AND-PANEL MODEL TO DESCRIBE SHAPE                                             ####

# fit a circular spot with centre at panel midpoint
{
    spot.lm <- rlm(gv ~ poly(sqrt((x - 1024)^2), 2) + poly(sqrt((y - 1024)^2),2), 
                   setNames(melt(pw.m[,,"black", "160430"]), nm = c("x", "y", "gv")))
    fitted <- spot.res <- array(dim = dim(pw.m[,,"black", "160430"]))
    fitted[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- spot.lm$fitted.values
    spot.res[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- spot.lm$residuals
    
    pixel.image(spot.res)
    pixel.image(fitted)
    
    # curve profile
    lines(fitted[1024,], col = "red")
}

# split data into subpanels, fit individual linear models
{
    sp <- array(pw.m[,,"black", "160430"], dim = c(128, 16, 1024, 2))
    panel.lm <- array(dim = dim(sp))
    coeffs <- c()
    
    for (ul in 1:2) {
        for (p in 1:16) {
            lm.fit <- rlm(gv ~ x + y,
                          setNames(melt(sp[,p,,ul]), nm = c("x", "y", "gv")))
            panel.lm[,p,,ul][which(!is.na(sp[,p,,ul]), arr.ind = T)] <- lm.fit$fitted.values
            coeffs <- rbind(coeffs, coef(lm.fit))
        }
    }
    panel.lm <- array(panel.lm, dim = c(2048, 2048))
    panel.res <- pw.m[,,"black", "160430"] - panel.lm
    pixel.image(panel.res)
    pixel.image(panel.lm)
    draw.panels.2048()
}

coeffs

####################################################################################################

# CIRCULAR SPOT W/SP OFFSETS                                                                    ####

p <- panel.edges(left.crop = 0, upper.crop = 0, x.dim = 2048, y.dim = 2048)

df <- setNames(melt(pw.m[,,"black", "160430"]), nm = c("x", "y", "gv"))
df$sp <- factor(apply(cbind(c("L", "U")[findInterval(df$y, p$y-0.5)],
                            formatC(findInterval(df$x, p$x-0.5), width = 2, flag = "0")),
                      1, paste, collapse = ""))
    
# allow extra offset for subpanels
spot.lm <- rlm(gv ~ poly(sqrt((x - 1024)^2 + (y - 1024)^2), 2) + sp, df)
spot.fitted <- spot.res <- array(dim = dim(pw.m[,,"black", "160430"]))

spot.res[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- spot.lm$residuals
spot.fitted[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- spot.lm$fitted.values

pixel.image(spot.fitted)
pixel.image(spot.res)

o.plot(pw.m[1090,,"black", "160430"])
lines(spot.fitted[1090,], col = "red")

o.plot(diag(pw.m[,,"black", "160430"]))
lines(diag(spot.fitted), col = "red")

coef(spot.lm)

# fit panel model to residuals
{
    sp <- array(spot.res, dim = c(128, 16, 1024, 2))
    panel.lm <- array(dim = dim(sp))
    coeffs <- c()
    
    for (ul in 1:2) {
        for (p in 1:16) {
            lm.fit <- rlm(gv ~ x + y,
                          setNames(melt(sp[,p,,ul]), nm = c("x", "y", "gv")))
            panel.lm[,p,,ul][which(!is.na(sp[,p,,ul]), arr.ind = T)] <- lm.fit$fitted.values
            coeffs <- rbind(coeffs, coef(lm.fit))
        }
    }
    panel.lm <- array(panel.lm, dim = c(2048, 2048))
    panel.res <- spot.res - panel.lm
    
    pixel.image(panel.res)
    s.hist(panel.res)
    pixel.image(panel.lm)
}

####################################################################################################

# QUADRATIC TREND MODEL                                                                         ####

quad.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y), 
               setNames(melt(pw.m[,,"black", "160430"]), nm = c("x", "y", "gv")))

quad.fitted <- quad.res <- array(dim = dim(pw.m[,,"black", "160430"]))
quad.fitted[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- quad.lm$fitted.values
quad.res[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- quad.lm$residuals

pixel.image(quad.res)
pixel.image(quad.fitted)
coef(quad.lm)
summary(quad.lm)

hist(quad.res)

####################################################################################################

# CUBIC TREND MODEL                                                                             ####

# need to cut down image. Retain odd columns only
sparse.im <- setNames(melt(pw.m[,,"black", "160430"]), nm = c("x", "y", "gv"))
sparse.im <- sparse.im[sparse.im$x %% 2 == 1 & sparse.im$y %% 2 == 1,]

cube.lm <- rlm(gv ~ x + y + I(x^2) + I(y^2) + I(x *y) + I(x^3) + I(y^3), 
               setNames(melt(pw.m[,,"black", "160430"]), nm = c("x", "y", "gv")))

cube.fitted <- cube.res <- array(dim = dim(pw.m[,,"black", "160430"]))
cube.fitted[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- cube.lm$fitted.values
cube.res[which(!is.na(pw.m[,,"black", "160430"]), arr.ind = T)] <- cube.lm$residuals

pixel.image(cube.res)
pixel.image(cube.fitted)
coef(cube.lm)
summary(cube.lm)

s.hist(cube.res)

sd(cube.res, na.rm = T)
sd(quad.res, na.rm = T)
sd(pw.m[,,"black", "160430"], na.rm = T)

####################################################################################################

# UNDERSTANDING PARAMETERS OF QUADRATIC TREND MODEL                                             ####

# do parameters have natural interpretation?