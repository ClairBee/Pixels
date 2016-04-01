
library("IO.Pixels")
library(MASS)               # needed for robust linear regression
library(SuppDists)          # moments, JohnsonFit

load.pixel.means()

# COMPARISON OF DIFFERENT PARAMETRIC MODEL FITS ACROSS ALL POWER SETTINGS

# linear spot seems less effective in white & grey images
# change in ratio of SD to MAD may give measure of change in dispersion of data.
# If SD = MAD, close to normal dist? 

# most effective models (in terms of residual MAD & largest number of 'bad' pixels identified):
# poly(2) in x  y with interaction & 2-spot

# 1-spot seems slightly less effective in white & grey channels?
#---------------------------------------------------------------------------------------------------

# compare histogram of all points vs histogram of healthy points
# histogram of all points vs (recentred) histogram of residuals after model fitting


####################################################################################################
# One-step model fitting trial                                                                  ####

# 

df <- data.frame(x = c(1:1996^2) %/% 1996 + 1,
                 y = c(1:1996^2) %% 1996,
                 z = sqrt(((c(1:1996^2) %/% 1996 + 1) - 1023.5)^2 + ((c(1:1996^2) %% 1996) - 992.5)^2),
                 p = factor(join.panels(array(sort(rep(c(1:32), 128*1024)), dim = c(128, 1024, 32)))),
                 p.x = c(join.panels(array(rep(rep(c(1:128), 1024), 32), dim = c(128, 1024, 32)))),
                 p.y = c(join.panels(array(rep(sort(rep(c(1:1024), 128)), 32), dim = c(128, 1024, 32)))),
                 value = c(pw.m[,,"white", "150828"]))

lm.test <- lm(value ~ z + p + p.x + p.y, data = df)
summary(lm.test)

pixel.image(matrix(lm.test$fitted.values, ncol = 1996))
pixel.image(matrix(lm.test$residuals, ncol = 1996))

# further panelwise fitting clearly needed. Keep stepwise fitting approach - more control.


####################################################################################################

# 1 - BASE MODEL                                                                                ####

# 'base' model for comparison: quadratic circular model, linear panel fit

circ.lm <- spot.lm(pw.m[ , , "white", "150828"], order = 2)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)

pixel.image(circ.res)

panel.lm <- panel.lm(circ.res, "x + y")
panel.res <- circ.res - panel.lm$fitted.values

pixel.image(panel.res, break.levels = sd.levels(pw.m[ , , "white", "150828"] - mean(pw.m[ , , "white", "150828"])))

# extra cross-panel effect fitted: minimal change to residuals
{
    y.lm <- panel.lm(panel.res, "poly(y, 2)")
    y.res <- panel.res - y.lm$fitted.values
    
    pixel.image(y.res)
    
    s.hist((pw.m[ , , "white", "150828"] - mean(pw.m[ , , "white", "150828"])), 
           col = "black", main = "Residuals at each stage of model fitting")
    
    hist(circ.res, add = T, breaks = "fd", col = adjustcolor("blue", alpha = 0.5), border = adjustcolor("blue", alpha = 0.3))
    hist(panel.res, add = T, breaks = "fd", col = adjustcolor("green3", alpha = 0.5), border = adjustcolor("green3", alpha = 0.3))
    hist(y.res, add = T, breaks = "fd", col = adjustcolor("gold", alpha = 0.3), border = adjustcolor("gold", alpha = 0.1))
    legend("topleft", pch = 15, col = c("black", "blue", "green3", "gold"), legend = c("Raw data (mean translated to 0)", "Circular only", "Panels added", "Panel.Y added"), bty = "n")
    
}

####################################################################################################

# 2 - ROBUST REGRESSION                                                                         ####
#---------------------------------------------------------------------------------------------------
# doesn't seem to make an enormous difference here, but maybe compare over old (more damaged) data
#---------------------------------------------------------------------------------------------------
# same model fitted, but using robust regression
circ.lm.rob <- spot.lm(pw.m[ , , "white", "150828"], order = 2, robust = T)
circ.res.rob <- matrix(circ.lm.rob$residuals, ncol = 1996)

pixel.image(circ.res.rob)

panel.lm.rob <- panel.lm(circ.res, "x + y")(circ.res.rob, "x + y", robust = T)
panel.res.rob <- circ.res.rob - panel.lm.rob$fitted.values

pixel.image(panel.res.rob, break.levels = sd.levels(pw.m[ , , "white", "150828"] - mean(pw.m[ , , "white", "150828"])))

# histogram comparison
{
    s.hist((pw.m[ , , "white", "150828"] - mean(pw.m[ , , "white", "150828"])), 
           col = "black", main = "Residuals at each stage of model fitting")
    hist(circ.res, add = T, breaks = "fd", col = adjustcolor("blue", alpha = 0.4), border = adjustcolor("blue", alpha = 0.3))
    hist(circ.res.rob, add = T, breaks = "fd", col = adjustcolor("slateblue1", alpha = 0.4), border = adjustcolor("slateblue1", alpha = 0.3))
    
    hist(panel.res, add = T, breaks = "fd", col = adjustcolor("green3", alpha = 0.4), border = adjustcolor("green3", alpha = 0.3))
    hist(panel.res.rob, add = T, breaks = "fd", col = adjustcolor("gold", alpha = 0.4), border = adjustcolor("gold", alpha = 0.3))
    legend("topleft", pch = 15, bty = "n", 
           col = c("black", "blue", "slateblue1", "green3", "gold"),
           legend = c("Raw data (mean translated to 0)", "Circular only", "Robust circular",
                      "Panels added", "Robust panels"))
}


bp <- bp.parametric.quantiles(150828, "white", spot.order = 2, panel.terms = "x + y", details = T)


####################################################################################################

# 3 - COMPARE MODEL EFFECTIVENESS BY QUANTILES                                                  ####

model.summ <- read.csv("./Other-data/Model-fit-results.csv", as.is = T)

# Most stringent models (& lowest residual SD) are those with most parameters: 
# 2-spot, 2-poly panels with interaction.

# identify bad pixels in each image over two dates, using this model
# choose two adjacent dates so that spots in screen will be as similar as possible
# 150108 and 150113 are the two closest batches, so compare these
{
    spot.b.150108 <- spot.lm(pw.m[ , , "black", "150108"], order = 2)
    spot.g.150108 <- spot.lm(pw.m[ , , "grey", "150108"], order = 2)
    spot.w.150108 <- spot.lm(pw.m[ , , "white", "150108"], order = 2)
    
    # plot residuals
    {
        pixel.image(matrix(spot.b.150108$residuals, ncol = 1996), title = "Residuals after circular spot - black")
        pixel.image(matrix(spot.g.150108$residuals, ncol = 1996), title = "Residuals after circular spot - grey")
        pixel.image(matrix(spot.w.150108$residuals, ncol = 1996), title = "Residuals after circular spot - white")
    }

    panel.b.150108 <- panel.lm(matrix(spot.b.150108$residuals, ncol = 1996), "poly(x, 2) * poly(y, 2)")
    panel.g.150108 <- panel.lm(matrix(spot.g.150108$residuals, ncol = 1996), "poly(x, 2) * poly(y, 2)")
    panel.w.150108 <- panel.lm(matrix(spot.w.150108$residuals, ncol = 1996), "poly(x, 2) * poly(y, 2)")
    
    res.b.150108 <- matrix(spot.b.150108$residuals, ncol = 1996) - panel.b.150108$fitted.values
    res.g.150108 <- matrix(spot.g.150108$residuals, ncol = 1996) - panel.g.150108$fitted.values
    res.w.150108 <- matrix(spot.w.150108$residuals, ncol = 1996) - panel.w.150108$fitted.values
    
    # plot panels & new residuals
    {
        pixel.image(panel.b.150108$fitted.values, title = "Panels fitted - black")
        pixel.image(panel.g.150108$fitted.values, title = "Panels fitted - grey")
        pixel.image(panel.w.150108$fitted.values, title = "Panels fitted - white")
        
        pixel.image(res.b.150108, title = "Residuals - black")
        pixel.image(res.g.150108, title = "Residuals - grey")
        pixel.image(res.w.150108, title = "Residuals - white")

        s.hist(spot.b.150108$residuals, col = "black", main = "Residuals before and after panel fitting: black")
        hist(res.b.150108, add = T, breaks = "fd", col = adjustcolor("blue", alpha = 0.3),
             border = adjustcolor("blue", alpha = 0.3))
        
        s.hist(spot.g.150108$residuals, col = "black", main = "Residuals before and after panel fitting: grey")
        hist(res.g.150108, add = T, breaks = "fd", col = adjustcolor("blue", alpha = 0.3),
             border = adjustcolor("blue", alpha = 0.3))
        
        s.hist(spot.w.150108$residuals, col = "black", main = "Residuals before and after panel fitting: white")
        hist(res.w.150108, add = T, breaks = "fd", col = adjustcolor("blue", alpha = 0.3),
             border = adjustcolor("blue", alpha = 0.3))
    }
    
    # unconvinced by fit of this model to white image, although errors in black & grey seem less systematic.
    spot.tmp <- spot.lm(pw.m[ , , "white", "150828"], order = 2)
    panel.tmp <- panel.lm(matrix(spot.tmp$residuals, ncol = 1996), "poly(x, 2) * poly(y, 2)")
    res.tmp <- matrix(spot.tmp$residuals, ncol = 1996) - panel.tmp$fitted.values
    


    }

# compare each date's bad pixels to system bad pixel map and each other
{
    
}

# run same image with robust regression, compare to previous model's bad pixel map
{
    
}

bp.parametric.quantiles(150828, "black", spot.order = NA, panel.terms = "poly(x, 2) * poly(y, 2)"),
                   bp.parametric.quantiles(150828, "black", spot.order = 1, panel.terms = "poly(x, 2) * poly(y, 2)"),
                   bp.parametric.quantiles(150828, "black", spot.order = 2, panel.terms = "poly(x, 2) * poly(y, 2)"))

# 3.5 - EFFECT OF BAD CORNER ON PANEL FIT                                                       ####

# top-left corner panel is always odd - possibly because model is skewed by 
cp <- subpanels(pw.m[,,"white", "150828"])[3:128,1:1004,1]
pixel.image(cp)

terms <- "X1 + X2";     form <- as.formula(paste0("value ~ ", terms))

{
    cp.lm <- lm(form, data = melt(cp))
    cp.rlm <- rlm(form, data = melt(cp))
    
    par(mfrow = c(1,2))
    pixel.image(matrix(cp.lm$fitted.values, nrow = 126), title = terms)
    pixel.image(matrix(cp.rlm$fitted.values, nrow = 126), title = "robust")
    par(mfrow = c(1,1))
}   # run model & plot

{
    cp.lm.2 <- lm(form, data = melt(cp[,-1]))
    cp.rlm.2 <- rlm(form, data = melt(cp[,-1]))
    
    par(mfrow = c(1,2))
    pixel.image(matrix(cp.lm.2$fitted.values, nrow = 126), title = paste0("Col1 removed: ",terms))
    pixel.image(matrix(cp.rlm.2$fitted.values, nrow = 126), title = "robust")
    par(mfrow = c(1,1))
}   # now without column 1


# 4 - JOHNSON DIST VS QUANTILES                                                                 ####

# fit a Johnson distribution and plot residuals (esp. in tails).
# Is there a cutoff point at which point are consistently higher than expected?

# find Johnson quantiles that match to quantile cutpoints: consistent?