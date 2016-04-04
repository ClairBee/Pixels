
library("IO.Pixels")
library(SuppDists)          # needed for JohnsonFit, dJohnson
load.pixel.means()

im <- pw.m[,,"white", "150828"]

# very little change over first few iterations (circular spot is perhaps minutely more accurate)
# better to concentrate on improving model design instead.

####################################################################################################
# COMPARE MODEL AT EACH POWER SETTING                                                           ####

circ.o <- 2; panel.o <- 1

#---------------------------------------------------------------------------------------------------

circ.b <- fit.circular.lm.poly(bp.b[,,"black", "150828"], o = 2)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)

panel.lm <- fit.panel.lm(circ.res)
panel.res <- circ.res - panel.lm$fitted.values


###################################################################################################
# FITTING WITHOUT REMOVING EXTREME VALUES                                                      ####

# fit circular model
circ.lm <- fit.circular.lm.poly(im, 2)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)
mad(circ.res); mean(circ.res); sd(circ.res)     # 545, 0, 706

# inspect plots: further fitting obviously needed? 
{
    par(mfrow = c(1,2)
        s.hist(circ.res, prob = T)
        lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(circ.res), sd = sd(circ.res)), lwd = 3, col = "cornflowerblue")
        lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(circ.res), sd = mad(circ.res)), lwd = 3, col = "orange")
        lines(c(-1500:1500), dJohnson(c(-1500:1500), JohnsonFit(c(circ.res), moment = "quant")), lwd = 3, col = "green3")
        pixel.image(circ.res)
        par(mfrow = c(1,1))
}
# clear pattern across panels. Further fitting definitely required.

#====================================================================================
# start with linear per-panel model
panel.lm <- fit.panel.lm(circ.res)
panel.res <- circ.res - panel.lm$fitted.values
mad(panel.res); mean(panel.res); sd(panel.res)     # 275, 0, 525

# inspect plots: further fitting needed?
{
    par(mfrow = c(1,2))
    s.hist(panel.res, prob = T)
    lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(panel.res), sd = sd(panel.res)), lwd = 3, col = "cornflowerblue")
    lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(panel.res), sd = mad(panel.res)), lwd = 3, col = "orange")
    lines(c(-1500:1500), dJohnson(c(-1500:1500), JohnsonFit(c(panel.res), moment = "quant")), lwd = 3, col = "green3")
    pixel.image(panel.res)
    par(mfrow = c(1,1))
}
# still slight panelwise pattern.

# get panel-only fitted values
panels <- fit.panels(panel.lm$models)
{
    # check offsets & gradients
    pixel.image(offset)
    pixel.image(grad)
    pixel.image(offset + grad)
}


#====================================================================================
# new starting data is original data - fitted panels; fit new circular model

im2 <- im - panels$offset - panels$grad
pixel.image(im2)    # v. clear circular pattern

# fit circular model to remainder
circ.lm.2 <- fit.circular.lm.poly(im2, 2)
circ.res.2 <- matrix(circ.lm.2$residuals, ncol = 1996)
mad(circ.res.2); mean(circ.res.2); sd(circ.res.2)     # 278, 0, 523

# plot spot profile & residuals - v similar
{
    # plot profile of circular model for comparison
    plot(c(0:1433), predict(circ.lm.2, data.frame(z = (c(0:1433)))), type = "l",
         xlab = "Distance from centre of panel", ylab = "Fitted value",
         main = "Profile of fitted circular spot")
    points(c(0:1433), predict(circ.lm, data.frame(z = (c(0:1433)))), type = "l", col = "blue")
    
    
    par(mfrow = c(1,2))
    s.hist(circ.res.2, prob = T)
    lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(circ.res.2), sd = sd(circ.res.2)), lwd = 3, col = "cornflowerblue")
    lines(c(-1500:1500), dnorm(c(-1500:1500), mean = mean(circ.res.2), sd = mad(circ.res.2)), lwd = 3, col = "orange")
    lines(c(-1500:1500), dJohnson(c(-1500:1500), JohnsonFit(c(circ.res.2), moment = "quant")), lwd = 3, col = "green3")
    
    hist(panel.res, breaks = "fd", xlim = c(-2000,2000),
         col = adjustcolor("green3", alpha = 0.4), border = adjustcolor("green3", alpha = 0.3))
    hist(circ.res.2, breaks = "fd", add = T,
         col = adjustcolor("blue", alpha = 0.4), border = adjustcolor("blue", alpha = 0.3))    
    par(mfrow = c(1,1))
}

#====================================================================================
# subtract circular model from original data, fit panels again
im3 <- im - circ.lm.2$fitted.values
pixel.image(im3)

# fit panel model to remainder
panel.lm.2 <- fit.panel.lm(im3)
panel.res.2 <- im3 - panel.lm.2$fitted.values
mad(panel.res.2); mean(panel.res.2); sd(panel.res.2)     # 271, 0, 521

panels2 <- fit.panels(panel.lm.2$models)

# plot panels & residuals
{
    org.par <- par()
    par(mfrow = c(2,2), mar = c(1,1,3,1))
    pixel.image(panels$offset, title = "First LM: offsets")
    pixel.image(panels2$offset, title = "Second LM: offsets")
    
    pixel.image(panels$grad, title = "First LM: gradients")
    pixel.image(panels2$grad, title = "Second LM: gradients")
    par(mfrow = c(1,1), mar = org.par$mar)
    
    hist(panel.res, breaks = "fd", xlim = c(-2000,2000),
         col = adjustcolor("green3", alpha = 0.4), border = adjustcolor("green3", alpha = 0.3))
    hist(panel.res.2, breaks = "fd", add = T,
         col = adjustcolor("blue", alpha = 0.4), border = adjustcolor("blue", alpha = 0.3))  
}
# very little change. Points appear to be fractionally more concentrated in centre, but could easily by simple random variation.
}
###################################################################################################

###################################################################################################

