
# EXPECT TO FIND CIRCULAR SPOT & LINEAR PANEL GRADIENT.
# WHAT REMAINS AFTER FITTING THIS MODEL?
# HOW WELL ARE BAD PIXELS IDENTIFIED USING THIS APPROACH?

library("IO.Pixels")
par(pch = 20)

load.pixel.means()
bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)

# use image date for which we have 'official' bad pixel map
dt <- "160314"    

####################################################################################################

# FIT SIMPLE MODEL                                                                              ####

# white model
{
    # quadratic spot
    w.spot <- spot.lm(pw.m[,,"white", dt], o = 2, robust = T)
    w.spot.res <- matrix(w.spot$residuals, ncol = 1996)
    
    # linear panels
    w.panel <- panel.lm(w.spot.res, "x + y", robust = T)
    
        # per-panel coefficients
        plot(w.panel$models[,"x"], w.panel$models[,"y"], pch = 20, col = c(rep("black", 16), rep("red", 16)))
    
        # fitted panels
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-white-panels.pdf"))
            pixel.image(w.panel$fitted.values, panels = T, title = paste0("Panels fitted to white image after spot removed: ", dt))
        dev.off()
    
    w.res <- w.spot.res - w.panel$fitted.values
    
        # residuals
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-white-residuals.pdf"))
            pixel.image(w.res, title = paste0("Residuals after fitting spot and linear panels to white image: ", dt), panels = T)
            s.hist(w.res, col = "black", main = "Histogram of residuals", xlab = "Residual value")
        dev.off()
        
    # how well does a Johnson distribution fit the data?
        w.JF <- JohnsonFit(w.res, moment = "quant")
        w.JF.cropped <- JohnsonFit(w.res[51:1946, 51:1946], moment = "quant")
        
        
    # QQ plot shows left (negative) skew - left tail is heavier than that of Johnson dist
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-white-residual-QQ.pdf"))
        plot(qJohnson((c(1:999)/1000), w.JF), quantile(w.res,(c(1:999)/1000)), 
             pch = 20, ylab = "Observed quantile", xlab = "Johnson quantile",
             main = paste0("Johnson Q-Q plot: white, ", dt), col = adjustcolor("blue", alpha = 0.5), asp = T)
        abline(0,1,col = "red")
        abline(h = quantile(w.res, c(0.01, 0.99)), col = "skyblue", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), w.JF), col = "skyblue", lty = 2)
        
        points(qJohnson((c(1:999)/1000), w.JF.cropped), quantile(w.res[51:1946, 51:1946],(c(1:999)/1000)), 
               pch = 20, col = adjustcolor("orange", alpha = 0.5))
        abline(h = quantile(w.res[51:1946, 51:1946], c(0.01, 0.99)), col = "gold", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), w.JF.cropped), col = "gold", lty = 2)

        legend("bottomright", legend = c("Quantiles of all data", "Quantiles of cropped data",
                                     "Q.01 and Q.99 of all data", "Q.01 and Q.99 of cropped data"),
               pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 2), 
               col = c(adjustcolor(c("blue", "orange"), alpha = 0.5), "skyblue", "gold"))
        dev.off()
        }
}

# grey model
{
    # quadratic spot
    g.spot <- spot.lm(pw.m[,,"grey", dt], o = 2, robust = T)
    g.spot.res <- matrix(g.spot$residuals, ncol = 1996)
    
    # linear panels
    g.panel <- panel.lm(g.spot.res, "x + y", robust = T)
    
        # per-panel coefficients
        plot(g.panel$models[,"x"], g.panel$models[,"y"], pch = 20, col = c(rep("black", 16), rep("red", 16)))
    
        # fitted panels
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-grey-panels.pdf"))
            pixel.image(g.panel$fitted.values, panels = T, title = paste0("Panels fitted to grey image after spot removed: ", dt))
        dev.off()
    
    g.res <- g.spot.res - g.panel$fitted.values
    
    # residuals
    pdf(paste0("./Plots/Simple-model-fitting/", dt, "-grey-residuals.pdf"))
        pixel.image(g.res, title = paste0("Residuals after fitting spot and linear panels to grey image: ", dt), panels = T)
        s.hist(g.res, col = "black", main = "Histogram of residuals", xlab = "Residual value")
    dev.off()
    
    # how well does a Johnson distribution fit the data?
    g.JF <- JohnsonFit(g.res, moment = "quant")
    g.JF.cropped <- JohnsonFit(g.res[51:1946, 51:1946], moment = "quant")
    
    
    # QQ plot shows heavy tails on both ends of distribution
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-grey-residual-QQ.pdf"))
        plot(qJohnson((c(1:999)/1000), g.JF), quantile(g.res,(c(1:999)/1000)), 
             pch = 20, ylab = "Observed quantile", xlab = "Johnson quantile",
             main = paste0("Johnson Q-Q plot: grey, ", dt), col = adjustcolor("blue", alpha = 0.5), asp = T)
        abline(0,1,col = "red")
        abline(h = quantile(g.res, c(0.01, 0.99)), col = "skyblue", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), g.JF), col = "skyblue", lty = 2)
        
        points(qJohnson((c(1:999)/1000), g.JF.cropped), quantile(g.res[51:1946, 51:1946],(c(1:999)/1000)), 
               pch = 20, col = adjustcolor("orange", alpha = 0.5))
        abline(h = quantile(g.res[51:1946, 51:1946], c(0.01, 0.99)), col = "gold", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), g.JF.cropped), col = "gold", lty = 2)
        
        legend("bottomright", legend = c("Quantiles of all data", "Quantiles of cropped data",
                                         "Q.01 and Q.99 of all data", "Q.01 and Q.99 of cropped data"),
               pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 2), 
               col = c(adjustcolor(c("blue", "orange"), alpha = 0.5), "skyblue", "gold"))
        dev.off()
    }
}

# black model
{
    # quadratic spot
    b.spot <- spot.lm(pw.m[,,"black", dt], o = 2, robust = T)
    b.spot.res <- matrix(b.spot$residuals, ncol = 1996)
    
    # linear panels
    b.panel <- panel.lm(b.spot.res, "x + y", robust = T)
    
        # per-panel coefficients
        plot(b.panel$models[,"x"], b.panel$models[,"y"], pch = 20, col = c(rep("black", 16), rep("red", 16)))
    
        # fitted panels
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-black-panels.pdf"))
            pixel.image(b.panel$fitted.values, panels = T, title = paste0("Panels fitted to black image after spot removed: ", dt))
        dev.off()
    
    b.res <- b.spot.res - b.panel$fitted.values
    
    # residuals
    pdf(paste0("./Plots/Simple-model-fitting/", dt, "-black-residuals.pdf"))
        pixel.image(b.res, title = paste0("Residuals after fitting spot and linear panels to black image: ", dt), panels = T)
        s.hist(b.res, col = "black", main = "Histogram of residuals", xlab = "Residual value")
    dev.off()
    
    # how well does a Johnson distribution fit the data?
    b.JF <- JohnsonFit(b.res, moment = "quant")
    b.JF.cropped <- JohnsonFit(b.res[51:1946, 51:1946], moment = "quant")
    
    
    # QQ plot shows heavy tails on both ends of distribution
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-black-residual-QQ.pdf"))
        plot(qJohnson((c(1:999)/1000), b.JF), quantile(b.res,(c(1:999)/1000)), 
             pch = 20, ylab = "Observed quantile", xlab = "Johnson quantile",
             main = paste0("Johnson Q-Q plot: black, ", dt), col = adjustcolor("blue", alpha = 0.5), asp = T)
        abline(0,1,col = "red")
        abline(h = quantile(b.res, c(0.01, 0.99)), col = "skyblue", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), b.JF), col = "skyblue", lty = 2)
        
        points(qJohnson((c(1:999)/1000), b.JF.cropped), quantile(b.res[51:1946, 51:1946],(c(1:999)/1000)), 
               pch = 20, col = adjustcolor("orange", alpha = 0.5))
        abline(h = quantile(b.res[51:1946, 51:1946], c(0.01, 0.99)), col = "gold", lty = 2)
        abline(v = qJohnson(c(0.01, 0.99), b.JF.cropped), col = "gold", lty = 2)
        
        legend("topleft", legend = c("Quantiles of all data", "Quantiles of cropped data",
                                         "Q.01 and Q.99 of all data", "Q.01 and Q.99 of cropped data"),
               pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 2), 
               col = c(adjustcolor(c("blue", "orange"), alpha = 0.5), "skyblue", "gold"))
        dev.off()
    }
}

####################################################################################################

# BAD PIXELS BY QUANTILE                                                                        ####

# use same threshold for all images - 0.1% + 1.5*IQR
{bp <- rbind(reset.bp(dt),
             data.frame(which(w.res > (quantile(w.res, 0.999) + (1.5*IQR(w.res))), arr.ind = T), src = "white", type = "bright"),
             data.frame(which(w.res < (quantile(w.res, 0.001) - (1.5*IQR(w.res))), arr.ind = T), src = "white", type = "dim"),
             data.frame(which(g.res > (quantile(g.res, 0.999) + (1.5*IQR(g.res))), arr.ind = T), src = "grey", type = "bright"),
             data.frame(which(g.res < (quantile(g.res, 0.001) - (1.5*IQR(g.res))), arr.ind = T), src = "grey", type = "dim"),
             data.frame(which(b.res > (quantile(b.res, 0.999) + (1.5*IQR(b.res))), arr.ind = T), src = "black", type = "bright"),
             data.frame(which(b.res < (quantile(b.res, 0.001) - (1.5*IQR(b.res))), arr.ind = T), src = "black", type = "dim"))

write.csv(bp, paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.csv"), row.names = F)
}

bp <- read.csv(paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.csv"))
table(bp$src, bp$type)

####################################################################################################

# POSSIBLE IMPROVEMENTS ETC                                                                     ####

# visible line of dead pixels with blooming in black residual image for 160314, panel U4
{
    pixel.image(subpanels(b.res)[,150:300, "U4"], break.levels = sd.levels(b.res), title = "Dead pixels - black residuals")
    pixel.image(subpanels(g.res)[,150:300, "U4"], break.levels = sd.levels(g.res), title = "Dead pixels - grey residuals")
    pixel.image(subpanels(w.res)[,150:300, "U4"], break.levels = sd.levels(w.res), title = "Dead pixels - white residuals")

    pixel.image(subpanels(pw.m[,,"black", dt])[,150:300, "U4"], title = "Dead pixels - black image")
    pixel.image(subpanels(pw.m[,,"grey", dt])[,150:300, "U4"], title = "Dead pixels - grey image")
    pixel.image(subpanels(pw.m[,,"white", dt])[,150:300, "U4"], title = "Dead pixels - white image")
}

# fit an 'edge curve' to account for dropoff in last 100-200px (single curve for upper, single for lower)
{
    lower <- rlm(value ~ poly(X2, 4), melt(w.res[,1:200]))
    upper <- rlm(value ~ poly(X2, 4), melt(w.res[,1796:1996]))

    o.plot(matrix(lower$fitted.values, nrow = 1996)[1,])
    o.plot(matrix(upper$fitted.values, nrow = 1996)[1,], add = T, col = "blue")

    w.edged <- w.res
    w.edged[,1:200] <- w.res[,1:200] - matrix(lower$fitted.values, nrow = 1996)
    w.edged[,1796:1996] <- w.res[,1796:1996] - matrix(upper$fitted.values, nrow = 1996)
    
    pixel.image(w.edged)
    # looks suspiciously over-smoothed at the edges now - less variable border
    s.hist(w.res)
    hist(w.edged, add = T, breaks = "fd", col = adjustcolor("cornflowerblue", alpha = 0.5),
         border = adjustcolor("cornflowerblue", alpha = 0.5))
    # lower tail definitely reduced, upper barely affected 
}


    