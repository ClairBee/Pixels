
# EXPECT TO FIND CIRCULAR SPOT & LINEAR PANEL GRADIENT.
# WHAT REMAINS AFTER FITTING THIS MODEL?
# HOW WELL ARE BAD PIXELS IDENTIFIED USING THIS APPROACH?

library("IO.Pixels")
library(scatterplot3d)      # for 3d scatterplots (eg. for SVM)

load.pixel.means()

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
        par(mfrow = c(1,2))
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
    par(mfrow = c(1,2))
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
    par(mfrow = c(1,2))
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
    bp <- rbind(reset.bp(dt),
                 data.frame(which(w.res > (quantile(w.res, 0.999) + (1.5*IQR(w.res))), arr.ind = T), src = "white", type = "bright"),
                 data.frame(which(w.res < (quantile(w.res, 0.001) - (1.5*IQR(w.res))), arr.ind = T), src = "white", type = "dim"),
                 data.frame(which(g.res > (quantile(g.res, 0.999) + (1.5*IQR(g.res))), arr.ind = T), src = "grey", type = "bright"),
                 data.frame(which(g.res < (quantile(g.res, 0.001) - (1.5*IQR(g.res))), arr.ind = T), src = "grey", type = "dim"),
                 data.frame(which(b.res > (quantile(b.res, 0.999) + (1.5*IQR(b.res))), arr.ind = T), src = "black", type = "bright"),
                 data.frame(which(b.res < (quantile(b.res, 0.001) - (1.5*IQR(b.res))), arr.ind = T), src = "black", type = "dim"))

    write.csv(bp, paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.csv"), row.names = F)


# combine bad pixel maps & compare
    bp <- read.csv(paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.csv"))
    table(bp$src, bp$type)

    # import 'official' bad pixel map
    bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)

    # combine into map of all bad pixels identified by both methods
    bp.all <- combine.maps(bp, bpm)
    bp.all$type <- ordered(bp.all$type, levels = c("hot", "dead", "bright", "dim"))
    
    saveRDS(bp.all, paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.rds"))
    
    bp.healthy <- sample.healthy(bp.all)

    # plot bad vs healthy pixels
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt,"-bad-vs-healthy.pdf"), width = 12, height = 4)
        par(mfrow = c(1,3), pch = 20, mar = c(2,2,3,1))
            plot(pw.m[,,"black",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - black, ", dt))
            points(pw.m[,,"black",dt][as.matrix(bp.all[,1:2])], col = adjustcolor("blue", alpha = 0.2))
            
            plot(pw.m[,,"grey",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - grey, ", dt))
            points(pw.m[,,"grey",dt][as.matrix(bp.all[,1:2])], col = adjustcolor("green3", alpha = 0.2))
            
            plot(pw.m[,,"white",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - white, ", dt))
            points(pw.m[,,"white",dt][as.matrix(bp.all[,1:2])], col = adjustcolor("turquoise3", alpha = 0.2))
        dev.off()
    }
    
    # plot bad pixels by means of identification
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt,"-bad-px-by-source.pdf"), width = 12, height = 4)
        par(mfrow = c(1,3), pch = 20, mar = c(2,2,3,1))
        
            plot(pw.m[,,"black",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - black, ", dt))
            points(pw.m[,,"black",dt][as.matrix(bp.all[,1:2])],
                   col = adjustcolor(c("slateblue1", "chartreuse3", "red")[bp.all$map], alpha = 0.2))
            legend("topleft", pch = 20, col = c("slateblue1", "chartreuse3", "red", "gold"), bty = "n",
                   legend = c("CB map only", "Both maps", "Auto map only", "Healthy sample"))
        
            plot(pw.m[,,"grey",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - grey, ", dt))
            points(pw.m[,,"grey",dt][as.matrix(bp.all[,1:2])],
                   col = adjustcolor(c("slateblue1", "chartreuse3", "red")[bp.all$map], alpha = 0.2))
            legend("topleft", pch = 20, col = c("slateblue1", "chartreuse3", "red", "gold"), bty = "n",
                   legend = c("CB map only", "Both maps", "Auto map only", "Healthy sample"))
            
            plot(pw.m[,,"white",dt][bp.healthy], ylim = c(0,65535), col = adjustcolor("gold", alpha = 0.2),
                 ylab = "Pixel value", main = paste0("Bad pixels - white, ", dt))
            points(pw.m[,,"white",dt][as.matrix(bp.all[,1:2])],
                   col = adjustcolor(c("slateblue1", "chartreuse3", "red")[bp.all$map], alpha = 0.2))
            legend("topleft", pch = 20, col = c("slateblue1", "chartreuse3", "red", "gold"), bty = "n",
                   legend = c("CB map only", "Both maps", "Auto map only", "Healthy sample"))
        dev.off()
    }
    
    # check line root: zoom to xlim = c(420,440), ylim = c(1180,1220),
    plot(bp.all[bp.all$map != "old",1:2], asp = T, pch = 20, 
         col = c("red", "blue", "gold", "green3")[bp.all$type[bp.all$map != "old"]])
    points(bp.all[bp.all$map != "new",1:2], pch = 1)

####################################################################################################

    
    
    
    
# DEVELOPMENT OF BAD PIXELS OVER TIME                                                           ####
    bp <- readRDS(paste0("./Plots/Simple-model-fitting/", dt, "-bpm-IQR.rds"))
    gp <- sample.healthy(bp)
    dimnames(gp)[[2]] <- c("row", "col")
    px <- as.matrix(rbind(bp[,1:2], gp))
    

    # plot development of bad pixels
    focus <- "hot"
    {
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-", focus, "-px.pdf"), width = 15, height = 5)
        px <- as.matrix(rbind(bp[bp$type == focus,1:2], gp[1:nrow(bp[bp$type == focus,1:2]),]))
        cols <- c("green4", "cyan3", "dodgerblue3", "blue2", "purple", "magenta3", "red",
                  "orange", "gold", "green", "black")
        par(mfrow = c(1,3))
        for (j in 1:3) {
            plot(0, type = "n", ylim = c(0, 65535), ylab = "Pixel value", xlim = c(0, nrow(px)),
                 main = paste0("Progression of ", focus, " pixels - ", dimnames(pw.m)[[3]][j]," images"))
            abline(v = nrow(bp[bp$type == focus,1:2]), col = "red")
            for (i in 1:11) {
                points(pw.m[,,j, i][px], pch = 20,
                       col = adjustcolor(cols[i], alpha = 0.2))
            }
        }
        par(mfrow = c(1,1))
        dev.off()
    }

    # transects over time
    {
        px <- as.matrix(rbind(bp[,1:2], gp))
        
        bp.b <- apply(pw.m[,,"black",], 3, "[", as.matrix(bp[,1:2]))
        bp.g <- apply(pw.m[,,"grey",], 3, "[", as.matrix(bp[,1:2])
        bp.w <- apply(pw.m[,,"white",], 3, "[", as.matrix(bp[,1:2])
                      
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-badpx-black.pdf"), width = 15, height = 4)
        par(mfrow = c(1,5))
        {
            matplot(t(apply(pw.m[,,"black",], 3, "[", as.matrix(gp))), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Healthy pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.b[bp$type == "dead",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dead pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.b[bp$type == "dim",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dim pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)

            matplot(t(bp.b[bp$type == "bright",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Bright pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.b[bp$type == "hot",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Hot pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            

        }
        dev.off()
        
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-badpx-grey.pdf"), width = 15, height = 4)
        par(mfrow = c(1,5))
        {
            matplot(t(apply(pw.m[,,"grey",], 3, "[", as.matrix(gp))), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Healthy pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.g[bp$type == "dead",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dead pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.g[bp$type == "dim",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dim pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.g[bp$type == "bright",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Bright pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.g[bp$type == "hot",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Hot pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
        }
        dev.off()
        
        pdf(paste0("./Plots/Simple-model-fitting/", dt, "-badpx-white.pdf"), width = 15, height = 4)
        par(mfrow = c(1,5))
        {
            matplot(t(apply(pw.m[,,"white",], 3, "[", as.matrix(gp))), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Healthy pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.w[bp$type == "dead",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dead pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.w[bp$type == "dim",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Dim pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.w[bp$type == "bright",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Bright pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
            
            matplot(t(bp.w[bp$type == "hot",]), ylim = c(0,65535), ylab = "",
                    type = "l", main = "Hot pixels", 
                    xaxt = "none")
            axis(1, at = c(1:11), labels = dimnames(pw.m)[[4]], las = 2)
        }
        dev.off()
    }

    
    
 

####################################################################################################
    
    
# MISCELLANEOUS                                                                                 ####

# visible line of dead pixels with blooming in black residual image for 160314, panel U4
# 'root' cluster at [44:46, 206:208, "U4"] == [426:428, 1198:1200]
    {
        pixel.image(subpanels(b.res)[,150:300, "U4"], break.levels = sd.levels(b.res), title = "Dead pixels - black residuals")
        pixel.image(subpanels(g.res)[,150:300, "U4"], break.levels = sd.levels(g.res), title = "Dead pixels - grey residuals")
        pixel.image(subpanels(w.res)[,150:300, "U4"], break.levels = sd.levels(w.res), title = "Dead pixels - white residuals")

        pixel.image(subpanels(pw.m[,,"black", dt])[,150:300, "U4"], title = "Dead pixels - black image")
        pixel.image(subpanels(pw.m[,,"grey", dt])[,150:300, "U4"], title = "Dead pixels - grey image")
        pixel.image(subpanels(pw.m[,,"white", dt])[,150:300, "U4"], title = "Dead pixels - white image")
    
        check.cluster <- function(dt, box = F) {
            dt <- toString(dt)
            par(mfrow = c(1,3))
            
            pixel.image(subpanels(pw.m[,,"black", dt])[, 150:300, "U4"], title = "Dead pixels - black image")
            if (box) {
                lines(matrix(c(30, 45, 30, 72, 60, 72, 60, 45, 30, 45), ncol = 2, byrow = T), lty = 2)
            }
            
            pixel.image(subpanels(pw.m[,,"grey", dt])[, 150:300, "U4"], title = "Dead pixels - grey image")
            if (box) {
                lines(matrix(c(30, 45, 30, 72, 60, 72, 60, 45, 30, 45), ncol = 2, byrow = T), lty = 2)
            }
            
            pixel.image(subpanels(pw.m[,,"white", dt])[, 150:300, "U4"], title = "Dead pixels - white image")
            if (box) {
                lines(matrix(c(30, 45, 30, 72, 60, 72, 60, 45, 30, 45), ncol = 2, byrow = T), lty = 2)
            }
            
            par(mfrow = c(1,1))
        }
        
        # root seems to have existed since the first images, but line only appeared in latest set.
        check.cluster(141009)
        check.cluster(150730)
        check.cluster(150828)
        check.cluster(151015)
        
        o.plot(pw.m[426, 993:1996,"black", dt])
        o.plot(pw.m[427, 993:1996,"black", dt], add = T, col = adjustcolor("red", alpha = 0.2))
        o.plot(pw.m[428, 993:1996,"black", dt], add = T, col = adjustcolor("blue", alpha = 0.2))
    }
    
# line isn't picked up - not (yet) sufficiently far from bulk of image.
# damaged line is currently around 250 higher than it should be
    {
        o.plot(pw.m[426, 993:1996,"black", dt])
        o.plot(pw.m[427, 993:1996,"black", dt], add = T, col = adjustcolor("red", alpha = 0.2))
        o.plot(pw.m[428, 993:1996,"black", dt], add = T, col = adjustcolor("blue", alpha = 0.2))
        
        o.plot(pw.m[427, 993:1996,"black", dt] - pw.m[426, 993:1996,"black", dt], ylim = c(-1000,1000))
            abline(h = median((pw.m[427, 993:1996,"black", dt] - pw.m[426, 993:1996,"black", dt])[203:1004]), col = "blue")
        o.plot(pw.m[428, 993:1996,"black", dt] - pw.m[427, 993:1996,"black", dt], add = T, col = adjustcolor("red", alpha = 0.3))
            abline(h = median((pw.m[428, 993:1996,"black", dt] - pw.m[427, 993:1996,"black", dt])[203:1004]), col = "red")

        sp <- subpanels(pw.m[,,"black", dt])[,1:1004,"U4"]
        diff <- sp[1:127,] - sp[2:128,]
        pixel.image(diff)
        s.hist(diff)
        o.plot(diff[,566]); abline(v = 44.5, col = "red")
        med.diff <- apply(diff, 1, median)
        o.plot(med.diff, title = "med. diff"); abline(v = 44.5, col = "red")
    }

# flat field (offset/gain) correction on damaged area: removes line effect
    {
        R <- pw.m[326:528, 1098:1300, "grey", dt]
        D <- pw.m[326:528, 1098:1300, "black", dt]
        FF <- pw.m[326:528, 1098:1300, "white", dt]
        m <- mean(FF - D)
        
        corr <- m * (R - D) / (FF - D)
        corr[is.na(corr)] <- 0      # otherwise get NA where FF == D
        
        pixel.image(D)
        s.hist(corr)
    }
    
####################################################################################################
# clustering
    # hierarchical clustering of all bad pixels (+ sample of healthy ones)
    {
        px <- as.matrix(rbind(bp[,1:2], gp))
        #---------------------------------------------------------------------------------
        
        bp.vals <- apply(pw.m[,,,dt], 3, "[", px)
        
        dist <- dist(bp.vals, method = "euclidean")
        clust <- hclust(dist)
        
        plot(clust, labels = F)
        abline(h = 40000, col = "red")
        clust.cat <- cutree(clust, h = 40000)
        
        px.match <- cbind(rbind(bp[,1:3], cbind(gp, type = "healthy")), clust = clust.cat)
        table(px.match$type, px.match$clust)
        #---------------------------------------------------------------------------------
        
    }
    
    # could SVM work? Try plotting first 3 images in 3d, or all batches in 3d
    {    
        for (i in 1:dim(pw.m)[4]) {
            nm <- dimnames(pw.m)[[4]][i]
            nm <- paste0(substring(nm,1,2), "-", substring(nm,3,4), "-", substring(nm,5,6))
            scatterplot3d(apply(pw.m[,,,i], 3, "[", px),
                          pch = 20, xlim = c(0,65535), ylim = c(0,65535), zlim = c(0,65535),
                          color = adjustcolor(c(rep("red", nrow(bp)), rep("gold", nrow(bp))), alpha = 0.2),
                          main = paste0("Bad (red) vs healthy points: ", nm))
        }
    }
    # only most extreme points are linearly separable - probably not subtle enough
    
####################################################################################################
    
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

####################################################################################################

# would a higher-order circular regression be more useful?
    {# spot profiles
        plot(c(0:1433), predict(spot.lm(pw.m[,,"black", dt], o = 1, robust = T),
                                data.frame(z = (c(0:1433)))), type = "l", xlab = "z", ylab = "",
             main = "Spot profiles (black)")
        
        for (i in 2:6) {
            points(c(0:1433), predict(spot.lm(pw.m[,,"black", dt], o = i, robust = T),
                                      data.frame(z = (c(0:1433)))), type = "l",  
                   col = c("blue", "red", "green4", "gold", "magenta3")[i-1])
        }
        
        # order 6 actually looks like it might be the best...
        b.spot <- spot.lm(pw.m[,,"black", dt], o = 2, robust = T)
        b.spot.res <- matrix(b.spot$residuals, ncol = 1996)
        
        # linear panels
        b.panel <- panel.lm(b.spot.res, "x + y", robust = T)
        b.res <- b.spot.res - b.panel$fitted.values
        
        pixel.image(b.res, title = paste0("Residuals after fitting 6-spot and linear panels to black image: ", dt), panels = T)
        s.hist(b.res, col = "black", main = "Histogram of residuals", xlab = "Residual value")
        mad(b.res)          # 98.93546 (101.2448 with 2-spot)
    } 
    # didn't give any improvement in residual MAD
    