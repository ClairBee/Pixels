
# CONSIDERING POSSIBLE IMPROVEMENTS TO MODEL

library("IO.Pixels")

load.pixel.means()

########################################################################################################
# DOUBLE SPOT IN WHITE/GREY IMAGES
#   No obvious improvement (second spot basically 0, even up to o = 4)
#   Panel gradients still erratic

# FIT HYPER-ELLIPSE
# Need to include constraint that terms sum to 1. Check with Julia.

# FIT MINI-PANELS
#   Not a big improvement on MAD (no change in SD), but changes spatial distribution of bad pixels
#   Actually makes fit worse in grey & white images because of edges 

# FIT PANELS TO CROPPED DATA
#   Panel gradients may actually be more erratic using this method
#   Again, no discernible improvement

# CHECK JOHNSON DISTRIBUTION FIT FOR SD
#   lognormal is v close to data in black images, if MAD is used instead of SD
#   however Johnson still looks closer in grey & white images
#   sJohnson also gives summary of skewness etc - may be useful
#   ks.test rejects all model fits at 5% level. To be expected due to sample size/known outliers?

# JOHNSON FITTING TO SD
#   15-01-08 black images 1:10 and 11:20 are identical. SD is skewed too low to fit a distribution.
#   Considered using gradient at quantiles to find 'flat' part of distribution. Nothing conclusive.

########################################################################################################

# FIT DOUBLE SPOT TO WHITE/GREY IMAGES                                                              ####

    spot1 <- spot.lm(pw.m[ , , "white", "160314"], o = 4, robust = T)
    spot1.res <- matrix(spot1$residuals, ncol = 1996)
    
        pixel.image(spot1.res)
    
    spot2 <- spot.lm(spot1.res, o = 4, robust = T)
    spot2.res <- matrix(spot2$residuals, ncol = 1996)
    
        pixel.image(spot1.res)

        par(mfrow = c(2,1))
            plot(c(0:1433), predict(spot1, data.frame(z = (c(0:1433)))), type = "l", xlab = "z", ylab = "")
            plot(c(0:1433), predict(spot2, data.frame(z = (c(0:1433)))), type = "l", xlab = "z", ylab = "")
        par(mfrow = c(1,1))
    
    panel <- panel.lm(spot2.res, "x + y", robust = T)
    res <- spot2.res - panel$fitted.values
    
        pixel.image(panel$fitted.values)
        pixel.image(res)

########################################################################################################
        
# REMOVE CONSTRAINT ON CENTRE POINT                                                                 ####
# FIT HYPER-ELLIPTICAL SPOTS (ellipse / squircle)                                                   ####
    
# Models fitted to images from 15-08-28
        
    # white, circular spot o = 2; MAD 265, SD 534
    {
        w.spot <- spot.lm(pw.m[,,"white", "150828"], o = 2, robust = T)
        w.spot.res <- matrix(w.spot$residuals, ncol = 1996)
        
        w.panel <- panel.lm(w.spot.res, "x + y", robust = T)
        w.res <- w.spot.res - w.panel$fitted.values
        mad(w.res); sd(w.res)   # 264.6319; 534.0722
    }
#--------------------------------------------------------------------------------------
    # white, n = 2, o = 2, x + y: MAD 262, SD 516
    {
        w.spot.2 <- he.spot.lm(pw.m[,,"white", "150828"], n = 2, order = 2, robust = T)
        w.spot.2.res <- matrix(w.spot.2$residuals, ncol = 1996)
        pixel.image(w.spot.2.res)
        plot(matrix(w.spot.2$fitted.values, ncol = 1996)[992,], type = "l")
        lines(matrix(w.spot.2$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
        
        w.panel.2 <- panel.lm(w.spot.2.res, "x + y", robust = T)
        w.2.res <- w.spot.2.res - w.panel.2$fitted.values
        mad(w.2.res); sd(w.2.res)   # 261.8727; 515.9726
    }
#-------------------------------------------------------------------------------------- 
    # white, n = 3, o = 2, x + y: MAD 288, SD 527
    {
        w.spot.3 <- he.spot.lm(pw.m[,,"white", "150828"], n = 3, order = 2, robust = T)
        w.spot.3.res <- matrix(w.spot.3$residuals, ncol = 1996)
        pixel.image(matrix(w.spot.3$fitted.values, ncol = 1996))
        
        plot(matrix(w.spot.3$fitted.values, ncol = 1996)[992,], type = "l")
        lines(matrix(w.spot.3$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
        
        w.panel.3 <- panel.lm(w.spot.3.res, "x + y", robust = T)
        w.3.res <- w.spot.3.res - w.panel.3$fitted.values
        mad(w.3.res); sd(w.3.res)   # 287.976; 527.1154
    }
#-------------------------------------------------------------------------------------- 
    # white, n = 4, o = 2, x + y: MAD 343, SD 553
    {
        w.spot.4 <- he.spot.lm(pw.m[,,"white", "150828"], n = 4, order = 2, robust = T)
        w.spot.4.res <- matrix(w.spot.4$residuals, ncol = 1996)
        pixel.image(matrix(w.spot.4$fitted.values, ncol = 1996))
        
        plot(matrix(w.spot.4$fitted.values, ncol = 1996)[992,], type = "l")
        lines(matrix(w.spot.4$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
        
        w.panel.4 <- panel.lm(w.spot.4.res, "x + y", robust = T)
        w.4.res <- w.spot.4.res - w.panel.4$fitted.values
        mad(w.4.res); sd(w.4.res)   # 342.3599; 552.5578
    }        
#-------------------------------------------------------------------------------------- 
    # plots
    {
            par(mfrow = c(2,2), mar = c(2,2,3,1))
                pixel.image(matrix(w.spot$fitted.values, ncol = 1996), title = paste0("Circular: MAD ", round(mad(w.res),0)))
                pixel.image(matrix(w.spot.2$fitted.values, ncol = 1996), title = paste0("n = 2: MAD ", round(mad(w.2.res),0)))
                pixel.image(matrix(w.spot.3$fitted.values, ncol = 1996), title = paste0("n = 3: MAD ", round(mad(w.3.res),0)))
                pixel.image(matrix(w.spot.4$fitted.values, ncol = 1996), title = paste0("n = 4: MAD ", round(mad(w.4.res),0)))
            par(mfrow = c(1,1))
            
            plot(matrix(w.spot$fitted.values, ncol = 1996)[992,], type = "l", main = "Fitted spot profiles (v)")
            lines(matrix(w.spot.2$fitted.values, ncol = 1996)[992,], type = "l", col = "blue")
            lines(matrix(w.spot.3$fitted.values, ncol = 1996)[992,], type = "l", col = "green3")
            lines(matrix(w.spot.4$fitted.values, ncol = 1996)[992,], type = "l", col = "red")
            legend("bottom", lty = 1, col = c("black","blue","green3","red"),
                   legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
            
            plot(matrix(w.spot$fitted.values, ncol = 1996)[,992], type = "l", main = "Fitted spot profiles (h)")
            lines(matrix(w.spot.2$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
            lines(matrix(w.spot.3$fitted.values, ncol = 1996)[,992], type = "l", col = "green3")
            lines(matrix(w.spot.4$fitted.values, ncol = 1996)[,992], type = "l", col = "red")
            legend("bottom", lty = 1, col = c("black","blue","green3","red"),
                   legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
        }
        
#============================================================================================
    # grey, circular spot o = 2; MAD 140, SD 383
    {
            g.spot <- spot.lm(pw.m[,,"grey", "150828"], o = 2, robust = T)
            g.spot.res <- matrix(g.spot$residuals, ncol = 1996)
            
            g.panel <- panel.lm(g.spot.res, "x + y", robust = T)
            g.res <- g.spot.res - g.panel$fitted.values
            mad(g.res); sd(g.res)   # 139.7952; 383.1336
        }
#--------------------------------------------------------------------------------------
    # grey, n = 2, o = 2, x + y: MAD 140, SD 380
    {
            g.spot.2 <- he.spot.lm(pw.m[,,"grey", "150828"], n = 2, order = 2, robust = T)
            g.spot.2.res <- matrix(g.spot.2$residuals, ncol = 1996)
           
            g.panel.2 <- panel.lm(g.spot.2.res, "x + y", robust = T)
            g.2.res <- g.spot.2.res - g.panel.2$fitted.values
            mad(g.2.res); sd(g.2.res)   # 140.1967; 380.3197
        }
#-------------------------------------------------------------------------------------- 
    # grey, n = 3, o = 2, x + y: MAD 145, SD 382
    {
            g.spot.3 <- he.spot.lm(pw.m[,,"grey", "150828"], n = 3, order = 2, robust = T)
            g.spot.3.res <- matrix(g.spot.3$residuals, ncol = 1996)
             
            g.panel.3 <- panel.lm(g.spot.3.res, "x + y", robust = T)
            g.3.res <- g.spot.3.res - g.panel.3$fitted.values
            mad(g.3.res); sd(g.3.res)   # 145.7017; 381.7986
        }
#-------------------------------------------------------------------------------------- 
    # grey, n = 4, o = 2, x + y: MAD 155, SD 385
    {
            g.spot.4 <- he.spot.lm(pw.m[,,"grey", "150828"], n = 4, order = 2, robust = T)
            g.spot.4.res <- matrix(g.spot.4$residuals, ncol = 1996)
            
            g.panel.4 <- panel.lm(g.spot.4.res, "x + y", robust = T)
            g.4.res <- g.spot.4.res - g.panel.4$fitted.values
            mad(g.4.res); sd(g.4.res)   # 154.732; 384.7382
        }        
#-------------------------------------------------------------------------------------- 
    # plots
    {
            par(mfrow = c(2,2), mar = c(2,2,3,1))
                pixel.image(matrix(g.spot$fitted.values, ncol = 1996), title = paste0("Circular: MAD ", round(mad(g.res),0)))
                pixel.image(matrix(g.spot.2$fitted.values, ncol = 1996), title = paste0("n = 2: MAD ", round(mad(g.2.res),0)))
                pixel.image(matrix(g.spot.3$fitted.values, ncol = 1996), title = paste0("n = 3: MAD ", round(mad(g.3.res),0)))
                pixel.image(matrix(g.spot.4$fitted.values, ncol = 1996), title = paste0("n = 4: MAD ", round(mad(g.4.res),0)))
            par(mfrow = c(1,1))
            
            plot(matrix(g.spot$fitted.values, ncol = 1996)[992,], type = "l", main = "Fitted spot profiles (v)")
            lines(matrix(g.spot.2$fitted.values, ncol = 1996)[992,], type = "l", col = "blue")
            lines(matrix(g.spot.3$fitted.values, ncol = 1996)[992,], type = "l", col = "green3")
            lines(matrix(g.spot.4$fitted.values, ncol = 1996)[992,], type = "l", col = "red")
            legend("bottom", lty = 1, col = c("black","blue","green3","red"),
                   legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
            
            plot(matrix(g.spot$fitted.values, ncol = 1996)[,992], type = "l", main = "Fitted spot profiles (h)")
            lines(matrix(g.spot.2$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
            lines(matrix(g.spot.3$fitted.values, ncol = 1996)[,992], type = "l", col = "green3")
            lines(matrix(g.spot.4$fitted.values, ncol = 1996)[,992], type = "l", col = "red")
            legend("bottom", lty = 1, col = c("black","blue","green3","red"),
                   legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
        }

        
#============================================================================================
    # black, circular spot o = 2; MAD 92, SD 413
    {
        b.spot <- spot.lm(pw.m[,,"black", "150828"], o = 2, robust = T)
        b.spot.res <- matrix(b.spot$residuals, ncol = 1996)
        
        b.panel <- panel.lm(b.spot.res, "x + y", robust = T)
        b.res <- b.spot.res - b.panel$fitted.values
        mad(b.res); sd(b.res)   # 92.26396; 412.6968
    }
#--------------------------------------------------------------------------------------
    # black, n = 2, o = 2, x + y: MAD 92, SD 413
    {
        b.spot.2 <- he.spot.lm(pw.m[,,"black", "150828"], n = 2, order = 2, robust = T)
        b.spot.2.res <- matrix(b.spot.2$residuals, ncol = 1996)
        
        b.panel.2 <- panel.lm(b.spot.2.res, "x + y", robust = T)
        b.2.res <- b.spot.2.res - b.panel.2$fitted.values
        mad(b.2.res); sd(b.2.res)   # 92.43132; 412.6959
    }
#-------------------------------------------------------------------------------------- 
    # black, n = 3, o = 2, x + y: MAD 91, SD 412
    {
        b.spot.3 <- he.spot.lm(pw.m[,,"black", "150828"], n = 3, order = 2, robust = T)
        b.spot.3.res <- matrix(b.spot.3$residuals, ncol = 1996)
        
        b.panel.3 <- panel.lm(b.spot.3.res, "x + y", robust = T)
        b.3.res <- b.spot.3.res - b.panel.3$fitted.values
        mad(b.3.res); sd(b.3.res)   # 91.26837; 412.3268
    }
#-------------------------------------------------------------------------------------- 
    # black, n = 4, o = 2, x + y: MAD 91, SD 412
    {
        b.spot.4 <- he.spot.lm(pw.m[,,"black", "150828"], n = 4, order = 2, robust = T)
        b.spot.4.res <- matrix(b.spot.4$residuals, ncol = 1996)
        
        b.panel.4 <- panel.lm(b.spot.4.res, "x + y", robust = T)
        b.4.res <- b.spot.4.res - b.panel.4$fitted.values
        mad(b.4.res); sd(b.4.res)   # 91.41425; 412.3623
    }        
#-------------------------------------------------------------------------------------- 
    # plots
    {
        par(mfrow = c(2,2), mar = c(2,2,3,1))
            pixel.image(matrix(b.spot$fitted.values, ncol = 1996), title = paste0("Circular: MAD ", round(mad(b.res),0)))
            pixel.image(matrix(b.spot.2$fitted.values, ncol = 1996), title = paste0("n = 2: MAD ", round(mad(b.2.res),0)))
            pixel.image(matrix(b.spot.3$fitted.values, ncol = 1996), title = paste0("n = 3: MAD ", round(mad(b.3.res),0)))
            pixel.image(matrix(b.spot.4$fitted.values, ncol = 1996), title = paste0("n = 4: MAD ", round(mad(b.4.res),0)))
        par(mfrow = c(1,1))
        
        plot(matrix(b.spot$fitted.values, ncol = 1996)[992,], type = "l", main = "Fitted spot profiles (v)")
        lines(matrix(b.spot.2$fitted.values, ncol = 1996)[992,], type = "l", col = "blue")
        lines(matrix(b.spot.3$fitted.values, ncol = 1996)[992,], type = "l", col = "green3")
        lines(matrix(b.spot.4$fitted.values, ncol = 1996)[992,], type = "l", col = "red")
        legend("top", lty = 1, col = c("black","blue","green3","red"),
               legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
        
        plot(matrix(b.spot$fitted.values, ncol = 1996)[,992], type = "l", main = "Fitted spot profiles (h)")
        lines(matrix(b.spot.2$fitted.values, ncol = 1996)[,992], type = "l", col = "blue")
        lines(matrix(b.spot.3$fitted.values, ncol = 1996)[,992], type = "l", col = "green3")
        lines(matrix(b.spot.4$fitted.values, ncol = 1996)[,992], type = "l", col = "red")
        legend("top", lty = 1, col = c("black","blue","green3","red"),
               legend = c("Circular", "n = 2", "n = 3", "n = 4"), bty = "n")
    }      
        
########################################################################################################
         
# FIT 32 x 4 SUBPANELS                                                                              ####

    spot <- spot.lm(pw.m[,,"black", "141009"], order = 2, robust = T)
    pixel.image(matrix(spot$fitted.values, ncol = 1996))
    spot.res <- matrix(spot$residuals, ncol = 1996)
    pixel.image(spot.res)
        
    panel <- panel.lm(spot.res, "x + y", robust = T)
    panel.res <- spot.res - panel$fitted.values
    pixel.image(panel.res)
    
    s.hist(panel.res)
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(panel.res, moment = "quant")), col = "red", lty = 2)
    qJohnson(c(0.001, 0.999), JF)       # -271.4574,  274.5265
    mad(panel.res); sd(panel.res)       # 89.41521, 407.9843
    
    mini.panel <- minipanel.lm(panel.res, terms = "x + y", robust = T)
    pixel.image(mini.panel$fitted.values)
    res <- panel.res - mini.panel$fitted.values
    mad(res); sd(panel.res)             # 77.55518, 405.022
    pixel.image(res)
    s.hist(res)
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(res, moment = "quant")), col = "red", lty = 2)
    
    # will actually pick up more points as bad pixels: over-fitting.
    
#-----------------------------------------------------------------------------------------
# compare to grey & white images...
    spot <- spot.lm(pw.m[,,"grey", "141009"], order = 2, robust = T)
    spot.res <- matrix(spot$residuals, ncol = 1996)
    panel <- panel.lm(spot.res, "x + y", robust = T)
    panel.res <- spot.res - panel$fitted.values
    mini.panel <- minipanel.lm(panel.res, terms = "x + y", robust = T)
    res <- panel.res - mini.panel$fitted.values
    mad(panel.res); sd(panel.res)       # 147.8663, 377.3608
    mad(res); sd(res)                   # 131.8415, 377.3608
    
    s.hist(panel.res)
    hist(res, add = T, breaks = "FD", col = adjustcolor("skyblue", alpha = 0.2), border = adjustcolor("skyblue", alpha = 0.2))
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(panel.res, moment = "quant")), col = "grey", lty = 2)
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(res, moment = "quant")), col = "skyblue", lty = 2)
    
#-----------------------------------------------------------------------------------------   
    spot <- spot.lm(pw.m[,,"white", "141009"], order = 2, robust = T)
    spot.res <- matrix(spot$residuals, ncol = 1996)
    panel <- panel.lm(spot.res, "x + y", robust = T)
    panel.res <- spot.res - panel$fitted.values
    mini.panel <- minipanel.lm(panel.res, terms = "x + y", robust = T)
    res <- panel.res - mini.panel$fitted.values
    mad(panel.res); sd(panel.res)       # 310.528, 542.4823
    mad(res); sd(res)             # 263.991, 542.4823
    
    s.hist(panel.res)
    hist(res, add = T, breaks = "FD", col = adjustcolor("skyblue", alpha = 0.2), border = adjustcolor("skyblue", alpha = 0.2))
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(panel.res, moment = "quant")), col = "grey", lty = 2)
    abline(v = qJohnson(c(0.001, 0.999), JohnsonFit(res, moment = "quant")), col = "skyblue", lty = 2)
    
########################################################################################################
    
# COMPARE MINIPANEL FIT WITH LINEAR VS QUADRATIC PANELS                                             ####
    
    spot <- spot.lm(pw.m[,,"black", "141009"], order = 2, robust = T)
    spot.res <- matrix(spot$residuals, ncol = 1996)
    panel <- panel.lm(spot.res, "x + y", robust = T)
    panel.res <- spot.res - panel$fitted.values
    mini.panel <- minipanel.lm(panel.res, terms = "x + y", robust = T)
    mp.res <- panel.res - mini.panel$fitted.values
    
    mad(panel.res); sd(panel.res)
    mad(mp.res); sd(mp.res)
    
    jpeg("./Plots/Minipanels/Minipanels-after-linear-panels-black-141009.jpg", width = 9000, height = 2400)
    par(mfrow = c(1,4), mar = c(2,2,4,1))
        pixel.image(panel$fitted.values, title = "Panels: x + y", cex.main = 10)
        pixel.image(panel.res, title = paste0("Panel residuals: MAD ", round(mad(panel.res),0), ", SD ", round(sd(panel.res),0)), cex.main = 10)
        pixel.image(mini.panel$fitted.values, title = paste0("Minipanels: x + y (SD ", round(sd(mini.panel$fitted.values),0), ")"), cex.main = 10)
        pixel.image(mp.res, title = paste0("Minipanel residuals: MAD ", round(mad(mp.res),0), ", SD ", round(sd(mp.res),0)), cex.main = 10)
    dev.off()

    # results (over 14-10-09 images)
    #           panel x + y      minipanel      panel poly(2)    minipanel
    #           sd      mad     sd      mad      sd      mad    sd      mad
    # black     408      89     405      78      406      82    405      77
    # grey      377     148     369     132      370     138    367     120
    # white     542     311     500     263      500     275    486     255
    

########################################################################################################
    
# FIT PANELS TO CROPPED DATA (EXCLUDE 100PX AROUND PANEL EDGES)                                     ####
    spot <- spot.lm(pw.m[ , , "white", "160314"], o = 2, robust = T)
    spot.res <- matrix(spot$residuals, ncol = 1996)
    
    cropped <- spot.res
    cropped[c(1:100, 1896:1996), c(1:100, 1896:1996)] <- NA
    
    panel <- panel.lm(cropped, "x + y", robust = T)
    colnames(panel$models) <- c("offset", "x", "y")
    fitted <- fit.panels(panel$models)

    pixel.image(fitted$offset + fitted$grad)

    res <- spot.res - fitted$offset - fitted$grad
    pixel.image(res)
    
########################################################################################################
# SD: JOHNSON DISTRIUTION OR LOGNORMAL?                                                             ####
fpath <- "./Models/Simple-parametric/"
    
res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))
load.pixel.sds()
    
JF.res <- read.csv(paste0(fpath, "JF-residuals.csv"), as.is = T)
JF.res <- JF.res[order(JF.res$batch),]
rownames(JF.res) <- apply(JF.res[,2:1], 1, paste, collapse = ".")
JF.res <- setNames(split(JF.res[,3:7], seq(nrow(JF.res))), rownames(JF.res[,3:7]))

JF.sd <- read.csv(paste0(fpath, "JF-sd.csv"), as.is = T)
JF.sd <- JF.sd[order(JF.sd$batch),]
rownames(JF.sd) <- apply(JF.sd[,2:1], 1, paste, collapse = ".")
JF.sd <- setNames(split(JF.sd[,3:7], seq(nrow(JF.sd))), rownames(JF.sd[,3:7]))

# SD came out as lognormal type. Try fitting LN model instead
p.log <- log(pw.sd[,,"black", "150828"])
p.log[pw.sd[,,"black", "150828"] == 0] <- 0

hist(p.log, breaks = "fd", prob = T)
lines(c(0:1000)/100, dnorm(c(0:1000)/100, mean = mean(p.log), sd = sd(p.log)), col = "orange", lwd = 2)
lines(c(0:1000)/100, dnorm(c(0:1000)/100, mean = mean(p.log), sd = mad(p.log)), col = "red", lwd = 2)

hist(pw.sd[,,"black", "150828"], breaks = "fd", prob = T, xlim = c(0,50))
lines(c(0:5000)/100, dlnorm(c(0:5000)/100, mean = mean(p.log), sd = sd(p.log)), col = "orange", lwd = 2)
lines(c(0:5000)/100, dlnorm(c(0:5000)/100, mean = mean(p.log), sd = mad(p.log)), col = "red", lwd = 2)
lines(c(0:5000)/100, dJohnson(c(0:5000)/100, JF.sd$"150828.black"), col = "green3", lwd = 2)
# log(SD) is fairly well fitted by normal distribution using MAD instead of standard deviation

# check in white & grey images also...
p.log.g <- log(pw.sd[,,"grey", "150828"])
p.log.g[pw.sd[,,"grey", "150828"] == 0] <- 0

hist(pw.sd[,,"grey", "150828"], breaks = "fd", prob = T, xlim = c(0,500))
lines(c(0:50000)/100, dlnorm(c(0:50000)/100, mean = mean(p.log.g), sd = sd(p.log.g)), col = "orange", lwd = 2)
lines(c(0:50000)/100, dlnorm(c(0:50000)/100, mean = mean(p.log.g), sd = mad(p.log.g)), col = "red", lwd = 2)
lines(c(0:50000)/100, dJohnson(c(0:50000)/100, JF.sd$"150828.grey"), col = "green3", lwd = 2)

# check in white & grey images also...
p.log.w <- log(pw.sd[,,"white", "150828"])
p.log.w[pw.sd[,,"white", "150828"] == 0] <- 0

hist(pw.sd[,,"white", "150828"], breaks = "fd", prob = T, xlim = c(0,500))
lines(c(0:50000)/100, dlnorm(c(0:50000)/100, mean = mean(p.log.w), sd = sd(p.log.w)), col = "orange", lwd = 2)
lines(c(0:50000)/100, dlnorm(c(0:50000)/100, mean = mean(p.log.w), sd = mad(p.log.w)), col = "red", lwd = 2)
lines(c(0:50000)/100, dJohnson(c(0:50000)/100, JF.sd$"150828.white"), col = "green3", lwd = 2)

ks.test(pw.sd[,,"black", "150828"], "pJohnson", parms = JF.sd$"150828.black")
ks.test(pw.sd[,,"grey", "150828"], "pJohnson", parms = JF.sd$"150828.grey")
ks.test(pw.sd[,,"white", "150828"], "pJohnson", parms = JF.sd$"150828.white")

ks.test(pw.sd[,,"black", "150828"], "plnorm", mean = mean(p.log), sd = sd(p.log))
ks.test(pw.sd[,,"black", "150828"], "plnorm", mean = mean(p.log), sd = mad(p.log))

ks.test(pw.sd[,,"grey", "150828"], "plnorm", mean = mean(p.log.g), sd = sd(p.log.g))
ks.test(pw.sd[,,"grey", "150828"], "plnorm", mean = mean(p.log.g), sd = mad(p.log.g))

ks.test(pw.sd[,,"white", "150828"], "plnorm", mean = mean(p.log.w), sd = sd(p.log.w))
ks.test(pw.sd[,,"white", "150828"], "plnorm", mean = mean(p.log.w), sd = mad(p.log.w))

# none of these appear to be a good fit. Need to understand exactly why

########################################################################################################
# PLOTS OF JOHNSON DISTRIBUTION OVER SD & MEAN: NATURAL CUTOFF?                                     ####
    fpath <- "./Models/Simple-parametric/"

    res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))
    load.pixel.sds()

    JF.res <- read.csv(paste0(fpath, "JF-residuals.csv"), as.is = T)
    JF.res <- JF.res[order(JF.res$batch),]
    rownames(JF.res) <- apply(JF.res[,2:1], 1, paste, collapse = ".")
    JF.res <- setNames(split(JF.res[,3:7], seq(nrow(JF.res))), rownames(JF.res[,3:7]))

    JF.sd <- read.csv(paste0(fpath, "JF-sd.csv"), as.is = T)
    JF.sd <- JF.sd[order(JF.sd$batch),]
    rownames(JF.sd) <- apply(JF.sd[,2:1], 1, paste, collapse = ".")
    JF.sd <- setNames(split(JF.sd[,3:7], seq(nrow(JF.sd))), rownames(JF.sd[,3:7]))

    JF.sim.res <- do.call("rbind", lapply(lapply(JF.res, rJohnson, n = 10000), sort))
    JF.sim.sd <- do.call("rbind", lapply(lapply(JF.sd, rJohnson, n = 10000), sort))
    
    #150108 not right - model predicts < 0 SD...
    {
        load.images(150108, "black")
        all(b.150108[,,1] == b.150108[,,10])
        sd <- apply(b.150108, c(1:2), sd)
        summary(c(sd))
        
        all(sd == pw.sd[,,"black", "150108"])
        
        s.hist(sd)
        
        length(which(sd == 0))
    }
    
    #---------------------------------------------------------------------------------------
    
    cols <- c("green4", "cyan3", "dodgerblue3", "blue2", "purple", "magenta3", "red",
              "orange", "gold", "green", "black")
    
    # plot shape of Johnson distribution
    {
        par(mfrow = c(5,2), mar = c(1,2,1,1))
        for (i in c(1:11)[-4]) {
            plot(JF.sim.sd[i,] - mean(JF.sim.sd[i,]), type = "l", col = cols[i], lwd = 2, ylim = c(-20,30))
            abline(h = zz[i,] - mean(JF.sim.sd[i,]), lty = 3, col = cols[i])
        }
        par(mfrow = c(1,1), mar = org.mar)
        # location changes - shape is basically same
    }


    # quick estimate of gradient around each point
    sd.quants <- do.call("rbind", lapply(JF.sd, qJohnson, p = sort(rep(c(0.0001, 0.01, 0.1, 0.9, 0.99, 0.999),3)) + rep(c(-.00001, 0, 0.00001), 6)))
    colnames(sd.quants) <- apply(cbind("q.", format(sort(rep(c(0.0001, 0.01, 0.1, 0.9, 0.99, 0.999),3)) + rep(c(-.00001, 0, 0.00001), 6), scientific = F)), 1, paste, collapse = "")

    plot(c(0.1:999.9)/1000, qJohnson(c(0.1:999.9)/1000, JF.sd[[1]]), pch = 20)
    sd.q.grad <- data.frame("g.0001" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0001)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0002)),
                            "g.001" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0010)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0011)),
                            "g.01" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0100)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.0101)),
                            "g.1" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.1000)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.1001)),
                            "g.9" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9000)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9001)),
                            "g.09" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9900)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9901)),
                            "g.009" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9990)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9991)),
                            "g.0009" = do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9998)) - 
                                do.call("rbind", lapply(JF.sd, qJohnson, p = 0.9999)))
    
    plot(0, type = "n", xlim = c(0,1), ylim = c(-2.5, 0))
    for (i in c(1:11)[-4]) {
        points(c(0.0001, 0.001, 0.01, 0.1, 0.9, 0.99, 0.999, 0.9999), sd.q.grad[i,], pch = 20, col = cols[i])
    }
    range(sd.q.grad[1:11,])
########################################################################################################
