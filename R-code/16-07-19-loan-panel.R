
library("IO.Pixels"); library("CB.Misc")

acq <- readRDS("./02_Objects/images/pwm-loan.rds")

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

####################################################################################################

# IMPORT & RE-SAVE DATA                                                                         ####

acq <- import.acq("/home/clair/Documents/Pixels/Image-data/loan")
saveRDS(acq, "./02_Objects/images/pwm-loan.rds")

df <- summarise.profiles()
write.csv(df, "./Other-data/xml-profiles.csv", row.names = F, quote = F)

####################################################################################################

# EXPLORATORY                                                                                   ####

pixel.image(acq[,,"grey"] - acq[,,"black"])

# mean value in each image vs power/exposure required
apply(acq, 3, mean, na.rm = T)
apply(acq, 3, sd, na.rm = T)

# contour plots of each image
{
    contour(acq[,,"white"], nlevels = 30)
}
####################################################################################################

# DARK COLUMNS                                                                                  ####

# columns 407L, 413U, 2007L
px <- as.data.frame(which(acq[,,"white"] - acq[,,"black"] < 10000, arr.ind = T))
pixel.plot(px, cex = 0.4)

ddply(px, .(row), summarise, length = length(row))

# 407
{
    plot(acq[407,,"white"], type = "l", xlim = c(0,1024), ylim = c(4500, 5500), col = "gold")
    lines(acq[407,,"grey"], col = "green3")
    lines(acq[407,,"black"])
    
    # neighbouring columns
    plot(acq[406,,"white"] - acq[405,,"white"], type = "l", xlim = c(0,1024))   # 406 ok
    plot(acq[407,,"white"] - acq[406,,"white"], type = "l", xlim = c(0,1024))   # 407 dark
    plot(acq[408,,"white"] - acq[409,,"white"], type = "l", xlim = c(0,1024))   # 408 ok
}

# 413
{
    plot(acq[413,,"white"], type = "l", xlim = c(1024, 2048), ylim = c(4500, 5500), col = "gold")
    lines(acq[413,,"grey"], col = "green3")
    lines(acq[413,,"black"])
    
    # neighbouring columns
    plot(acq[412,,"white"] - acq[411,,"white"], type = "l", xlim = c(1024, 2048))   # 412 ok
    plot(acq[413,,"white"] - acq[412,,"white"], type = "l", xlim = c(1024, 2048))   # 413 dark
    plot(acq[414,,"white"] - acq[415,,"white"], type = "l", xlim = c(1024, 2048))   # 414 ok
}

# 2007
{
    plot(acq[2007,,"white"], type = "l", xlim = c(0, 1024), ylim = c(4500, 5500), col = "gold")
    lines(acq[2007,,"grey"], col = "green3")
    lines(acq[2007,,"black"])
    
    # neighbouring columns
    plot(acq[2006,,"white"] - acq[2005,,"white"], type = "l", xlim = c(0,1024))   # 412 ok
    plot(acq[2007,,"white"] - acq[2006,,"white"], type = "l", xlim = c(0,1024))   # 2007 dark
    plot(acq[2008,,"white"] - acq[2009,,"white"], type = "l", xlim = c(0,1024))   # 414 ok
}

# all together
{
    range(acq[407,1:1024, ], acq[413,1025:2048, ], acq[2007,1:1024, ], na.rm = T)
    
    plot(acq[407,1:1024, "black"], ylim = c(4500, 5500), type = "l", col = "forestgreen")
    lines(acq[407,1:1024, "grey"], col = adjustcolor("green3", alpha = 0.4))
    lines(acq[407,1:1024, "white"], col = adjustcolor("green", alpha = 0.4))
    
    lines(acq[413,2048:1025, "black"], col = adjustcolor("darkblue", alpha = 0.4))
    lines(acq[413,2048:1025, "grey"], col = adjustcolor("blue", alpha = 0.4))
    lines(acq[413,2048:1025, "white"], col = adjustcolor("cornflowerblue", alpha = 0.4))
    
    lines(acq[2007,1:1024, "black"], col = adjustcolor("darkred", alpha = 0.4))
    lines(acq[2007,1:1024, "grey"], col = adjustcolor("red", alpha = 0.4))
    lines(acq[2007,1:1024, "white"], col = adjustcolor("orange", alpha = 0.4))
}

# DARK ROW!!! (row 77 in black images)                                                          ####

pdf("./dark-row.pdf")
pixel.image(acq[,,"black"], xlim = c(0,512), ylim = c(0,100))
dev.off()

pdf("./dark-row-transect.pdf")
plot(acq[,77,"black"] - acq[,76,"black"], type = "l")
lines(acq[,76,"black"] - acq[,75,"black"], col = adjustcolor("blue", alpha = 0.4))
dev.off()

plot(acq[,77,"white"] - acq[,76,"white"], type = "l", col = "gold")
lines(acq[,77,"grey"] - acq[,76,"grey"], col = adjustcolor("green3", alpha = 0.4))
lines(acq[,77,"black"] - acq[,76,"black"])
abline(h = 0, col = adjustcolor("darkred", alpha = 0.7))

# transect plots
{
    plot(acq[,77,"black"] - acq[,76,"black"], type = "l")
    lines(acq[,76,"black"] - acq[,75,"black"], col = adjustcolor("blue", alpha = 0.4))
    
    plot(acq[,77,"black"], type = "l")
    lines(acq[,76, "black"], col = "green3")
    lines(acq[,78, "black"], col = "cyan3")
    
    plot(acq[,77,"grey"], type = "l")
    lines(acq[,76, "grey"], col = "green3")
    lines(acq[,78, "grey"], col = "cyan3")
    
    plot(acq[,77,"white"], type = "l")
    lines(acq[,76, "white"], col = "green3")
    lines(acq[,78, "white"], col = "cyan3")
}

# check whether it appears in all exposures
{
    # import all 20 black acquisitions
    ims <- list.files("/home/clair/Documents/Pixels/Image-data/loan/black",
                          pattern = "\\.tif$", full.names = T)
    tmp <- array(dim = c(2048, 2048, 20))
    tmp[25:2024, 25:2024,] <- abind(lapply(ims, function(im) t(readTIFF(im, as.is = T)[2000:1,,drop = F])), along = 3)


    apply(tmp, 3, function(im) plot(im[,77] - im[,76], type = "l"))
    
}



####################################################################################################

# PARAMETRIC MODELLING - QUADRATIC TREND                                                        ####

# black images - gradient per subpanel
{
    sp <- array(acq, dim = c(128, 16, 1024, 2))
    
    sp.lm <- apply(sp, 4, apply, 2, 
                   function(s) {
                       lm(value ~  x + y + I(x^2) + I(y^2) + I(x * y), 
                          setNames(melt(s), nm = c("x", "y", "value")))
                   })
    sp.fv <- array(dim = c(2048, 2048))
    sp.fv[25:2024, 25:2024] <- cbind(do.call("rbind", lapply(sp.lm[[1]], function(mm) matrix(fitted.values(mm), ncol = 1000))),
                                     do.call("rbind", lapply(sp.lm[[2]], function(mm) matrix(fitted.values(mm), ncol = 1000))))
    pixel.image(sp.fv)
    sp.res <- acq[,,"black"] - sp.fv
    .smoothScatter(sp.fv, sp.res)
    pixel.plot(which(sp.res > 1000, arr.ind = T), cex = 0.4)
    pixel.plot(which(sp.res < -1000, arr.ind = T), cex = 0.4)
    hist(sp.res, breaks = "fd", ylim = c(0,30))
}

# black images - universal x-gradient, y-gradient per subpanel
{

}

# grey image - spot + per-panel gradient
{
    # correct offset by subtracting black image
    os.g <- acq[,,"grey"] - acq[,,"black"]
    pixel.image(os.g)
    
    # fit quadratic trend
    spot.lm <- lm(value ~ x + y + I(x^2) + I(y^2) + I(x * y), 
                  setNames(melt(os.g), nm = c("x", "y", "value")))
    spot.fv <- array(dim = dim(os.g))
    spot.fv[25:2024, 25:2024] <- spot.lm$fitted.values
    pixel.image(spot.fv)
    
    pixel.image(os.g - spot.fv)
    hist(os.g - spot.fv, )
}

####################################################################################################

# PARAMETRIC MODELLING - 2D GAUSSIAN                                                            ####

# correct offset by subtracting black image
os.g <- acq[,,"grey"] - acq[,,"black"]
gv <- setNames(melt(os.g), nm = c("x", "y", "z"))

o.plot(os.g[1024,])
o.plot(os.g[1025,])
o.plot(os.g[,1024])
o.plot(os.g[,1025])

# grey image - 2d gaussian spot by least squares
{
    g2d.ls <- function(gv, par) {
        gv <- gv[!is.na(gv$z),]
        c <- par[1]; x0 <- par[2]; y0 <- par[3]
        sig.x <- par[4]; sig.y <- par[5]
        
        est <- exp(-((gv$x - x0)^2 / (2 * sig.x^2) + (gv$y-y0)^2 / (2*sig.y^2)))
        
        # get square of difference - this is what we minimize
        sum((est - gv$z)^2)
    }
    
    # ~211 seconds to complete
    system.time(zz <- optim(c(c = 1, x0 = 1024.5, y0 = 1024.5, sig.x = 1000, sig.y = 1000), g2d.ls, gv = gv))
    zz$par
    #          c         x0         y0      sig.x      sig.y 
    # 15866.6004   838.2912  3684.1408  2993.5229 10677.9436
    
    # c          x0          y0       sig.x       sig.y 
    # -128797.335    1014.828    1876.670  140988.389  122480.528 
    
    # plot fitted model
    g2d.fv <- setNames(melt(os.g)[,1:2], nm = c("x", "y"))
    g2d.fv$fv <- (zz$par[1] / (zz$par[4] * zz$par[5])) * exp(-((g2d.fv$x - zz$par[2])^2 / (2 * zz$par[4]^2) + (g2d.fv$y-zz$par[3])^2 / (2*zz$par[5]^2)))
    g2d.fitted <- array(g2d.fv$fv, dim = c(2048, 2048))
    pixel.image(g2d.fitted)
    pixel.image(os.g)
}

# grey image - ellipse by least squares
{
    ellipse.ls <- function(values, par) {
        values <- values[!is.na(values$z),]
        x0 <- par[1]; y0 <- par[2]; a <- par[3]; b <- par[4]; c <- par[5]

        est <- c + ((values$x - x0) / a)^2 + ((values$y-y0) / b)^2
        
        # get square of difference - this is what we minimize
        sum((est - values$z)^2)
    }
 
    zz <- optim(c(x0 = 1024.5, y0 = 1024.5, a = 1, b = 1, c = mean(gv$x, na.rm = T)), ellipse.ls, values = gv)
    ellipse.fv <- array(((gv$x - zz$par["x0"]) / zz$par["a"])^2 + ((gv$y-zz$par["y0"]) / zz$par["b"])^2, dim = c(2048, 2048))
    
    #   x0          y0           a           b           c 
    # 1054.862398  913.471142    7.847946    8.464464 1174.991253 
    
    pixel.image(ellipse.fv)
    pixel.image(os.g)
    ellipse.res <- os.g - ellipse.fv
    pixel.image(ellipse.res)
    
    mean(ellipse.res, na.rm = T)
    o.plot(os.g[1024,])
    o.plot(ellipse.fv[1024,], col = "blue")
    
}

# toy examples for practice
{
    min.RSS <- function(data, par) {
        with(data, sum((par[1] + par[2] * x - y)^2))
    }
    
    dat <- data.frame(x = c(1,2,3,4,5,6), y = c(1,3,5,6,8,12))
    
    result <- optim(par = c(0, 1), min.RSS, data = dat)
    # I find the optimised parameters in result$par
    # the minimised RSS is stored in result$value
    result
    ## $par
    ## [1] -1.267  2.029
    ## 
    ## $value
    ## [1] 2.819
    ## 
    ## $counts
    ## function gradient 
    ##       89       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
}

####################################################################################################

# LARGE CLUSTER DEFECTS                                                                         ####

