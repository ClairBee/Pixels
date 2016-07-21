
# parametric modelling of expected spot shape

library("IO.Pixels"); library("CB.Misc")

acq <- readRDS("./02_Objects/images/pwm-loan.rds")
os.g <- acq[,,"grey"] - acq[,,"black"]
gv <- setNames(melt(os.g), nm = c("x", "y", "z"))

pw.m <- abind(sapply(c("131122", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)

# should probably normalise panels first - BUT this can assign variation to subpanels erroneously.
# keep simplest approach: fit circular spot & report deviations

####################################################################################################

# ELLIPTICAL SPOT                                                                               ####

# mystery hybrid function...
# elliptical, no eccentricity (independent x & y), exponential curve; manual function (ls, optim)
{
    ellipse.exp.ls <- function(obs, param) {
        A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
        sig.x <- param["sig.x"]; sig.y <- param["sig.y"]
        
        est <- A * exp(- 0.5 * ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2))
        sum((est - obs$z)^2, na.rm = T)
    }
    
    zz <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700),
                g2d.ls, obs = gv)
    gv$fv.g2d.ls <- zz$par["A"] * exp(-((((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + ((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2))
    fv.g2d.ls <- array(gv$fv.g2d.ls, dim = c(2048, 2048))
    res.g2d.ls <- os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048))
    
    # plot fitted model
    {
        pixel.image(array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Fitted values")
        pixel.image(os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Residuals")
        
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1:1024], ylim = c(0,1000), breaks = "fd", main = "Residuals after circular spot")
        hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1025:2048], breaks = "fd", 
             add = T, border = adjustcolor("red", alpha = 0.3))
        legend("topright", col = c("black", "red"), pch = 15, legend = c("Lower", "Upper"), bty = "n")
        
        image(1:2048, 1:2048, os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), asp = T,
              breaks = c(-20000, -500, 500, 20000), col = c("blue", NA, "red"), xlim = c(1900,2048), ylim = c(0,100))
        
        plot(os.g[1024,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1024,], col = "blue")
        
        plot(os.g[1025,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1025,], col = "green3")
        
        plot(os.g[, 1024], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[, 1024], col = "green3")
        
        plot(os.g[,1025], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[,1025], col = "green3")
        
        plot(os.g[512,], type = "l")
        lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[512,], col = "green3")
    }
    
    # now fit subpanels to model
    {
        panel.os <- lm(fv.g2d.ls ~ tt, data = gv[!is.na(gv$z),])
        pp <- array(dim = c(2048, 2048))
        pp[25:2024, 25:2024] <- panel.os$fitted.values - coef(panel.os)[1]
        pixel.image(pp)
        pixel.image(res.g2d.ls - pp)
        pp <- predict(panel.os, newdata = gv$tt)
    }
    
    # function to apply directly to image
    fit.ellipse <- function(im, A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700) {
        
        surf.dat <- setNames(melt(im), nm = c("x", "y", "z"))
        
        optim(c(A = A, x0 = x0, y0 = y0, sig.x = sig.x, sig.y = sig.y), ellipse.exp.ls, obs = surf.dat)
    }
    
    zz <- fit.ellipse(im = pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])
    fv.g2d.ls <- array(zz$par["A"] * exp(-((((gv$x - zz$par["x0"]) / zz$par["sig.x"])^2) + ((gv$y - zz$par["y0"]) / zz$par["sig.y"])^2)),
                       dim = c(2048, 2048))
    
    pixel.image(fv.g2d.ls)
    pixel.image(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"] - fv.g2d.ls)
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1024,], type = "l")
    lines(fv.g2d.ls[1024,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[1025,], type = "l")
    lines(fv.g2d.ls[1025,], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1024], type = "l")
    lines(fv.g2d.ls[,1024], col = "blue")
    
    plot((pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"])[,1025], type = "l")
    lines(fv.g2d.ls[,1025], col = "blue")
}


####################################################################################################

# 2D GAUSSIAN SPOT                                                                              ####
# allows for eccentricity through change in covariance matrix
bvn <- function(x, y, x0, y0, sig.x, sig.y, rho) {
    const <- (2 * pi * sig.x * sig.y * sqrt(1 - rho^2))^-1
    rr <- (2 * (1 - rho^2))^-1
    
    z.x <- (x - x0) / sig.x
    z.y <- (y - y0) / sig.y
    
    const * exp(- rr * (z.x^2 + z.y^2 - (2 * rho * z.x * z.y)))  
}

gauss.2d.ls <- function(param, x, y, z) {
    A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- A * bvn(x, y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - z)^2, na.rm = T)
}


# try MLE using built-in formula
gauss.2d.ml <- function(obs, param) {
    A <- param["A"]; x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- A + bvn(x = obs$x, y = obs$y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - obs$z)^2, na.rm = T)
    
    LL <- function(mu, sigma) {
        +     R = dnorm(x, mu, sigma)
        +     #
            +     -sum(log(R))
        + }
}


    
g2d <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5,
               sig.x = sd(gv$z, na.rm = T), sig.y = sd(gv$z, na.rm = T), rho = 0),
             gauss.2d.ls, x = gv$x, y = gv$y, z = gv$z, method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 2048, 2048, Inf, Inf, 1))
g2d$par
}

hh <- bvn(gv$x, gv$y, g2d$par["x0"], g2d$par["y0"], g2d$par["sig.x"], g2d$par["sig.y"], g2d$par["rho"])
pixel.image(array(hh, dim = c(2048, 2048)))

####################################################################################################

# MWE FOR STACK EXCHANGE                                                                        ####

# bivariate normal density function
bvn <- function(x, y, x0, y0, sig.x, sig.y, rho) {
    const <- 1 / (2 * pi * sig.x * sig.y * sqrt(1 - rho^2))
    rr <- 1 / (2 * (1 - rho^2))^-1
    
    z.x <- (x - x0) / sig.x
    z.y <- (y - y0) / sig.y
    
    const * exp(- rr * (z.x^2 + z.y^2 - (2 * rho * z.x * z.y)))  
}

# create some fake data
x <- seq(0, 500, 25); y <- seq(0, 500, 25)
z <- outer(x, y, bvn, x0 = 250, y0 = 250, sig.x = 100, sig.y = 100, rho = 0)

fdat <- setNames(melt(z), nm = c("x", "y", "z"))

gauss.2d.ls <- function(param, x, y, z) {
    x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- bvn(x, y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - z)^2, na.rm = T)
}

g2d <- optim(c(x0 = 240, y0 = 240, sig.x = 80, sig.y = 80, rho = 0),
             gauss.2d.ls, x = fdat$x, y = fdat$y, z = fdat$z, method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 250, 250, Inf, Inf, 1))
g2d$par

#           x0           y0        sig.x        sig.y          rho 
# 2.400000e+02 2.400000e+02 8.000000e+01 8.000000e+01 1.985233e-20


####################################################################################################

# QUADRATIC TREND SURFACE BY LS                                                                 ####

library(spatial)
s.ls <- surf.ls(2, gv[!is.na(gv$z), c("x", "y", "z")])
trsurf <- trmat(s.ls, 1, 2048, 1, 2048, 50)
pixel.image(trsurf$z)

# does my bivariate normal function fit successfully here?
g2d.ls <- optim(c(A = 15000, x0 = 25, y0 = 25,
               sig.x = sd(trsurf$z), sig.y = sd(trsurf$z), rho = 0),
             gauss.2d.ls, obs = setNames(melt(trsurf$z), nm = c("x", "y", "z")), method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 25, 25, Inf, Inf, 1))
g2d.ls$par

hh <- bvn(trsurf$x, trsurf$y, g2d.ls$par["x0"], g2d.ls$par["y0"], g2d.ls$par["sig.x"], g2d.ls$par["sig.y"], g2d.ls$par["rho"])
pixel.image(array(hh, dim = c(2048, 2048)))

# residual plot
pixel.image(os.g - trsurf$z)
o.plot(os.g[1024,])
lines(trsurf$z[1024,], col = "green3")

####################################################################################################

# 2D EXAMPLES                                                                                   ####

# Rosenbrock surface
{
    f <- function(x1,y1) (1-x1)^2 + 100*(y1 - x1^2)^2
    x <- seq(-2,2,by=.15)
    y <- seq(-1,3,by=.15)
    z <- outer(x,y,f)
    persp(x,y,z,phi=45,theta=-45,col="yellow",shade=.00000001,ticktype="detailed")
    
    # When using optim for multidimensional optimization, the input in your function definition must be a single vector
    f <- function(x) (1-x[1])^2 + 100*(x[2]-x[1]^2)^2
    
    # starting values must be a vector now
    optim( c(0,0), f )$par      # [1] 0.9999564 0.9999085
    
}

# Himmelblau's function
{
    f <- function(x1,y1) (x1^2 + y1 - 11)^2 + (x1 + y1^2 - 7)^2
    x <- seq(-4.5,4.5,by=.2)
    y <- seq(-4.5,4.5,by=.2)
    z <- outer(x,y,f)
    persp(x,y,z,phi=-45,theta=45,col="yellow",shade=.65 ,ticktype="detailed")
    
    f <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
    optim(c(-4,-4), f)$par
    optim(c(2,-2), f)$par
    optim(c(2,2), f)$par
    optim(c(-4,4),f)$par
}

# solving a maximisation problem
{
    f <- function(x,y) 2*x*(y**2)+2*(x**2)*y+x*y
    
    # set x & y values, get z
    x <- seq(-0.5,0.5, len=200)
    y <- seq(-0.5,0.5, len=200)
    z <- outer(x,y,f)
    image(x,y,z)
    
    # partial derivative functions
    fx <- function(x,y,h=0.001) (f(x+h,y)-f(x,y))/h
    fy <- function(x,y,h=0.001) (f(x,y+h)-f(x,y))/h
    
    zfx <- outer(x,y,fx)
    zfy <- outer(x,y,fy)
    
    contour(x,y,zfx,level=0)
    contour(x,y,zfy,level=0, add=T, col="red")
    
    # create new function for optim
    fbb<-function(x) f(x[1],x[2])
    optim(c(0.5,0.5),fbb,control=list(fnscale=-1))
    
    fxb <- function(x) fx(x[1],x[2])
    fyb <- function(x) fy(x[1],x[2])
    
    sumssq <- function(x) fxb(x)**2+fyb(x)**2       # sum of squares
    
    optim(c(0.1,0.1),sumssq)
    optim(c(-0.4,0.1),sumssq)
    
    optim(c(0.4,0.4),fbb,control=list(fnscale=-1))
}
####################################################################################################

# UPPER VS LOWER PANEL OFFSET                                                                   ####

# create factor to fit offset
gv$upper <- gv$y > 1024.5

ul.os <- lm(z ~ upper, data = gv[!is.na(gv$z),])
coef(ul.os)["upperTRUE"]

gv$p.adj <- gv$z - (gv$upper * coef(ul.os)["upperTRUE"])

pixel.image(array(gv$p.adj, dim = c(2048, 2048)))
o.plot(array(gv$z, dim = c(2048, 2048))[1024,])
o.plot(array(gv$z, dim = c(2048, 2048))[512,])

# or, maybe manually will be less mad - aim to smooth centre line by offsetting upper panel
# use median to remove effect of dead lines on either side
d <- apply(os.g[,1025 + c(0:19)], 1, median) - apply(os.g[,1024 -  + c(0:19)], 1, median)
median(d, na.rm = T)

gv$p.adj <- gv$z - (gv$upper * median(d, na.rm = T))

pixel.image(array(gv$p.adj, dim = c(2048, 2048)))
o.plot(array(gv$p.adj, dim = c(2048, 2048))[1024,])
o.plot(array(gv$p.adj, dim = c(2048, 2048))[512,])

gv2 <- setNames(gv[,c("x", "y", "p.adj")], nm = c("x", "y", "z"))
zz.adj <- optim(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 700, sig.y = 700),
            g2d.ls, obs = gv2)
zz.adj$par
zz$par
# plot fitted model
{
    gv$fv.g2d.ls <- zz.adj$par["A"] * exp(-((((gv$x - zz.adj$par["x0"]) / zz.adj$par["sig.x"])^2) + ((gv$y - zz.adj$par["y0"]) / zz.adj$par["sig.y"])^2))
    pixel.image(array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Fitted values")
    pixel.image(os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)), title = "Residuals")
    
    hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1:1024], breaks = "fd", ylim = c(0,30))
    hist((os.g - array(gv$fv.g2d.ls, dim = c(2048, 2048)))[,1025:2048], breaks = "fd", 
         add = T, border = adjustcolor("red", alpha = 0.3))
    
    plot(os.g[1024,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1024,], col = "green3")
    
    plot(os.g[1025,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[1025,], col = "green3")
    
    plot(os.g[, 1024], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[, 1024], col = "green3")
    
    plot(os.g[,1025], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[,1025], col = "green3")
    
    plot(os.g[512,], type = "l")
    lines(array(gv$fv.g2d.ls, dim = c(2048, 2048))[512,], col = "green3")
}

####################################################################################################

# PER-PANEL OFFSET                                                                              ####

# linear regression (simpler approach)
{
    d <- apply(os.g[,1025 + c(0:19)], 1, mean) - apply(os.g[,1024 -  + c(0:19)], 1, mean)
    o.plot(d)
    
    pixel.image(os.g[,1025:2048])
    
    sp <- array(dim = c(128, 16, 1024, 2), dimnames = list(NULL, NULL, NULL, c("L", "U")))
    
    for (ul in dimnames()) {
        for (p in 1:16) {
            sp[,p,,ul] <- paste0(ul, formatC(p, width = 2, flag = "0"))
        }
    }
    tt <- array(sp, dim = c(2048, 2048))
    
    gv$tt <- c(tt)
    
    panel.os <- lm(z ~ tt, data = gv[!is.na(gv$z),])
    
    pp <- array(predict(panel.os), dim = c(2000,2000))
    os.g.pa <- os.g
    os.g.pa[25:2024, 25:2024] <- os.g[25:2024, 25:2024] - pp
    pixel.image(os.g.pa)
}