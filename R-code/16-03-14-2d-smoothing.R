
library("IO.Pixels")
library(fields)

pw.m.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
n <- dim(pw.m.b)[3]
p <- panel.edges()

#########################################################################################
#                       TRY SMOOTHING DATA ACCROSS IMAGES                               #
#########################################################################################

# define panel start
x <- c(1023, 1150); y <- c(1, 992)

panel <- pw.m.b[x[1]:x[2], y[1]:y[2], 1]


#=====================================================================================
# spatstat::blur
par(mfrow = c(1,2))
pixel.image(as.im(panel))
pixel.image(blur(as.im(panel), bleed = T), title = "blur")
par(mfrow = c(1,1))

# need to find out how to handle edges. May be useful for quick visualisation later though


#=====================================================================================
# spatial::surf.ls
library(spatial)
library(MASS)       # for eqsplot
ls <- surf.ls(2, x = x[1]:x[2], y = y[1]:y[2], z = panel)
ls.tr <- trmat(ls, x[1], x[2], y[1], y[2], 1000)
eqscplot(ls.tr, type = "n")
contour(ls.tr, add = TRUE)

contour(ls.tr)

# example from documentation
data(topo, package="MASS")
topo.kr <- surf.ls(2, topo)
trsurf <- trmat(topo.kr, 0, 6.5, 0, 6.5, 50)
eqscplot(trsurf, type = "n")
contour(trsurf, add = TRUE)
points(topo)

eqscplot(trsurf, type = "n")
contour(trsurf, add = TRUE)
plot(topo.kr, add = TRUE)
title(xlab= "Circle radius proportional to Cook's influence statistic")


#=====================================================================================
# loess

# slow: takes several minutes even to run over single panel

library(rgl)    # needed for surface3d plot
library(reshape)    # needed for melt

l <- sd.levels(pw.m.b[,,1])

system.time(res <- loess(value ~ X1 + X2, data = melt(panel)))

zz <- matrix(res$fitted, ncol = ncol(panel))

par(mfrow = c(1,2))
pixel.image(panel, break.levels = l, title = "original")
pixel.image(zz,  break.levels = l, title = "loess")
par(mfrow = c(1,1))

# 3d representation not very clear - too steep
col.3d <- sd.colours()[(t(zz) - min(zz)) * (20/(max(zz)-min(zz))) + 1]
surface3d(x[1]: x[2], y[1]: y[2], zz, col = col.3d)

persp(c(x[1]:x[2]), c(y[1]:y[2]), zz, theta = 45, phi = 30, col = "green2", scale = FALSE,
      ltheta = -30, shade = 0.75, border = NA, box = FALSE)

# try some transects
tr.cols <- c("blue", "purple", "darkred", "red", "orangered", "orange", "gold", 
             "yellowgreen", "green", "cyan", "darkseagreen", "cornflowerblue")
c <- 4
s.cols <- matrix(ncol = 2, nrow = 13)

plot(zz[c,], type = "l", ylim = c(min(zz), max(zz)),
     main = "Transects across columns of Loess-smoothed panel",
     ylab = "Smoothed value", col = adjustcolor("darkblue", alpha = 0.5))
abline(coef(line(zz[c,])), col =  adjustcolor("darkblue", alpha = 0.3))
s.cols[1,] <- coef(line(zz[c,]))
       
for (i in 1:12) {
    points(zz[c + (10 * i),], type = "l", col = adjustcolor(tr.cols[i], alpha = 0.5))
    abline(coef(line(zz[c + (10 * i),])), col =  adjustcolor(tr.cols[i], alpha = 0.3))
    s.cols[i+1,] <- coef(line(zz[c + (10 * i),]))
}

#=====================================================================================
# could also use regression - simple linear model to get planar image?
# much quicker than loess - use while proving out concept
lm <- lm(value ~ X1 + X2, data = melt(panel))     # 0.171 elapsed
linear.smoothed <- predict(lm, zm[,1:2])   # 0.051 elapsed

z.lm <- matrix(linear.smoothed, ncol = ncol(panel))
par(mfrow = c(1,2))
pixel.image(panel, break.levels = l, title = "original")
pixel.image(z.lm,  break.levels = l, title = "linear regression")
spar(mfrow = c(1,1))

#--------------------------------------------------------------------------
# v quick & not terrible results on this panel. Test on whole image:
k <- 1
l <- sd.levels(pw.m.b[,,k])

panel.lm <- list()
smoothed.panels <- array(dim = c(1996, 1996))
lm.coeffs <- array(dim = c(32, 3), dimnames = list(NULL, c("Intercept", "X", "Y")))

for (i in c(2,1)) {
    for (j in c(1:16)) {
        m <- length(panel.lm) + 1
        melted.panel <- melt(pw.m.b[p$x[j] : (p$x[j+1]-1),
                                    p$y[i] : (p$y[i+1]-1), 
                                    k])
        
        panel.lm[[m]] <- lm(value ~ X1 + X2, data = melted.panel)
        lm.coeffs[m,] <- panel.lm[[m]]$coefficients
        
        smoothed.panels[p$x[j] : (p$x[j+1]-1), p$y[i] : (p$y[i+1]-1)] <- 
            predict(panel.lm[[m]], melted.panel[,1:2]) 
    }
}

pixel.image(smoothed.panels, break.levels = l)
panels.removed <- pw.m.b[,,1] - smoothed.panels
pixel.image(panels.removed)

o.plot(pw.m.b[x[1]+c, y[1]:y[2], k])
o.plot(panels.removed[x[1]+c, y[1]:y[2]], add = F, col = "red")


#=====================================================================================
# fields::smooth.2d
look <- smooth.2d( RMprecip$y,  x=RMprecip$x, theta=.25)
look3<-smooth.2d( RMprecip$y, x=RMprecip$x, theta=.25, nrow=256, 
                  ncol=256,Nwidth=32,
                  Mwidth=32)

image( look)
points( RMprecip$x, pch=".")

image( look3)
points( RMprecip$x, pch=".")

look2<- smooth.2d( RMprecip$y, x=RMprecip$x, cov.function=Exp.cov,theta=.25)
image(look2)
points( RMprecip$x, pch=".")
