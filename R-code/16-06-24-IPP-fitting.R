
# IN WHICH I SHALL FIT AN INHOMOGENEOUS POISSON PROCESS TO THE BAD PIXEL MAP

library("IO.Pixels"); library("CB.Misc")

# convert latest bad pixel map to point process
bp <- readRDS("./Notes/Final-classifications/fig/bad-px-by-feature-incl-local.rds")
bpx <- bp$"160430"
bpx <- bpx[bpx$f.type %in% c("cl.root", "singleton") & !bpx$type %in% c("l.bright", "l.dim"),]
bp.ppp <- ppp(bpx$row, bpx$col, c(1,1996), c(1,1996))

fpath <- "./Notes/Spatial/fig/"
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "blue", "green3")

####################################################################################################

# FUNCTIONS                                                                                     ####

env.plot <- function(px.ppp, px.ppm, dist.fun, normalise = F, ...) {
    
    if (normalise) {
        trans <- expression((. - pi * r ** 2))
    } else {
        trans <- NULL
    }
    
    plot(envelope(px.ppm, dist.fun, nsim = 99, nrank = 2, transform = trans, verbose = F), 
         col = "blue", legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), ...)
    plot(envelope(px.ppp, dist.fun, nsim = 99, nrank = 2, transform = trans, verbose = F), 
         add = T, legend = F, shadecol = adjustcolor("orange", alpha = 0.2))

    legend("topleft", bty = "n",
           pt.bg = adjustcolor(c("orange", "cyan3", NA, NA, NA), alpha = 0.2), 
           col = c(NA, NA, "red", "blue", "black"),
           pch = c(22, 22, NA, NA, NA), 
           lty = c(NA, NA, 2, 2, 1),
           legend = c("Envelope for data", "Envelope for model", "Theoretical (data)", "Theoretical (model)", "Observed"))
}


####################################################################################################

# NONPARAMETRIC - QUARTIC SMOOTHING                                                             ####

nonpara <- density(bp.ppp, sigma = bw.diggle(bp.ppp))    # bandwidth est using MSE approach
nonpara$v <- nonpara$v * 1996^2       # rescale intensity for easier interpretation
plot(nonpara)
contour(nonpara)

####################################################################################################

#  PARAMETRIC - QUADRATIC TREND                                                                ####

qt.ppm <- ppm(bp.ppp ~ x + y + I(x^2) + I(y^2) + I(x *y))

contour(predict(qt.ppm))
points(bp.ppp, pch = 20, cex = 0.8)
draw.panels(col = "grey", lty = 2)

env.plot(bp.ppp, qt.ppm, Fest, normalise = F, main = "")
env.plot(bp.ppp, qt.ppm, Gest, normalise = F, main = "")
env.plot(bp.ppp, qt.ppm, Kest, normalise = T, main = "")


