
library("IO.Pixels"); library("beepr")

pml <- gsub(".rds", "", gsub("pixel-map-","", list.files("./02_Objects/pixel-maps", pattern = "pixel-map")))

img.nm <- "160705"

############################################################################################################
# CHECKING FOR CSR                                                                                      ####

# load bad pixel map, remove all larger features & convert to Poisson point process
bpm <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", img.nm, ".rds"))
px.p <- bpm[bpm$f.type %in% c("cl.root", "singleton", "dense.region"),1:2]
px <- ppp(px.p$row, px.p$col, c(1,2048), c(1,2048))

pxp.s <- bpm[bpm$f.type %in% c("cl.root", "singleton"),1:2]      # remove dense regions
px.s <- ppp(pxp.s$row, pxp.s$col, c(1,2048), c(1,2048))

# hypothesis test: are quadrat counts in line with those we would expect under CSR?
qt <- quadrat.test(px, xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)

############################################################################################################
# INHOMOGENEOUS POINT PROCESSES                                                                         ####
#   Functions                                                                                           ####

# Poisson point process with quadratic trend
quadratic.ppm <- function(px) {
    ppm(px ~ x + y + I(x^2) + I(y^2) + I(x *y))
}

# convert point process coefficients into matrix of values
ppmfit.matrix <- function(ppm) {
    
    zz <- setNames(melt(array(dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))), nm = c("x", "y", "z"))
    cc <- coef(ppm)
    
    zz$z <- cc[1] + (zz$x * cc["x"]) + (zz$y * cc["y"]) + 
        (zz$x^2 * cc["I(x^2)"]) + (zz$y^2 * cc["I(y^2)"]) + (zz$x * zz$y * cc["I(x * y)"]) 
    
    array(zz$z, dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit point processes                                                                                   ####

# distance functions of raw point process
px.F <- Fest(px); px.G <- Gest(px); px.K <- Kest(px)
pxs.F <- Fest(px.s); pxs.G <- Gest(px.s); pxs.K <- Kest(px.s)

# quadratic IPP
qpp <- quadratic.ppm(px)
qpp.K <- Kinhom(px, qpp); qpp.G <- Ginhom(px, qpp); qpp.F <- Finhom(px, qpp)

qpp.s <- quadratic.ppm(px.s)
qpps.K <- Kinhom(px.s, qpp.s); qpps.G <- Ginhom(px.s, qpp.s); qpps.F <- Finhom(px.s, qpp.s)

############################################################################################################
# PDF SUMMARY                                                                                           ####

pdf(paste0("./00_June-2017/plots/", img.nm, "-spatial.pdf"), height = 4, width = 4 * 4); {
    
    print(paste0("Quadrat test statistic: ", round(qt$statistic, 0), " (p ~ ", round(qt$p.value, 4), ")"))
    
    par(mfrow = c(1,4), mar = c(2,2,2,1), oma = c(0,0,2,0))
    {
        plot(px.F, type = "n", main = "F-function")
        lines(px.F$r, px.F$rs, lwd = 2)
        lines(px.F$r, px.F$theo, lwd = 2, col = "red")
        lines(qpp.F$r, qpp.F$bord, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        plot(px.G, type = "n", main = "G-function")
        lines(px.G$r, px.G$rs, lwd = 2)
        lines(px.G$r, px.G$theo, lwd = 2, col = "red")
        lines(qpp.G$r, qpp.G$bord, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        plot(Kest(px), legend = F, lty = 1, lwd = 2, main = "K-function")
        lines(qpp.K$r, qpp.K$bord.modif, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        pixel.plot(px.p, main = "Quadratic trend")
        contour(1:2048, 1:2048, 2048^2 * exp(ppmfit.matrix(qpp)), nlevels = 30, col = "blue", add = T, lwd = 2)
    }
    mtext(paste0(img.nm, " quadrat test statistic: ", round(qt$statistic, 0), " (p ~ ", round(qt$p.value, 4), ")"), outer = T)

    {
        plot(pxs.F, type = "n", main = "F-function")
        lines(pxs.F$r, pxs.F$rs, lwd = 2)
        lines(pxs.F$r, pxs.F$theo, lwd = 2, col = "red")
        lines(qpps.F$r, qpps.F$bord, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        plot(pxs.G, type = "n", main = "G-function")
        lines(pxs.G$r, pxs.G$rs, lwd = 2)
        lines(pxs.G$r, pxs.G$theo, lwd = 2, col = "red")
        lines(qpps.G$r, qpps.G$bord, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        plot(Kest(px.s), legend = F, lty = 1, lwd = 2, main = "K-function")
        lines(qpps.K$r, qpps.K$bord.modif, col = "blue", lwd = 2)
        
        legend("topleft", lty = 1, col = c("black", "red", "blue"), lwd = 2, bty = "n",
               legend = c("Observed", "CSR", "Quadratic-trend IPP"))
    }
    {
        pixel.plot(pxp.s, main = "Quadratic trend")
        contour(1:2048, 1:2048, 2048^2 * exp(ppmfit.matrix(qpp.s)), nlevels = 30, col = "blue", add = T, lwd = 2)
    }
    mtext("Dense regions removed", outer = T)
}; dev.off()

