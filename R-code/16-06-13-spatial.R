
library("IO.Pixels"); library("CB.Misc")

#bp <- readRDS("./Notes/Final-classifications/fig/bad-px-by-feature.rds")
bp <- readRDS("./Notes/Final-classifications/fig/bad-px-by-feature-incl-local.rds")

fpath <- "./Notes/Spatial/fig/"
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "blue", "green3")

set.seed(24747)

####################################################################################################

# ENVELOPE FUNCTIONS                                                                            ####

plot.dist.ests <- function(bpx, dt, excl = c("line.b", "line.d"), file.id, cc = Cat.cols,
                           cl.excl = c("cl.body", "line.body"), r = NULL, ...) {

    dt <- toString(dt)
    
    px <- bpx[[dt]]
    
    px <- px[!(px$f.type %in% cl.excl),]

    qt <- quadrat.test(ppp(px$row[!(px$type %in% excl)], px$col[!(px$type %in% excl)], c(1,1996), c(1,1996)),
                 xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)
    
    if (qt$statistic < 999) {dp <- 2} else {dp <- 0}
    if (round(qt$p.value, 5) == 0) {ss <- "\\approx"} else {ss <- "="}
    
    qt.out <- paste0("$\\chi^2 = ", round(qt$statistic, dp), ", p ", ss, " ", round(qt$p.value, 5), "$")
    write(qt.out, paste0(fpath, file.id, "-", dt, ".txt"))
    
    # plot all bad pixels
    pdf(paste0(fpath, file.id, "-", dt, "-plot.pdf")); {
        plot(px[!(px$type %in% excl),1:2], 
             col = adjustcolor(cc[px[!(px$type %in% excl),"type"]], alpha = 0.5), 
             pch = 20, asp = T, xlab = "", ylab = "")
        dev.off()
    }
    cat("Bad pixel plot created.")
    
    # various intensity functions
    pdf(paste0(fpath, file.id, "-", dt, "-Kenv.pdf")); {
        plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                          px[!(px$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Kest, nsim = 99, nrank = 2, r = r), 
             main = "")
        dev.off()
    }
    cat("K-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Fenv.pdf")); {
        plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                          px[!(px$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Fest, nsim = 99, nrank = 2, r = r), 
             main = "", ...)
        dev.off()
    }
    cat("F-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Genv.pdf")); {
        plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                          px[!(px$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Gest, nsim = 99, nrank = 2, r = r), 
             main = "", ...)
        dev.off()
    }
    cat("G-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Henv.pdf")); {
        plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                          px[!(px$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Hest, nsim = 99, nrank = 2, r = r), 
             main = "", ...)
        dev.off()
    }
    cat("H-function plot created.")
    
  }


####################################################################################################

# ENVELOPES & QUADRAT TESTS                                                                     ####

r <- c(0:500)
# plots of latest image
{
    # including 'locally bright'/'locally dim' pixels, all cluster types
    plot.dist.ests(bp, 160430, excl = c("line.b", "line.d"), file.id = "incl-l", cl.excl = "", r = r, xlim = c(0,500))
    
    # excluding 'locally bright'/'locally dim' pixels
    plot.dist.ests(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "excl-l", cl.excl = "")
    
    # plots of latest image, including ONLY 'locally bright'/'locally dim' pixels
    plot.dist.ests(bp, 160430, excl = c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "v.dim", "dim"), file.id = "only-l", cl.excl = "")

    # plots of latest image, including 'locally bright'/'locally dim' pixels
    plot.dist.ests(bp, 160430, excl = c(""), file.id = "even-lines", cl.excl = "")
}

# cluster roots vs singletons
{
    plot.dist.ests(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "cl-roots", cl.excl = c("singleton", "cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "cl-singles", cl.excl = c("cl.root", "cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "cl-roots-and-singles", cl.excl = c("cl.body", "line.body"))
}

# plot per bad pixel type
{
    Cat <- c("no.resp" , "dead" , "hot" , "v.bright" , "bright" , "line.b" , "edge" , "l.bright" , "line.d" , "screen.spot" , "v.dim" , "dim" , "l.dim")
    
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "hot"], file.id = "cl-hot", cl.excl = c("cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "v.bright"], file.id = "cl-vbright", cl.excl = c("cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "bright"], file.id = "cl-bright", cl.excl = c("cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "no.resp"], file.id = "cl-no-resp", cl.excl = c("cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "l.dim"], file.id = "cl-ldim", cl.excl = c("cl.body", "line.body"))
    plot.dist.ests(bp, 160430, excl = Cat[Cat != "l.bright"], file.id = "cl-lbright", cl.excl = c("cl.body", "line.body"))
}

# try subpanels + minipanels (each subpanel divided vertically into 2) - df 63
{
    mp <- c(1, 993 + c(-512, 0, 512), 1997)
    
    quadrat.test(bp.ppp(bp, 160430, excl = c("line.b", "line.d")),
                 xbreaks = panel.edges()$x - 0.5, ybreaks = mp - 0.5)       # x2 = 25938, p = 0
    
    quadrat.test(bp.ppp(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim")),
                 xbreaks = panel.edges()$x - 0.5, ybreaks = mp - 0.5)       # x2 = 281.21, p = 0
    
    quadrat.test(bp.ppp(bp, 141009, excl = c("line.b", "line.d")),
                 xbreaks = panel.edges()$x - 0.5, ybreaks = mp - 0.5)       # x2 = 4443.8, p = 0
    
    # expected counts are only 7 - chi^2 approximation likely to be inaccurate
    quadrat.test(bp.ppp(bp, 141009, excl = c("line.b", "line.d", "l.bright", "l.dim")),
                 xbreaks = panel.edges()$x - 0.5, ybreaks = mp - 0.5)       # x2 = 177.39, p = 0
    
}
# also test for homogeneity within subpanels?

# focal plots of interesting pixels - particularly the triangular shapes around column 1000/900
{
    bp$"160430"[bp$"160430"$type == "hot" & bp$"160430"$f.type == "singleton",]
    plot(bpx[,1:2], pch = c(20, 15, 0, 20)[bpx$f.type], col = Cat.cols[bpx$type], 
         xlim = c(950, 1050), ylim = c(900, 1000))
    focal.plot(pw.m[ , , "black", "160430"], centre = c(1000, 943), surround = 20, bad.px = bp$"160430", pt.cex = 0.5, cex.main = 0.8)
    
}

####################################################################################################

# PARAMETRIC MODELLING - LOG-LAMBDA                                                             ####

bp.ppp <- function(px, im.dim = c(1996, 1996)) {
    
    # can't cope with ordered factor as input
    ppp(px$row, px$col, c(1,im.dim[1]), c(1,im.dim[2]), marks = factor(as.character(px$type)))
}

bpx.pp <- bp.ppp(bp$"160430"[bp$"160430"$type == "hot",])

summary(bpx.pp)

df <- rbind(data.frame(type = "hot", 
                       t(coef(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "hot",]) ~ x + y + I(x^2) + I(y^2) + I(x *y))))),
            data.frame(type = "v.bright", 
                       t(coef(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "v.bright",]) ~ x + y + I(x^2) + I(y^2) + I(x *y))))))

# convert coefficients to matrix
ppmfit.matrix <- function(ppm) {
    
    zz <- setNames(melt(array(dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))), nm = c("x", "y", "z"))
    cc <- coef(ppm)
    
    zz$z <- cc[1] + (zz$x * cc["x"]) + (zz$y * cc["y"]) + 
        (zz$x^2 * cc["I(x^2)"]) + (zz$y^2 * cc["I(y^2)"]) + (zz$x * zz$y * cc["I(x * y)"]) 
    
    array(zz$z, dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))
}


nl <- 20
# plot lambda
{
    pdf(paste0(fpath, "lambda-contour-hot.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "hot",]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "magenta3")
        points(bp$"160430"[bp$"160430"$type == "hot",1:2], pch = 20)
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    pdf(paste0(fpath, "lambda-contour-vbright.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "v.bright",]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "red")
        points(bp$"160430"[bp$"160430"$type == "v.bright",1:2], pch = 20)
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    pdf(paste0(fpath, "lambda-contour-bright.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "bright",]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "orange")
        points(bp$"160430"[bp$"160430"$type == "bright",1:2], pch = 20)
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    pdf(paste0(fpath, "lambda-contour-lbright.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "l.bright",]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "gold")
        #points(bp$"160430"[bp$"160430"$type == "l.bright",1:2], pch = 20, cex = 0.2)
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    
    pdf(paste0(fpath, "lambda-contour-clusters.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$f.type == "cl.root",]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "slateblue1")
        points(bp$"160430"[bp$"160430"$f.type == "cl.root",1:2], pch = 20,
               col = Cat.cols[bp$"160430"[bp$"160430"$f.type == "cl.root","type"]])
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    pdf(paste0(fpath, "lambda-contour-singletons.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$f.type == "singleton" & !(bp$"160430"$type %in% c("l.bright", "l.dim")),]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "cyan3")
        points(bp$"160430"[bp$"160430"$f.type == "singleton" & !(bp$"160430"$type %in% c("l.bright", "l.dim")),1:2],
               pch = 20, col = Cat.cols[bp$"160430"[bp$"160430"$f.type == "singleton" & !(bp$"160430"$type %in% c("l.bright", "l.dim")),"type"]])
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
    
    pdf(paste0(fpath, "lambda-contour-clusters-and-singles.pdf")); {
        par(mar = c(2,2,1,1))
        contour(1:1996, 1:1996, 
                1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$f.type %in% c("singleton", "cl.root") & !(bp$"160430"$type %in% c("l.bright", "l.dim")),]) ~ 
                                                   x + y + I(x^2) + I(y^2) + I(x *y)))),
                nlevels = nl, col = "blue")
        points(bp$"160430"[bp$"160430"$f.type  %in% c("singleton", "cl.root") & !(bp$"160430"$type %in% c("l.bright", "l.dim")),1:2],
               pch = 20, col = Cat.cols[bp$"160430"[bp$"160430"$f.type  %in% c("singleton", "cl.root") & !(bp$"160430"$type %in% c("l.bright", "l.dim")),"type"]])
        draw.panels(lty = 3, col = "grey")
        dev.off()
    }
}





plot(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "v.bright",]) ~ 
             x + y + I(x^2) + I(y^2) + I(x *y)), pch = 4, se = F, col = topo.colors(51))

plot(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "bright",]) ~ 
             x + y + I(x^2) + I(y^2) + I(x *y)), pch = 4, se = F, col = topo.colors(51))

plot(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "no.resp",]) ~ 
             x + y + I(x^2) + I(y^2) + I(x *y)), pch = 4, se = F, col = topo.colors(51))

plot(ppm(bp.ppp(bp$"160430"[bp$"160430"$f.type == "cl.root",]) ~ 
             x + y + I(x^2) + I(y^2) + I(x *y)), pch = 4, se = F, col = topo.colors(51))

# need to work out which of these is most appropriate way to display the data
plot(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "hot",]) ~ 
             x + y + I(x^2) + I(y^2) + I(x *y)), pch = 4, se = F, col = topo.colors(51),
     main = "Hot pixels")

tt <- ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "hot",]), 
           ~ x + y + I(x^2) + I(y^2) + I(x *y), DiggleGatesStibbard(10))
plot(tt, trend = T, cif = F, se = F, main = "DiggleGatesStibbard(10)")

####################################################################################################

# NONPARAMETRIC MODELLING                                                                       ####

# taken from p174
{
    data(bodmin)
    Mse2d <- mse2d(as.points(bodmin), bodmin$poly, nsmse=50, range=8)
    plot(Mse2d$h[5:50],Mse2d$mse[5:50], type="l")
    points(Mse2d$h[which.min(Mse2d$mse)], Mse2d$mse[which.min(Mse2d$mse)], col = "red")
    
    plot(bodmin)
    image(kernel2d(as.points(bodmin), bodmin$poly, h0=2, nx=100, ny=100))
    image(kernel2d(as.points(bodmin), bodmin$poly, h0=3, nx=100, ny=100))
    points(bodmin, pch = 20)
}

bpx <- bp$"160430"[bp$"160430"$f.type %in% c("cl.root", "singleton"),]
{
    # find bandwidth candidates using MSE approach
    bp.mse <- mse2d(as.points(list(x = bpx$row, y = bpx$col)),
                    as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                    nsmse = 50, range = 500)
    
    plot(bp.mse$h,bp.mse$mse, type="l")
    # v. flat, could use quite a wide range of bandwidths (minimised at 10)
    points(bp.mse$h[which.min(bp.mse$mse)], bp.mse$mse[which.min(bp.mse$mse)], col = "red")
    title(paste0("Minimised at ", bp.mse$h[which.min(bp.mse$mse)]))
    
    image(kernel2d(as.points(list(x = bpx$row, y = bpx$col)),
                   as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                   h0 = 7, nx = 100, ny = 100))
    # basically a bad pixel map - dominated by lines & local brightness/dimness
    
    image(kernel2d(as.points(list(x = bpx$row, y = bpx$col)),
                   as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                   h0 = 50, nx = 100, ny = 100))
    
    image(kernel2d(as.points(list(x = bpx$row, y = bpx$col)),
                   as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                   h0 = 100, nx = 100, ny = 100))
    # v. similar to contour fitted to locally bright pixels
    contour(1:1996, 1:1996, 
            1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bp$"160430"[bp$"160430"$type == "l.bright",]) ~ 
                                               x + y + I(x^2) + I(y^2) + I(x *y)))),
            nlevels = 20, add = T)
}


# remove locally bright/dim pixels
bpx <- bp$"160430"[bp$"160430"$f.type %in% c("cl.root", "singleton") & 
                       !(bp$"160430"$type %in% c("l.bright", "l.dim")),]
{
    bp.mse <- mse2d(as.points(list(x = bpx$row, y = bpx$col)),
                    as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                    nsmse = 50, range = 1000)
    
    plot(bp.mse$h,bp.mse$mse, type="l")
    # v. flat, could use quite a wide range of bandwidths (minimised at 120)
    points(bp.mse$h[which.min(bp.mse$mse)], bp.mse$mse[which.min(bp.mse$mse)], col = "red")
    title(paste0("Minimised at ", bp.mse$h[which.min(bp.mse$mse)]))
    
    image(kernel2d(as.points(list(x = bpx$row, y = bpx$col)),
                   as.points(list(x = c(0,0,1996,1996), y = c(0,1996,1996,0))),
                   h0 = 120, nx = 100, ny = 100))
    # overlay parametric contour plot
    contour(1:1996, 1:1996, 
            1996^2 * exp(ppmfit.matrix(ppm(bp.ppp(bpx) ~ x + y + I(x^2) + I(y^2) + I(x *y)))),
            nlevels = 20, col = "blue", add = T)
}

# OLD CODE #################################################################################### ####

# 7.3 PRELIMINARY ANALYSIS OF A POINT PATTERN
{
library(spatstat)

pdf(paste0(fpath, "plot-jpines.pdf")); {
    par(mar = c(0,0,0,0))
    plot(japanesepines, pch = 16, main = "")     # pattern compatible with CSR
    dev.off()
}
pdf(paste0(fpath, "plot-cells.pdf")); {
    par(mar = c(0,0,0,0))
    plot(cells, pch = 16, main = "")     # regular
    dev.off()
}
pdf(paste0(fpath, "plot-redwoods.pdf")); {
    par(mar = c(0,0,0,0))
    plot(redwoodfull, pch = 16, main = "")     # clustered
    dev.off()
}

# G function: distance to the nearest event
{
    envjap <- envelope(as(japanesepines, "ppp"), fun = Gest,
                       nrank = 2, nsim = 99)
    envred <- envelope(as(redwoodfull, "ppp"), fun = Gest, 
                       nrank = 2, nsim = 99)
    envcells <- envelope(as(cells, "ppp"), fun = Gest,
                         nrank = 2, nsim = 99)
    
    pdf(paste0(fpath, "Gf-jpines.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(envjap, main = "")     # pattern compatible with CSR
        dev.off()
    }
    
    plot(envjap$obs, envjap$theo, type = "l")
    lines(envjap$obs, envjap$lo, lty = 2)
    lines(envjap$obs, envjap$hi, lty = 2)
    
    pdf(paste0(fpath, "Gf-cells.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(envcells, main = "")     # pattern compatible with CSR
        dev.off()
    }
    pdf(paste0(fpath, "Gf-redwoods.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(envred, main = "")     # pattern compatible with CSR
        dev.off()
    }
}

# F function: distance from a point to the nearest event
{
    r <- seq(0, sqrt(2)/6, by = 0.001)
    Fenvjap <- envelope(as(japanesepines, "ppp"), fun = Fest,
                        nrank = 2, nsim = 99)
    Fenvred <- envelope(as(redwoodfull, "ppp"), fun = Fest, 
                        nrank = 2, nsim = 99)
    Fenvcells <- envelope(as(cells, "ppp"), fun = Fest,
                          nrank = 2, nsim = 99)
    
    pdf(paste0(fpath, "Ff-jpines.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Fenvjap, main = "")     # pattern compatible with CSR
        dev.off()
    }
    pdf(paste0(fpath, "Ff-cells.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Fenvcells, main = "")     # pattern compatible with CSR
        dev.off()
    }
    pdf(paste0(fpath, "Ff-redwoods.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Fenvred, main = "")     # pattern compatible with CSR
        dev.off()
    }
}

# K function: distance from a point to the nearest event
{
    r <- seq(0, sqrt(2)/6, by = 0.001)
    Kenvjap <- envelope(as(japanesepines, "ppp"), fun = Kest,
                        r = r, nrank = 2, nsim = 99)
    Kenvred <- envelope(as(redwoodfull, "ppp"), fun = Kest, r = r,
                        nrank = 2, nsim = 99)
    Kenvcells <- envelope(as(cells, "ppp"), fun = Kest,
                          r = r, nrank = 2, nsim = 99)
    
    pdf(paste0(fpath, "Kf-jpines.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Kenvjap, main = "")     # pattern compatible with CSR
        dev.off()
    }
    pdf(paste0(fpath, "Kf-cells.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Kenvcells, main = "")     # pattern compatible with CSR
        dev.off()
    }
    pdf(paste0(fpath, "Kf-redwoods.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Kenvred, main = "")     # pattern compatible with CSR
        dev.off()
    }
}

# remove curve from K function to give clearer representation
{
    Kjap <- Kenvjap; {
        Kjap$obs <- Kjap$obs - pi * r^2
        Kjap$theo <- Kjap$theo - pi * r^2
        Kjap$lo <- Kjap$lo - pi * r^2
        Kjap$hi <- Kjap$hi - pi * r^2
    }
    Krw <- Kenvred; {
        Krw$obs <- Krw$obs - pi * r^2
        Krw$theo <- Krw$theo - pi * r^2
        Krw$lo <- Krw$lo - pi * r^2
        Krw$hi <- Krw$hi - pi * r^2
    }
    Kcells <- Kenvcells; {
        Kcells$obs <- Kcells$obs - pi * r^2
        Kcells$theo <- Kcells$theo - pi * r^2
        Kcells$lo <- Kcells$lo - pi * r^2
        Kcells$hi <- Kcells$hi - pi * r^2
    }
    
    pdf(paste0(fpath, "Kf-jpines-flat.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Kjap, main = "")     # pattern compatible with CSR
        dev.off()
    }    
    pdf(paste0(fpath, "Kf-cells-flat.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Kcells, main = "")     # pattern compatible with CSR
        dev.off()
    }    
    pdf(paste0(fpath, "Kf-redwoods-flat.pdf")); {
        par(mar = c(4.5, 4.5,0,0.5))
        plot(Krw, main = "")     # pattern compatible with CSR
        dev.off()
    }    
}
}

####################################################################################################

# 2.5 QUADRAT TESTS 
{
quadrat.test(as(japanesepines, "ppp"), nx = 3, ny = 3)                              # X2 = 15.169, p = 0.11
quadrat.test(as(japanesepines, "ppp"), nx = 3, ny = 3, alternative = "clustered")   # X2 = 15.169, p = 0.05594



}

####################################################################################################

# TESTS APPLIED TO BAD PIXEL MAP
{
    
    # manual calculation & plotting
    {
        bpx <- bp$"141009"
        # remove lines & any other categories to be excluded
        {
            bpx <- bpx[bpx$type != "line.b", ]
        }
        bp.pp <- ppp(bpx$row, bpx$col, c(1,1996), c(1,1996))
        
        plot(bpx[,1:2], pch = 20, col = adjustcolor(Cat.cols[bpx$type], alpha = 0.5), asp = T)
        
        # tests of CSR
        {
            env.k <- envelope(bp.pp, Kest, nsim = 99, nrank = 2)
            plot(env.k); summary(env.k)
            
            env.g <- envelope(bp.pp, Gest, nsim = 99, nrank = 2)
            plot(env.g); summary(env.g)
            
            env.f <- envelope(bp.pp, Fest, nsim = 99, nrank = 2)
            plot(env.f); summary(env.f)
            
            env.h <- envelope(bp.pp, Hest, nsim = 99, nrank = 2)
            plot(env.h); summary(env.h)
        }
    }
}

####################################################################################################