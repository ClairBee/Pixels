
library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Spatial/fig/"
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", NA, "yellow", "grey", NA, "blue", "blue", "green3")

set.seed(24747)

####################################################################################################

# 7.3 PRELIMINARY ANALYSIS OF A POINT PATTERN                                                   ####

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

####################################################################################################

# 2.5 QUADRAT TESTS                                                                             ####

quadrat.test(as(japanesepines, "ppp"), nx = 3, ny = 3)                              # X2 = 15.169, p = 0.11
quadrat.test(as(japanesepines, "ppp"), nx = 3, ny = 3, alternative = "clustered")   # X2 = 15.169, p = 0.05594




####################################################################################################

# TESTS APPLIED TO BAD PIXEL MAP                                                                ####

bp <- readRDS(paste0(fpath, "bad-px.rds"))

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

####################################################################################################

# ENVELOPE FUNCTIONS                                                                            ####

# quick plotting function
plot.dist.ests <- function(bpx, dt, excl = c("line.b", "line.d"), file.id, cc = Cat.cols) {
    
    dt <- toString(dt)
    
    # plot all bad pixels
    pdf(paste0(fpath, file.id, "-", dt, "-plot.pdf")); {
        plot(bpx[[dt]][!(bpx[[dt]]$type %in% excl),1:2], 
             col = adjustcolor(cc[bpx[[dt]][!(bpx[[dt]]$type %in% excl),"type"]], alpha = 0.5), 
             pch = 20, asp = T, xlab = "", ylab = "")
        dev.off()
    }
    cat("Bad pixel plot created.")
    
    # various intensity functions
    pdf(paste0(fpath, file.id, "-", dt, "-Kenv.pdf")); {
        plot(envelope(ppp(bpx[[dt]][!(bpx[[dt]]$type %in% excl),"row"], 
                          bpx[[dt]][!(bpx[[dt]]$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Kest, nsim = 99, nrank = 2), 
             main = "")
        dev.off()
    }
    cat("K-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Fenv.pdf")); {
        plot(envelope(ppp(bpx[[dt]][!(bpx[[dt]]$type %in% excl),"row"], 
                          bpx[[dt]][!(bpx[[dt]]$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Fest, nsim = 99, nrank = 2), 
             main = "")
        dev.off()
    }
    cat("F-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Genv.pdf")); {
        plot(envelope(ppp(bpx[[dt]][!(bpx[[dt]]$type %in% excl),"row"], 
                          bpx[[dt]][!(bpx[[dt]]$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Gest, nsim = 99, nrank = 2), 
             main = "")
        dev.off()
    }
    cat("G-function plot created.")
    
    pdf(paste0(fpath, file.id, "-", dt, "-Henv.pdf")); {
        plot(envelope(ppp(bpx[[dt]][!(bpx[[dt]]$type %in% excl),"row"], 
                          bpx[[dt]][!(bpx[[dt]]$type %in% excl),"col"], 
                          c(1,1996), c(1,1996)),
                      Hest, nsim = 99, nrank = 2), 
             main = "")
        dev.off()
    }
    cat("H-function plot created.")
}

# plots of first image, excluding 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 141009, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "excl-l")

# plots of first image, including 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 141009, excl = c("line.b", "line.d"), file.id = "incl-l")

# plots of first image, including ONLY 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 141009, excl = c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "v.dim", "dim"), file.id = "only-l")

# plots of latest image, including 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 160430, excl = c("line.b", "line.d"), file.id = "incl-l")

# plots of latest image, excluding 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim"), file.id = "excl-l")

# plots of latest image, including ONLY 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 160430, excl = c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "v.dim", "dim"), file.id = "only-l")

# plots of latest image, including 'locally bright'/'locally dim' pixels
plot.dist.ests(bp, 160430, excl = c("banananana"), file.id = "even-lines")

####################################################################################################

# QUADRAT TESTS                                                                                 ####

bp.ppp <- function(bpx, dt, excl, im.dim = c(1996, 1996)) {
    
    dt <- toString(dt)
    
    ppp(bpx[[dt]][!(bpx[[dt]]$type %in% excl),"row"], 
        bpx[[dt]][!(bpx[[dt]]$type %in% excl),"col"], 
        c(1,im.dim[1]), c(1,im.dim[2])) 
}

# split by subpanels - df 31
quadrat.test(bp.ppp(bp, 160430, excl = c("line.b", "line.d")),
             xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)      # x2 = 22171, p = 0

quadrat.test(bp.ppp(bp, 160430, excl = c("line.b", "line.d", "l.bright", "l.dim")), 
             xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)      # x2 = 201.58, p = 0

quadrat.test(bp.ppp(bp, 141009, excl = c("line.b", "line.d")), 
             xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)      # x2 = 3887.9, p = 0

quadrat.test(bp.ppp(bp, 141009, excl = c("line.b", "line.d", "l.bright", "l.dim")), 
             xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)      # x2 = 85.671, p = 0.000001015


# try subpanels + minipanels (each subpanel divided vertically into 2) - df 63
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


# also test for homogeneity within subpanels?

####################################################################################################

# BAD PIXEL MAP WITH CLUSTERS & LINES REMOVED                                                   ####
# (bpx object currently taken fro 16-06-17-feature-roots.R)

px <- bpx[bpx$f.type %in% c("cl.root", "singleton"),]
excl <- c("l.bright", "l.dim")

plot(px[!(px$type %in% excl),1:2], 
     col = adjustcolor(Cat.cols[px[!(px$type %in% excl),"type"]], alpha = 0.5), 
     pch = 20, asp = T, xlab = "", ylab = "")

plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                  px[!(px$type %in% excl),"col"], 
                  c(1,1996), c(1,1996)),
              Kest, nsim = 99, nrank = 2), 
     main = "")

plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                  px[!(px$type %in% excl),"col"], 
                  c(1,1996), c(1,1996)),
              Gest, nsim = 99, nrank = 2), 
     main = "")

plot(envelope(ppp(px[!(px$type %in% excl),"row"], 
                  px[!(px$type %in% excl),"col"], 
                  c(1,1996), c(1,1996)),
              Fest, nsim = 99, nrank = 2), 
     main = "")