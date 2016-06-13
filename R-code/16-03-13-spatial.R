
library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Spatial/fig/"

bp <- readRDS(paste0(fpath, "bad-px.rds"))

####################################################################################################

# 7.3 PRELIMINARY ANALYSIS OF A POINT PATTERN                                                   ####

library(spatstat)
data(japanesepines)

plot(japanesepines, pch = 20)     # pattern compatible with CSR

# G function: distance to the nearest event
{
    r <- seq(0, sqrt(2)/6, by = 0.005)
    envjap <- envelope(as(japanesepines, "ppp"), fun = Gest,
                       r = r, nrank = 2, nsim = 99)
    envred <- envelope(as(redwoodfull, "ppp"), fun = Gest, r = r,
                       nrank = 2, nsim = 99)
    envcells <- envelope(as(cells, "ppp"), fun = Gest,
                         r = r, nrank = 2, nsim = 99)
    
    par(mfrow = c(1,3))
    plot(cells, pch = 20)
    plot(japanesepines, pch = 20)
    plot(redwoodfull, pch = 20)
    
    plot(envcells)
    plot(envjap)
    plot(envred)
    par(mfrow = c(1,1))
}

# F function: distance from a point to the nearest event
{
    r <- seq(0, sqrt(2)/6, by = 0.001)
    Fenvjap <- envelope(as(japanesepines, "ppp"), fun = Fest,
                        r = r, nrank = 2, nsim = 99)
    Fenvred <- envelope(as(redwoodfull, "ppp"), fun = Fest, r = r,
                        nrank = 2, nsim = 99)
    Fenvcells <- envelope(as(cells, "ppp"), fun = Fest,
                          r = r, nrank = 2, nsim = 99)
    
    par(mfrow = c(1,3))
        plot(Fenvcells)
        plot(Fenvjap)
        plot(Fenvred)
    par(mfrow = c(1,1))
}


