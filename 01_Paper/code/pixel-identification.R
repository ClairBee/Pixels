
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()
pw.sd <- load.pixel.sds()

sc <- apply(pw.m, 4, shading.corrected)

# need to get remaining pixelwise SDs

####################################################################################################

# METHOD FROM MANUAL                                                                            ####

# get offset-corrected images
os <- pw.m[,,c("grey", "white"), ] - pw.m[,,c("black", "black"), ] # offset-corrected images

# signal sensitivity
ss <- list(up.b = which(os[,,"grey", "160430"] > 1.5 * median(os[,,"grey", "160430"], na.rm = T) |
                            os[,,"white", "160430"] > 1.5 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T),
           n.g = which(os[,,"grey", "160430"] < 10000 | os[,,"white", "160430"] < 15000, arr.ind = T),
           up.d = which(os[,,"grey", "160430"] < 0.45 * median(os[,,"grey", "160430"], na.rm = T) |
                            os[,,"white", "160430"] < 0.45 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T))
ss[sapply(ss, length) == 0] <- NULL
ss <- rbind.fill(lapply(names(ss), function(tt) data.frame(ss[[tt]], type = tt)))

# noise


# uniformity
sc <- shading.corrected(pw.m[,,,"160430"])
hist(sc, breaks = "fd", ylim = c(0,30))
abline(v = median(sc, na.rm = T) * c(0.98, 1, 1.02), lty = c(2,1,2), col = "red")

sc.md <- r2m(focal(m2r(sc), matrix(rep(1, 81), nrow = 9), fun = median))
hist(sc.md, breaks = "fd", ylim = c(0,30))
abline(v = median(sc.md, na.rm = T) * c(0.99, 1, 1.01), lty = c(2,1,2), col = "red")

uu <- list(gl.u = which(sc > 1.02 * median(sc, na.rm = T) | sc < 0.98 * median(sc, na.rm = T), arr.ind = T),
           loc.u = which(sc > 1.01 * sc.md | sc < 0.99 * sc.md, arr.ind = T))

sapply(uu, length)
pixel.plot(uu$loc.u, cex = 0.4, col = "green3")
pixel.plot(uu$gl.u, cex = 0.4)

zz <- apply(sc, 3, function(im) r2m(focal(m2r(im), matrix(rep(1, 81), nrow = 9), fun = median)))

saveRDS(zz, "./02_Objects/sc-med-diffs/sc-md-9x9.rds")

####################################################################################################

# MEDIAN-DIFFERENCING OVER INDIVIDUAL IMAGES                                                    ####

####################################################################################################

# MEDIAN-DIFFERENCING OVER SHADING CORRECTION                                                   ####

####################################################################################################

# COMPARE RESULTS                                                                               ####

# check for constant offset in any bad pixels
{
    
}

# check correlation between pixelwise SD and median difference (per image/power)
{
    
}