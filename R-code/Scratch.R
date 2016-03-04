
#=============================================================================================
pw.sd.b <- pixelwise.sd(b.150828)

# 1319 pixels with high SD, but only 788 classified as 'noisy'.
# remainder already classified as bright/dim. Need to establish hierarchy.
noisy.sd <- pw.sd.b[as.matrix(b.150828.bad[b.150828.bad$type == "Noisy",c(1:2)])]
variable.sd <- pw.sd.b[pw.sd.b > (6 * median(pw.sd.b))]

w <- 5
z <- hist(pw.sd.b, breaks = c(0:ceiling((max(pw.sd.b)/w))*w), 
          main = "", xlab = "Pixelwise standard deviation", ylab = "# pixels",
          ylim = c(0, 500), xlim = c(0, 500), col = "black")
hist(variable.sd, breaks = z$breaks, col = "cornflowerblue", add = T, border = "cornflowerblue")
hist(noisy.sd, breaks = z$breaks, col = "gold", add = T, border = "gold")

legend("topright", pch = 15, bty = "n",
       col = c("black", "cornflowerblue", "gold"),
       legend = c("All pixels", "SD > 6x median SD, other classification", "Classified as 'Noisy'"))

#=============================================================================================
# LOCAL NON-UNIFORMITY

# simple filter: pixels that are max/min of their neighbourhood
pw.m.b.150828 <- pixelwise.mean(b.150828)

#=============================================================================================
# median filter: pixels that differ by more than 1% from local median
system.time({
    local.med.9.9 <- as.matrix(focal(raster(pw.m.b.150828), 
                                        w = matrix(c(rep(1,81)), nrow = 9), 
                                        fun = median, na.rm = T))
})      # 158.45 elapsed; now 77.386

local.low <- which(pw.m.b.150828 < (0.99 * local.med.9.9), arr.ind = T)
local.high <- which(pw.m.b.150828 > (1.01 * local.med.9.9), arr.ind = T)

plot.ppp(ppp(local.low[,1], local.low[,2], c(1,1996), c(1,1996)))
plot.ppp(ppp(local.high[,1], local.high[,2], c(1,1996), c(1,1996)))


#=============================================================================================
# simple filter: pixels that are max/min of their neighbourhood
local.max.9.9 <- as.matrix(focal(raster(pw.m.b.150828), 
                                 w = matrix(c(rep(1,40), NA, rep(1,40)), nrow = 9), 
                                 fun = max, na.rm = T))

local.min.9.9 <- as.matrix(focal(raster(pw.m.b.150828), 
                                 w = matrix(c(rep(1,40), NA, rep(1,40)), nrow = 9), 
                                 fun = min, na.rm = T))

dark.9.9 <- which(pw.m.b.150828 < local.min.9.9, arr.ind = T)
bright.9.9 <- which(pw.m.b.150828 > local.max.9.9, arr.ind = T)

ppp.bright <- ppp(bright.9.9[,1], bright.9.9[,2], c(1,1996), c(1,1996))
plot.ppp(ppp.bright)

ppp.dark <- ppp(dark.9.9[,1], dark.9.9[,2], c(1,1996), c(1,1996))
plot.ppp(ppp.dark)
# vertical panel edges are visible: anything using value of neighbours will need to take this into account!
