
library("IO.Pixels"); library(raster)
load.images(150828, "black")
load.images(150828, "grey")
load.images(150828, "white")


# Various approaches to identification of dead pixels

#=============================================================================================
#
# approach based on manual, but without performing offset:
#   dead pixel: daily mean value is exactly 0
#   hot pixel: daily mean value is exactly 65535
#   underperforming bright pixel: value > 1.5x median
#   unerperforming dark pixel: value < 0.45x median
#   noise: sigma > 6x median sigma
#   globally non-uniform: value is > 2% from median value
#   locally non-uniform: value is > 1% from median of its 9x9 neighbours



system.time(bp <- threshold.pixels(b.150828))     # 174.03 elapsed - now 107.320

table(bp$type)

pw.m.b.150828 <- pixelwise.mean(b.150828)

pixel.image(pw.m.b.150828, title = "Pixelwise mean value, black channel 15-08-28")
points(bp[,1], bp[,2], pch = 20, cex = 0.5, col = c("black", "white", "aquamarine", "grey", "darkgreen")[bp$type])


# how does bad pixel identification differ if we carry it out on a panel-by-panel basis?
p <- panel.edges()
panel.bp <- list()

for (i in 1:32) {
    j <- ceiling(i/16)
    k <- ((i-1) %% 16)+1
    
    panel <- b.150828[c(p$x[k]:(p$x[k+1]-1)),
                      c(p$y[j]:(p$y[j+1]-1)),]
    tmp <- get.bad.pixels(panel)
    panel.bp[[i]] <- cbind(panel = paste0(p$x[k], ":", (p$x[k+1]-1),", ", p$y[j],":",(p$y[j+1]-1)),
                           tmp)
}

panel.all <- do.call("rbind", panel.bp)

# save all as .csv to create neater tables later
write.csv(panel.all, "./Other-data/Panelwise-bad-pixels.csv", row.names = F)

bp$panel.no <- apply(cbind(cut(bp$col, p$y-1), "-", cut(bp$row, p$x-1)), 1, paste, collapse = "")
write.csv(bp, "../../Other-data/All-bad-pixels.csv", row.names = F)

#=============================================================================================

bp.grey <- threshold.pixels(g.150828)
bp.grey$panel <- apply(cbind(cut(bp.grey$col, p$y-1), "-",
                             cut(bp.grey$row, p$x-1)), 1, paste, collapse = "")
write.csv(bp.grey, "./Other-data/All-bad-pixels-grey.csv", row.names = F)

bp.white <- threshold.pixels(w.150828)
bp.white$panel <- apply(rbind(cut(bp.white$col, p$y-1), "-",
                              cut(bp.white$row, p$x-1)), 1, paste, collapse = "")
write.csv(bp.white, "./Other-data/All-bad-pixels-white.csv", row.names = F)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# get panelwise values as well
p <- panel.edges()

panel.bp.white <- list()
for (i in 1:32) {
    j <- ceiling(i/16)
    k <- ((i-1) %% 16)+1
    
    panel <- w.150828[c(p$x[k]:(p$x[k+1]-1)),
                      c(p$y[j]:(p$y[j+1]-1)),]
    tmp <- threshold.pixels(panel)
    panel.bp.white[[i]] <- cbind(panel = paste0(p$x[k], ":", (p$x[k+1]-1),", ", p$y[j],":",(p$y[j+1]-1)),
                           tmp)
}
bp.white.panels <- do.call("rbind", panel.bp.white)
write.csv(bp.white.panels, "./Other-data/Panelwise-bad-pixels-white.csv", row.names = F)


panel.bp.grey <- list()
for (i in 1:32) {
    j <- ceiling(i/16)
    k <- ((i-1) %% 16)+1
    
    panel <- g.150828[c(p$x[k]:(p$x[k+1]-1)),
                      c(p$y[j]:(p$y[j+1]-1)),]
    tmp <- threshold.pixels(panel)
    panel.bp.grey[[i]] <- cbind(panel = paste0(p$x[k], ":", (p$x[k+1]-1),", ", p$y[j],":",(p$y[j+1]-1)),
                                 tmp)
}
bp.grey.panels <- do.call("rbind", panel.bp.grey)
write.csv(bp.grey.panels, "./Other-data/Panelwise-bad-pixels-grey.csv", row.names = F)

#=============================================================================================
# simple csv of each base 'bad pixel' set, ready to produce ppp plots for comparison


load.images(150828, "black")
load.images(150828, "grey")
load.images(150828, "white")

b.150828.bad <- threshold.pixels(b.150828)
g.150828.bad <- threshold.pixels(g.150828)
w.150828.bad <- threshold.pixels(w.150828)

load.images(141009, "black")
load.images(141009, "grey")
load.images(141009, "white")

b.141009.bad <- threshold.pixels(b.141009)
g.141009.bad <- threshold.pixels(g.141009)
w.141009.bad <- threshold.pixels(w.141009)


# merge into one table and save as .csv file for easier conversion to knitr
df <- rbind(cbind(src = "b.141009", b.141009.bad[,c(1:2)], type = as.character(b.141009.bad[,3])),
            cbind(src = "g.141009", g.141009.bad[,c(1:2)], type = as.character(g.141009.bad[,3])),
            cbind(src = "w.141009", w.141009.bad[,c(1:2)], type = as.character(w.141009.bad[,3])),
            cbind(src = "b.150828", b.150828.bad[,c(1:2)], type = as.character(b.150828.bad[,3])),
            cbind(src = "g.150828", g.150828.bad[,c(1:2)], type = as.character(g.150828.bad[,3])),
            cbind(src = "w.150828", w.150828.bad[,c(1:2)], type = as.character(w.150828.bad[,3])))

write.csv(df, "../../Other-data/Bad-pixels-by-threshold.csv", row.names = F)

#=============================================================================================
# Check impact of threshold on noisy pixels
assess.noise <- function(data) {
    
    # get pixelwise standard deviation & median value for all pixels
    pw.sd <- pixelwise.sd(data)
    median.sd <- median(pw.sd)
    
    l <- 1
    i <- 1
    noise <- matrix(ncol = 2, nrow = 10)
    while ((l > 0) & (i < 11)) {
        noise[i, 1] <- i
        l <- length(which(pw.sd > (i * median.sd)))
        noise[i, 2] <- l
        i <- i + 1
    }
    
    noise
}

noise.black.1 <- assess.noise(b.141009)
noise.grey.1 <- assess.noise(g.141009)
noise.white.1 <- assess.noise(w.141009)

noise.black.2 <- assess.noise(b.150828)
noise.grey.2 <- assess.noise(g.150828)
noise.white.2 <- assess.noise(w.150828)

noise <- cbind(b.141009 = noise.black.1[,2], 
               b.150828 = noise.black.2[,2],
               g.141009 = noise.grey.1[,2],
               g.150828 = noise.grey.2[,2],
               w.141009 = noise.white.1[,2],
               w.150828 = noise.white.2[,2])
write.csv(noise, "../../Other-data/Noise.csv", row.names = F)

brightness <- cbind(b.141009 = brightness.black.1[,2],
                    b.150828 = brightness.black.2[,2],
                    g.141009 = brightness.grey.1[,2],
                    g.150828 = brightness.grey.2[,2],
                    w.141009 = brightness.white.1[,2],
                    w.150828 = brightness.white.2[,2])
rownames(brightness) <- brightness.black.1[,1]
write.csv(brightness, "../../Other-data/Brightness.csv", row.names = T)


plot(noise.black[c(2:9), 1], noise.black[c(2:9), 2], type = "l", ylim = c(0,10000))
points(noise.grey[c(2:9), 1], noise.grey[c(2:9), 2], type = "l")
points(noise.white[c(2:9), 1], noise.white[c(2:9), 2], type = "l")

assess.brightness <- function(data) {
    pw.mean <- pixelwise.mean(data)
    median.value <- mean(pw.mean)
    
    b <- matrix(ncol = 2, nrow = 10)
    d <- matrix(ncol = 2, nrow = 10)
    
    for (i in 1:10) {
        b[i,1] <- (i/10) + 1
        b[i,2] <- length(which(pw.mean > (((i/10)+1) * median.value)))
        
        d[i,1] <- i/10
        d[i,2] <- length(which(pw.mean < (i/10) * median.value))
        
    }
    rbind(d,b)
} 

brightness.black.1 <- assess.brightness(b.141009)
brightness.grey.1 <- assess.brightness(g.141009)
brightness.white.1 <- assess.brightness(w.141009)

brightness.black.2 <- assess.brightness(b.150828)
brightness.grey.2 <- assess.brightness(g.150828)
brightness.white.2 <- assess.brightness(w.150828)

plot(brightness.black[c(1:9),1], brightness.black[c(1:9),2], type = "l", ylim = c(0,5000), lwd = 2,
     xlab = "Threshold as prop. of median value", ylab = "Pixels classified as 'underperforming dark'")
points(brightness.grey[c(1:9),1], brightness.grey[c(1:9),2], type = "l", col = "grey", lwd = 2)
points(brightness.white[c(1:9),1], brightness.white[c(1:9),2], type = "l", col = "gold", lwd = 2)
abline(v = 0.45, col = "red", lty = 3)


plot(brightness.black[c(11:20),1], brightness.black[c(11:20),2], type = "l", ylim = c(0,10000), lwd = 2,
     xlab = "Threshold as prop. of median value", ylab = "Pixels classified as 'underperforming bright'",
     main = "10,000 pixels = 0.025% of total image")
points(brightness.grey[c(11:20),1], brightness.grey[c(11:20),2], type = "l", col = "grey", lwd = 2)
points(brightness.white[c(11:20),1], brightness.white[c(11:20),2], type = "l", col = "gold", lwd = 2)
abline(v = 1.5, col = "red", lty = 3)


#=============================================================================================
# local non-uniformity:
system.time({
local.median.9.9 <- as.matrix(focal(raster(pw.m.b), 
                                    w = matrix(c(rep(1,81)), nrow = 9), 
                                    fun = median, na.rm = T))
})      # 158.45 elapsed

tmp <- which(pw.m.b < (0.99 * local.median.9.9), arr.ind = T)
bad.pixels <- rbind(bad.pixels, 
                    data.frame(x = tmp[,1], y = tmp[,2], type = rep("Local low", nrow(tmp))))

tmp <- which(pw.m.b > (1.01 * local.median.9.9), arr.ind = T)
bad.pixels <- rbind(bad.pixels, 
                    data.frame(x = tmp[,1], y = tmp[,2], type = rep("Local high", nrow(tmp))))


# global non-uniformity:
# temporarily excluded: picks up 1/3 of all pixels each time under current method
# tmp <- which((pw.m.b > (1.02 * median.black)), arr.ind = T)       # 1341656
# tmp <- which((pw.m.b < (0.98 * median.black)), arr.ind = T)       # 1301899


# count of pixel types before removing duplicated
# Dead     Hot    Too bright   Too dark      Noisy  Local low   Local high 
# 5        125       1130          5         1319     821413      810268 

# remove duplicates from bad pixel list - should keep first instance as long as data isn't sorted
bad.pixels <- bad.pixels[!duplicated(bad.pixels[,c(1:2)]),]


# Dead     Hot    Too bright   Too dark      Noisy  Local low   Local high 
# 5        125       1005          0          788     821407       808417 