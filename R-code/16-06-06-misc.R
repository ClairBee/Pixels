

library("IO.Pixels"); library("CB.Misc")

d.141009 <- load.daily("141009")
bpx <- readRDS("./Notes/Standard-deviations/bad-px-maps-with-noise.rds")
fpath.fig <- "./Notes/Standard-deviations/fig/"

####################################################################################################

# RTS/'POPCORN' NOISE                                                                           ####
{
    # divide transects of noisy pixels at midpoint
    npx <- apply(d.141009[,,,"black"], 3, function(x) x[as.matrix(bpx$"141009"[!is.na(bpx$"141009"$noisy),1:2])])
    
    zz <- data.frame(mid.range = apply(npx, 1, function(x) (max(x) + min(x)) / 2),
                     sd.above = apply(npx, 1, function(x) sd(x[which(x > (max(x) + min(x)) / 2)])),
                     sd.below = apply(npx, 1, function(x) sd(x[which(x < (max(x) + min(x)) / 2)])),
                     max.above = apply(npx, 1, function(x) max(x[which(x > (max(x) + min(x)) / 2)])),
                     min.above = apply(npx, 1, function(x) min(x[which(x > (max(x) + min(x)) / 2)])),
                     max.below = apply(npx, 1, function(x) max(x[which(x < (max(x) + min(x)) / 2)])),
                     min.below = apply(npx, 1, function(x) min(x[which(x < (max(x) + min(x)) / 2)])),
                     total.sd = apply(npx, 1, sd))
    
    zz$mid.above <- ((zz$max.above + zz$min.above) / 2) - zz$mid.range
    zz$mid.below <- zz$mid.range - ((zz$max.below + zz$min.below) / 2)
    
    # get SD of pixels above & below midpoint
    zz$sd.diff <- zz$sd.above - zz$sd.below
    
    qq <- npx[which(zz$sd.above < 20 & zz$sd.below < 20),]
    
    zz[which(zz$sd.above < 10 & zz$sd.below < 10),]
    
    pdf(paste0(fpath.fig, "tmp.pdf"))
    par(mfrow = c(5,4), mar = c(2, 2, 1, 1))
    for(i in 1:nrow(qq)) {
        o.plot(qq[i,])
    }
    dev.off()
}
# possibly a number of points with RTS noise? However, hard to tell in a sequence this short.

# look for RTS noise in sequence of 100 images
{
    bl <- array(dim = c(1996, 1996, 100))
    for (i in 1:100) {
        bl[,,i] <- readTIFF(paste0("./Other-data/Other-images/Block-images/block_", i, ".tif"), as.is = T)
    }
    
    zz <- abind(mid.range = array(apply(bl, c(1:2), function(x) (max(x) + min(x)) / 2), dim = c(1996, 1996)),
                sd.above = array(apply(bl, c(1:2), function(x) sd(x[which(x > (max(x) + min(x)) / 2)])), dim = c(1996, 1996)),
                sd.below = array(apply(bl, c(1:2), function(x) sd(x[which(x < (max(x) + min(x)) / 2)])), dim = c(1996, 1996)),
                along = 3)
    
    qq <- which(zz[,,"sd.above"] < 100 & zz[,,"sd.below"] < 100, arr.ind = T)       
    
    for (i in 1:nrow(qq)) {
        o.plot(bl[qq[i,1], qq[i,2],])
    }
}
# no evidence of RTS in images with xray source

####################################################################################################

# PER-PIXEL PHOTON TRANSFER CURVES                                                              ####

{
    jay.load <- function(lvl) {
        
        lvl <- toString(lvl)
        jpath = paste0("/home/clair/Documents/Pixels/Other-data/Other-images/line_investigation/", lvl, "ua/")
        
        jay.tiffs <- list.files(jpath, pattern = "\\.tif$")
        
        # create array to hold loaded data
        ims <- array(dim = c(1996, 1996, length(jay.tiffs)))
        
        for (i in 1:length(jay.tiffs)) {
            tmp <- readTIFF(paste0(jpath, jay.tiffs[i]), as.is = T)
            
            # transpose & rotate to get image right way up
            ims[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
        }
        ims
    }
    
    lapply(c(20, 40, 60, 80, 100), jay.load)
    
    m <- abind(m.ua20, m.ua40, m.ua60, m.ua80, m.ua100, along = 3, 
               new.names = list(NULL, NULL, c("ua20", "ua40", "ua60", "ua80", "ua100")))
    s <- abind(sd.ua20, sd.ua40, sd.ua60, sd.ua80, sd.ua100, along = 3, 
               new.names = list(NULL, NULL, c("ua20", "ua40", "ua60", "ua80", "ua100")))
    
    ua <- abind(m, s, along = 4, 
                new.names = list(NULL, NULL, c("ua20", "ua40", "ua60", "ua80", "ua100"), c("mean", "sd")))}
ua <- readRDS("/home/clair/Documents/Pixels/Other-data/Other-images/ua-images.rds")

# load all images
ua <- abind(lapply(c(20, 40, 60, 80, 100), jay.load), along = 4,
            new.names = list(NULL, NULL, NULL, c("ua20", "ua40", "ua60", "ua80", "ua100")))

{
    # get per-pixel difference between ua20 and uaXX
    ua.diffs.40 <- ua[,,,"ua40"] - ua[,,,"ua20"]
    ua.diffs.60 <- ua[,,,"ua60"] - ua[,,,"ua20"]
    ua.diffs.80 <- ua[,,,"ua80"] - ua[,,,"ua20"]
    ua.diffs.100 <- ua[,,,"ua100"] - ua[,,,"ua20"]
    
    # get mean & sd of per-pixel differences
    zz <- abind(pixelwise.mean(ua.diffs.40),
                pixelwise.mean(ua.diffs.60),
                pixelwise.mean(ua.diffs.80),
                pixelwise.mean(ua.diffs.100), along = 3)
    
    qq <- abind(pixelwise.sd(ua.diffs.40),
                pixelwise.sd(ua.diffs.60),
                pixelwise.sd(ua.diffs.80),
                pixelwise.sd(ua.diffs.100), along = 3)
    
    hh <- abind(zz, qq, along = 4, new.names = list(NULL, NULL, c("ua40", "ua60", "ua80", "ua100"), c("mean", "sd")))
    
    # fit line of mean vs var for each pixel
    ll <- apply(hh, 1:2, function(x) line(x[,"mean"], x[,"sd"]^2)$coef[2])
    
    pixel.image(ll)
    
    s.hist(ll)
    
    ss <- sample(1996^2, 100000)
    
    plot(ll[ss], pch = 20, ylim = c(-50, 50))
    abline(h = mean(ll, na.rm = T) + c(2, -2) * sd(ll, na.rm = T), col = "red")
    
    plot(hh[,,"ua40", "mean"][ss], hh[,,"ua40", "sd"][ss]^2, pch = 20, xlim = c(0,50000), ylim = c(0,100000))
    points(hh[,,"ua60", "mean"][ss], hh[,,"ua60", "sd"][ss]^2, pch = 20, col = adjustcolor("blue", alpha = 0.5))
    points(hh[,,"ua80", "mean"][ss], hh[,,"ua80", "sd"][ss]^2, pch = 20, col = adjustcolor("purple", alpha = 0.5))
    points(hh[,,"ua100", "mean"][ss], hh[,,"ua100", "sd"][ss]^2, pch = 20, col = adjustcolor("green3", alpha = 0.5))
}

# get 100x100 subarray
{
    cc <- get.focus(data.frame(x = 959, y = 1200), surround = 50)
    
    sa <- apply(ua, 3:4, "[", cc)
    
    sa.m <- c(apply(sa[,, "ua40"] - sa[,, "ua20"], 2, mean), 
              apply(sa[,, "ua60"] - sa[,, "ua20"], 2, mean),
              apply(sa[,, "ua80"] - sa[,, "ua20"], 2, mean),
              apply(sa[,, "ua100"] - sa[,, "ua20"], 2, mean))
    
    sa.sd <- c(apply(sa[,, "ua40"] - sa[,, "ua20"], 2, sd), 
               apply(sa[,, "ua60"] - sa[,, "ua20"], 2, sd),
               apply(sa[,, "ua80"] - sa[,, "ua20"], 2, sd),
               apply(sa[,, "ua100"] - sa[,, "ua20"], 2, sd))
    
    lm <- lm(log(sa.sd) ~ log(sa.m))
    plot(log(sa.m), log(sa.sd), pch = 20)
    abline(line(log(sa.m), log(sa.sd)), col = "cyan3")
    abline(coef(lm), col = "red")
    
    lm2 <- lm(log(sa.sd) ~ poly(log(sa.m), 2))
    lines(log(sa.m), lm2$fitted, col = "red", pch = 20)   
    lines(log(sa.m), lm(log(sa.sd) ~ poly(log(sa.m), 3))$fitted, col = "blue")   
    
    # need more points. Doesn't cover enough of dynamic range.
}

