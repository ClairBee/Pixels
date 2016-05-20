
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Notes/Standard-deviations/sd-plots/"
load.pixel.sds()
bp <- readRDS(paste0(fpath, "bad-px-maps.rds"))

cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")
cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "gold", "yellow", "grey", "violet", "green", "green3", "green", "skyblue", "lightskyblue")
fancyCat <- c("No response", "Dead", "Hot", "V. bright", "Bright", "Bright line", "Locally bright", "Slightly bright", "Screen spot", "Dim line", "Edge", "V. dim", "Dim", "Locally dim", "Slightly dim")


library(abind)      # merge subpanel arrays

####################################################################################################
# EXPLORATORY                                                                                   ####

# no point in using median-differencing: pattern is not so systematic

# try Johnson quantiles instead
# 'official' definition is > 6 sigma from median

pdf(paste0(fpath, "sd-hist-black-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"black", "160430"], xlim = c(0, 100), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"black", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"black", "160430"]) + (6 * sd(pw.sd[,,"black", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"black", "160430"]) + (6 * mad(pw.sd[,,"black", "160430"])), col = "red", lty = 2)
    dev.off()
}

pdf(paste0(fpath, "sd-hist-grey-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"grey", "160430"], xlim = c(0, 400), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"grey", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"grey", "160430"]) + (6 * sd(pw.sd[,,"grey", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"grey", "160430"]) + (6 * mad(pw.sd[,,"grey", "160430"])), col = "red", lty = 2)
    dev.off()
}  

pdf(paste0(fpath, "sd-hist-white-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"white", "160430"], xlim = c(0, 700), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"white", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"white", "160430"]) + (6 * sd(pw.sd[,,"white", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"white", "160430"]) + (6 * mad(pw.sd[,,"white", "160430"])), col = "red", lty = 2)
    dev.off()
}  

length(which(pw.sd[,,"black", "160430"] > qJohnson(0.9999, JohnsonFit(pw.sd[,,"black", "160430"]))))
length(which(pw.sd[,,"black", "160430"] > median(pw.sd[,,"black", "160430"]) + (6 * mad(pw.sd[,,"black", "160430"]))))

# Normal & Johnson Q-Q plots of SD at each power setting
{
    QQ.norm <- function(data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99),...) {
        plot(qnorm(quantiles, mean(data), sd(data)),
             quantile(data, quantiles), pch = 20, asp = T, ylab = "Observed quantile", xlab = "Normal quantile", ...)
        abline(0,1,col = "red")
        
        abline(h = quantile(data, grid.quantiles), col = "skyblue", 
               lty = 2)
        abline(v = qnorm(grid.quantiles, mean(data), sd(data)), col = "skyblue", 
               lty = 2)
    }
    
    QQ.norm(pw.sd[,,"black", "160430"], grid.quantiles = c(0.999, 0.99))     # very very not normal
    QQ.norm(pw.sd[,,"grey", "160430"], grid.quantiles = c(0.999, 0.99))      # close to normal
    QQ.norm(pw.sd[,,"white", "160430"], grid.quantiles = c(0.999, 0.99))     # close to normal
    
    Johnson.QQ(pw.sd[,,"black", "160430"], grid.quantiles = c(0.999, 0.99))
    Johnson.QQ(pw.sd[,,"grey", "160430"], grid.quantiles = c(0.999, 0.99))
    Johnson.QQ(pw.sd[,,"white", "160430"], grid.quantiles = c(0.999, 0.99))
}

# summarise standard deviations
{
    sd.summ <- apply(pw.sd, 4, 
                     function(y) data.frame(mean = apply(y, 3, mean),
                                            median = apply(y, 3, median),
                                            sd = apply(y, 3, sd),
                                            mad = apply(y, 3, mad),
                                            th = apply(y, 3, function(x) median(x) + 6 * sd(x)),
                                            noisy1 = apply(y, 3, 
                                                           function(x) length(which(x > median(x) + 6 * sd(x)))),
                                            q999 = apply(y, 3, function(x) qJohnson(0.999, JohnsonFit(x))),
                                            noisy999 = apply(y, 3, 
                                                             function(x) length(which(x > qJohnson(0.999, JohnsonFit(x))))),
                                            q9999 = apply(y, 3, function(x) qJohnson(0.9999, JohnsonFit(x))),
                                            noisy9999 = apply(y, 3, 
                                                              function(x) length(which(x > qJohnson(0.9999, JohnsonFit(x)))))))
    
    sd.summ <- lapply(sd.summ, round, 1)
    
    for (dt in names(sd.summ)) {
        write.csv(sd.summ[[dt]], paste0(fpath, "sd-summary-", dt, ".csv"), quote = F)
    }
}

# ratio of normal grey/white SD to extreme black SD
unlist(lapply(sd.summ, function(x) x$median[2] / x$th[1]))
unlist(lapply(sd.summ, function(x) x$median[3] / x$th[1]))

# compare 'official' threshold with 
unlist(lapply(sd.summ, function(x) x$q9999[2] / x$th[2]))
unlist(lapply(sd.summ, function(x) x$q9999[3] / x$th[3]))

####################################################################################################

# SD CORRELATION BETWEEN POWER SETTINGS                                                         ####

th <- list(black = c(median(pw.sd[,,"black", "160430"]) + 6 * sd(pw.sd[,,"black", "160430"]),
                     qJohnson(0.9999, JohnsonFit(pw.sd[,,"black", "160430"]))),
           grey = c(median(pw.sd[,,"grey", "160430"]) + 6 * sd(pw.sd[,,"grey", "160430"]),
                    qJohnson(0.9999, JohnsonFit(pw.sd[,,"grey", "160430"]))),
           white = c(median(pw.sd[,,"white", "160430"]) + 6 * sd(pw.sd[,,"white", "160430"]),
                     qJohnson(0.9999, JohnsonFit(pw.sd[,,"white", "160430"]))))

sp <- abind(black = subpanels(pw.sd[ , , "black", "160430"]), 
            grey = subpanels(pw.sd[ , , "grey", "160430"]),
            white = subpanels(pw.sd[ , , "white", "160430"]), along = 4)

bpx <- data.frame(bp$"160430",
                  sd.b = pw.sd[ , , "black", "160430"][as.matrix(bp$"160430"[,1:2])],
                  sd.g = pw.sd[ , , "grey", "160430"][as.matrix(bp$"160430"[,1:2])],
                  sd.w = pw.sd[ , , "white", "160430"][as.matrix(bp$"160430"[,1:2])])

# plot SD comparison for each classification
pdf(paste0(fpath, "sd-by-bad-px-type.pdf")); {
    par(mfrow = c(4, 4))
    for(ct in cat[!cat %in% c("dead", "hot", "dim", "v.dim", "line.d")]) {
        plot(bpx[bpx$type == ct, c("sd.b", "sd.g")], pch = 20, col = cat.cols[cat == ct],
             ylab = "grey", xlab = "black", main = fancyCat[cat == ct],
             xlim = c(0,5000), ylim = c(0,5000))
        abline(v = th$black, lty = 2, col = c("red", "cyan3"))
        abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
        
        plot(bpx[bpx$type == ct, c("sd.g", "sd.w")], pch = 20, col = cat.cols[cat == ct],
             ylab = "white", xlab = "grey", main = fancyCat[cat == ct],
             xlim = c(0,5000), ylim = c(0,5000))
        abline(v = th$grey, lty = 2, col = c("red", "cyan3"))
        abline(h = th$white, lty = 2, col = c("red", "cyan3"))
    }
    dev.off()
}

# get unclassified pixels, plot for eacg subpanel
sd.b.other <- pw.sd[ , , "black", "160430"]
sd.b.other[as.matrix(bp$"160430"[,1:2])] <- NA
sd.b.other <- subpanels(sd.b.other)

sd.g.other <- pw.sd[ , , "grey", "160430"]
sd.g.other[as.matrix(bp$"160430"[,1:2])] <- NA
sd.g.other <- subpanels(sd.g.other)

sd.w.other <- pw.sd[ , , "white", "160430"]
sd.w.other[as.matrix(bp$"160430"[,1:2])] <- NA
sd.w.other <- subpanels(sd.w.other)

# produce .png images of all SDs
    for(p in dimnames(sd.b.other)[[3]]) {
        png(paste0(fpath, "sd-by-power-", p, ".png"), width = 1218, height = 683, pointsize = 28); {
            par(mfrow = c(1,2), mar = c(2, 2, 3, 0.5))
            plot(sd.b.other[,,p], sd.g.other[,,p], pch = 20, xlim = c(0,1300), ylim = c(0,1300), asp = T,
                 xlab = "black", ylab = "grey", main = paste0("Panel ", p))
            abline(v = th$black, lty = 2, col = c("red", "cyan3"))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            
            plot(sd.w.other[,,p], sd.g.other[,,p], pch = 20, xlim = c(0,1300), ylim = c(0,1300), asp = T,
                 ylab = "grey", xlab = "white", main = paste0("Panel ", p))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            abline(v = th$white, lty = 2, col = c("red", "cyan3"))
            dev.off()
        }
    }

####################################################################################################
# DB CLUSTERING USING SD AT EACH POWER SETTING                                                  ####

library(dbscan)

px <- data.frame(b = c(sd.b.other[,,"U13"]), g = c(sd.g.other[,,"U13"]), w = c(sd.w.other[,,"U13"]))
px <- as.matrix(px[!is.na(px$b),])

# test plots, with sample
{
    zz <- px[sample(1:nrow(px), 1000),]
    
    hh <- kNNdist(zz, k = 4)
    
    kk <- 7
    plot(sort(kNNdist(zz[,1:2], k = kk)), type = "l")
    lines(sort(kNNdist(zz[,2:3], k = kk)), col = "blue")
    lines(sort(kNNdist(zz[,c(1, 3)], k = kk)), col = "red")
    
    plot(sort(kNNdist(zz[,1:2], k = 4)), type = "l"); lapply(c(5:30), function(x) lines(sort(sample(kNNdist(zz[,1:2], k = x), 4000)), col = x))
    abline(h = 4, col = "gold")
    
    plot(sort(kNNdist(zz[,2:3], k = 4)), type = "l"); lapply(c(5:30), function(x) lines(sort(sample(kNNdist(zz[,2:3], k = x), 4000)), col = x))
    abline(h = 15, col = "gold")
    
    plot(sort(kNNdist(zz[,c(1,3)], k = 4)), type = "l"); lapply(c(5:30), function(x) lines(sort(sample(kNNdist(zz[,c(1,3)], k = x), 4000)), col = x))
    abline(h = 5, col = "gold")
}


plot(sort(kNNdist(px[,1:2], k = 4)), type = "l", ylim = c(0, 10))
abline(h = 0.5, col = "gold")

plot(px[,1:2], col = dbscan(px[,1:2], eps = 4, minPts = 4)$cluster + 1, pch = 20, main = "Black vs grey; Epsilon = 4", xlim = c(0,200))

plot(px[,2:3], col = dbscan(px[, 2:3], eps = 15, minPts = 4)$cluster + 1, main = "Grey vs white; Epsilon = 15", pch = 20)

plot(px[, c(1,3)], col = dbscan(px[, c(1,3)], eps = 5, minPts = 4)$cluster + 1, main = "Black vs white; Epsilon = 5", pch = 20)
