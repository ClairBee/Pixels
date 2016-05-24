
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Notes/Standard-deviations/sd-plots/"
fpath.fig <- "./Notes/Standard-deviations/fig/"

load.pixel.sds()
bp <- readRDS("./Other-data/bad-px-maps.rds")
bpx <- readRDS("./Notes/Standard-deviations/bad-px-maps-with-noise.rds")

Cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "gold", "yellow", "grey", "violet", "green", "green3", "green", "skyblue", "lightskyblue")
fancyCat <- c("No response", "Dead", "Hot", "V. bright", "Bright", "Bright line", "Locally bright", "Slightly bright", "Screen spot", "Dim line", "Edge", "V. dim", "Dim", "Locally dim", "Slightly dim")
headerCat <- gsub("[ ]", "", gsub("[.]", "", Cat))
longCat <- c("normal", Cat, "noisy")

load.pixel.means()

####################################################################################################
# EXPLORATORY                                                                                   ####

# no point in using median-differencing: pattern is not so systematic

# try Johnson quantiles instead?
# 'official' definition is > 6 sigma from median



# plot per-colour SD relation: bad px and unclassified px
png(paste0(fpath.fig, "sd-plot-bad-px.png"), width = 1200, height = 600, pointsize = 28); {
    par(mar = c(2, 2, 1, 1), mfrow = c(1, 2))
    plot(pw.sd[,,"black", "160430"][as.matrix(bp$"160430"[,1:2])], 
         pw.sd[,,"grey", "160430"][as.matrix(bp$"160430"[,1:2])],
         pch = 20, col = Cat.cols[bp$"160430"$type], xlim = c(0,5000), ylim = c(0,5000))
    abline(v = th$black, lty = 2)
    abline(h = th$grey, lty = 2)
    text(4800,0,"black"); text(0, 4500, "grey", srt = 90)
    
    plot(pw.sd[,,"white", "160430"][as.matrix(bp$"160430"[,1:2])], 
         pw.sd[,,"grey", "160430"][as.matrix(bp$"160430"[,1:2])],
         pch = 20, col = Cat.cols[bp$"160430"$type], xlim = c(0,5000), ylim = c(0,5000))
    abline(v = th$white, lty = 2)
    abline(h = th$grey, lty = 2)
    text(4800,0,"white"); text(0, 4500, "grey", srt = 90)
    
    dev.off()
}

sdn.b <- pw.sd[,,"black", "160430"]; sdn.b[as.matrix(bp$"160430"[,1:2])] <- NA
sdn.b <- pw.sd[,,"grey", "160430"]; sdn.b[as.matrix(bp$"160430"[,1:2])] <- NA
sdn.w <- pw.sd[,,"white", "160430"]; sdn.w[as.matrix(bp$"160430"[,1:2])] <- NA

png(paste0(fpath.fig, "sd-plot-unc-px.png"), width = 1200, height = 600, pointsize = 28); {
    par(mar = c(2, 2, 1, 1), mfrow = c(1, 2))
    plot(sdn.b, sdn.g, pch = 20, xlim = c(0,5000), ylim = c(0,5000))
    abline(v = th$black, lty = 2)
    abline(h = th$grey, lty = 2)
    text(4800,0,"black"); text(0, 4500, "grey", srt = 90)
    
    plot(sdn.w, sdn.g, pch = 20, xlim = c(0,5000), ylim = c(0,5000))
    abline(v = th$white, lty = 2)
    abline(h = th$grey, lty = 2)
    text(4800,0,"white"); text(0, 4500, "grey", srt = 90)
    
    dev.off()
}


#===================================================================================================


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

# plot SD comparison for each classifiCation
pdf(paste0(fpath, "sd-by-bad-px-type.pdf")); {
    par(mfrow = c(4, 4))
    for(ct in Cat[!Cat %in% c("dead", "hot", "dim", "v.dim", "line.d")]) {
        plot(bpx[bpx$type == ct, c("sd.b", "sd.g")], pch = 20, col = Cat.cols[Cat == ct],
             ylab = "grey", xlab = "black", main = fancyCat[Cat == ct],
             xlim = c(0,5000), ylim = c(0,5000))
        abline(v = th$black, lty = 2, col = c("red", "cyan3"))
        abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
        
        plot(bpx[bpx$type == ct, c("sd.g", "sd.w")], pch = 20, col = Cat.cols[Cat == ct],
             ylab = "white", xlab = "grey", main = fancyCat[Cat == ct],
             xlim = c(0,5000), ylim = c(0,5000))
        abline(v = th$grey, lty = 2, col = c("red", "cyan3"))
        abline(h = th$white, lty = 2, col = c("red", "cyan3"))
    }
    dev.off()
}

# get unclassified pixels, plot for each subpanel
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
                 xlab = "black", ylab = "grey", main = "")
            abline(v = th$black, lty = 2, col = c("red", "cyan3"))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            text(c(1200), c(00), "black"); text(00, 1200, "grey", srt = 90)
            
            plot(sd.w.other[,,p], sd.g.other[,,p], pch = 20, xlim = c(0,1300), ylim = c(0,1300), asp = T,
                 ylab = "grey", xlab = "white", main = "")
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            abline(v = th$white, lty = 2, col = c("red", "cyan3"))
            text(c(1200), c(00), "white"); text(00, 1200, "grey", srt = 90)
            dev.off()
        }
    }

# .png images of all SDs for bad pixels already identified
    for (t in Cat) {
        png(paste0(fpath, "sd-by-type-", headerCat[Cat == t], ".png"), width = 1218, height = 683, pointsize = 28); {
            par(mfrow = c(1,2), mar = c(2, 2, 3, 0.5))
            
            plot(bpx[bpx$type == t, c("sd.b", "sd.g")], col = Cat.cols[Cat == t], pch = 20, asp = T, xlim = c(0,5000), ylim = c(0,5000),
                 xlab = "black", ylab = "grey", main = fancyCat[Cat == t])
            abline(v = th$black, lty = 2, col = c("red", "cyan3"))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            
            plot(bpx[bpx$type == t, c("sd.w", "sd.g")], col = Cat.cols[Cat == t], pch = 20, asp = T, xlim = c(0,5000), ylim = c(0,5000),
                 xlab = "white", ylab = "grey", main = fancyCat[Cat == t])
            abline(v = th$white, lty = 2, col = c("red", "cyan3"))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            
            dev.off()
        }
    }

# .png images of all SDs for bad pixels already identified - per subpanel
{
    bpx <- data.frame(bpx, 
                      sp = apply(cbind(c("L", "U")[findInterval(bpx$row, panel.edges()$y - 0.5)],
                                       c(1:32)[findInterval(bpx$col, panel.edges()$x - 0.5)]), 1, paste, collapse = ""),
                      stringsAsFactors = F)
    
    for(p in apply(cbind(c(rep("U", 16), rep("L", 16)), 
                         rep(c(1:16), 2)), 1, paste, collapse = "")) {
        png(paste0(fpath, "sd-by-bpx-sp-", p, ".png"), width = 1218, height = 683, pointsize = 28); {
            par(mfrow = c(1,2), mar = c(2, 2, 3, 0.5))
            
            plot(bpx[bpx$sp == p, c("sd.b", "sd.g")], col = Cat.cols[bpx$type[bpx$sp == p]], pch = 20, asp = T, xlim = c(0,1300), ylim = c(0,1300),
                 xlab = "black", ylab = "grey", main = "")
            abline(v = th$black, lty = 2, col = c("red", "cyan3"))
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            text(c(1200), c(00), "black"); text(00, 1200, "grey", srt = 90)
            
            plot(bpx[bpx$sp == p, c("sd.w", "sd.g")], col = Cat.cols[bpx$type[bpx$sp == p]], pch = 20, asp = T, xlim = c(0,1300), ylim = c(0,1300),
                 ylab = "grey", xlab = "white", main = "")
            abline(h = th$grey, lty = 2, col = c("red", "cyan3"))
            abline(v = th$white, lty = 2, col = c("red", "cyan3"))
            text(c(1200), c(00), "white"); text(00, 1200, "grey", srt = 90)
            
            dev.off()
        }
    }
}

# px-legend
{
    pdf(paste0(fpath, "bpx-legend-1.pdf")); {
        plot.new()
        legend("center", pch = 20, col = c("black", Cat.cols[1:3]), legend = c("Unclassified", fancyCat[1:3]), horiz = T, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bpx-legend-1.pdf"))
    }
    pdf(paste0(fpath, "bpx-legend-2.pdf")); {
        plot.new()
        legend("center", pch = 20, col = c(Cat.cols[4:7]), legend = c(fancyCat[4:7]), horiz = T, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bpx-legend-2.pdf"))
    }
    pdf(paste0(fpath, "bpx-legend-3.pdf")); {
        plot.new()
        legend("center", pch = 20, col = Cat.cols[8:11], legend = fancyCat[8:11], horiz = T, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bpx-legend-3.pdf"))
    }
    pdf(paste0(fpath, "bpx-legend-4.pdf")); {
        plot.new()
        legend("center", pch = 20, col = Cat.cols[12:15], legend = fancyCat[12:15], horiz = T, bty = "n")
        
        dev.off()
        crop.pdf(paste0(fpath, "bpx-legend-4.pdf"))
    }
}

# overplot black vs grey vs white SD for unclassified pixels
{
    zz <- list("black" = join.panels(sd.b.other),
               "grey" = join.panels(sd.g.other),
               "white" = join.panels(sd.w.other))
    
    png(paste0(fpath, "sd-black-vs-white.png"), width = 1218, height = 1218, pointsize = 28); {
        plot(c(zz$white), c(zz$grey), pch = 20, col = "grey", xlim = c(0,1500), ylim = c(0,1500), asp = T,
             main = "", xlab = "", ylab = "")
        points(c(zz$black), c(zz$grey), pch = 20)
        
        abline(v = lapply(th, "[", 1), col = c("red", "cyan3", "gold"), lwd = 2)
        abline(h = lapply(th, "[", 1), col = c("red", "cyan3", "gold"), lwd = 2)
        
        dev.off()
    }
}

####################################################################################################

# BLACK SD VS GREY/WHITE SD                                                                     ####

# produce tables showing mean/median/sd/threshold/n identified.
# get ratio of black threshold to grey threshold
# get actual quantiles & compare to 

# effect on shading correction
{
    # get median value & upper limit by SD threshold
    r.b <- median(pw.m[,,"black", "160430"]) + (c(-1,0,1) * (median(pw.sd[,,"black","160430"]) + 6 * sd(pw.sd[,,"black","160430"])))
    r.g <- median(pw.m[,,"grey", "160430"]) + (c(-1,0,1) * (median(pw.sd[,,"grey","160430"]) + 6 * sd(pw.sd[,,"grey","160430"])))
    r.w <- median(pw.m[,,"white", "160430"]) + (c(-1,0,1) * (median(pw.sd[,,"white","160430"]) + 6 * sd(pw.sd[,,"white","160430"])))
    
    res <- array(dim = c(3, 3, 3), dimnames = list(c("b+", "b", "b-"), c("g+", "g", "g-"), c("w+", "w", "w-")))
    baseline <- r.g[2] - 60000 * (r.g[2] - r.b[2]) / (r.w[2] - r.b[2])   # use all medians
    
    # G - (G - B / W - B * 60000)
    for (w in 1:3) {
        for (g in 1:3) {
            for (b in 1:3) {
                res[b, g, w] <- r.g[g] - 60000 * (r.g[g] - r.b[b]) / (r.w[w] - r.b[b])
            }
        }
    }
    
    res.tab <- cbind(res[,,"w+"], res[,,"w"], res[,,"w-"])
    res.tab <- round(res.tab - baseline, 0)
    
    write.csv(res.tab, paste0(fpath.fig, "shading-corr-effect.csv"), quote = F)
    
    # summarise per power setting
    {
        round(apply(res, 2:3, mean), 0)
        round(apply(res, c(1,3), mean), 0)
        round(apply(res, 1:2, mean), 0)
        
        #           effect of b                    effect of g                 effect of w     
        #
        #          w+     w    w-                 w+    w    w-               g+    g    g-
        #   g+   -117   150   409           b+  -369  -93   174         b+    53  -96  -245
        #   g    -275     0   267           b   -275    0   267         b    147   -3  -153
        #   g-   -432  -150   124           b-  -181   93   359         b-   242   90   -61
        
        apply(cbind(res[,,"w+"], res[,,"w"], res[,,"w-"]), 2, function(x) round(max(x) - min(x),0))  # mean black effect
        apply(rbind(res[,,"w+"], res[,,"w"], res[,,"w-"]), 1, function(x) round(max(x) - min(x),0))  # mean grey effect
        
        apply(cbind(res[,"g+",], res[,"g",], res[,"g-",]), 2, function(x) round(max(x) - min(x),0))  # mean 
        apply(rbind(res[,"g+",], res[,"g",], res[,"g-",]), 1, function(x) round(max(x) - min(x),0))
        
        
        apply(cbind(res[1,,], res[2,,], res[3,,]), 2, function(x) round(max(x) - min(x),0))          # white effect
        
        apply(apply(res, 1, range), 2, function(x) max(x) - min(x))
        #              b+         b        b-
        #    [1,] -828.35127 -735.7367 -642.7060
        #    [2,]   12.74402  105.8993  199.4605
        #           841.0953  841.6360  842.1665    # range
        
        apply(apply(res, 2, range), 2, function(x) max(x) - min(x))
        #               g+          g         g-
        #    [1,] -515.3601 -671.85571 -828.35127
        #    [2,]  199.4605   55.73105  -87.99838
        #          714.8206   727.5868   740.3529    # range
        
        apply(apply(res, 3, range), 2, function(x) max(x) - min(x))
        #               w+          w        w-
        #    [1,] -828.3513 -545.19956 -270.5053
        #    [2,] -325.2416  -58.89546  199.4605
        #          503.1097   486.3041  469.9658    # range
    }
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

df <- data.frame(medb = apply(pw.sd[,,"black",], 3, median, na.rm = T),
                 medg = apply(pw.sd[,,"grey",], 3, median, na.rm = T),
                 medw = apply(pw.sd[,,"white",], 3, median, na.rm = T),
                 sigb = apply(pw.sd[,,"black",], 3, sd, na.rm = T),
                 sigg = apply(pw.sd[,,"grey",], 3, sd, na.rm = T),
                 sigw = apply(pw.sd[,,"white",], 3, sd, na.rm = T))

df$thb <- df$medb + 6 * df$sigb
df$thg <- df$medg + 6 * df$sigg
df$thw <- df$medw + 6 * df$sigw

df$thbmedg <- df$thb / df$medg
df$thbmedw <- df$thb / df$medw

df$nbl <- unlist(lapply(rownames(df), function(x) length(which(pw.sd[,,"black",x] > df$thb[rownames(df) == x]))))
df$ng <- unlist(lapply(rownames(df), function(x) length(which(pw.sd[,,"grey",x] > df$thg[rownames(df) == x]))))
df$nw <- unlist(lapply(rownames(df), function(x) length(which(pw.sd[,,"white",x] > df$thw[rownames(df) == x]))))

df$clb <- df$medb - (1.96 * df$sigb)
df$cub <- df$medb + (1.96 * df$sigb) 

df$clg <- df$medg - (1.96 * df$sigg)
df$cug <- df$medg + (1.96 * df$sigg) 

df$clw <- df$medw - (1.96 * df$sigw)
df$cuw <- df$medw + (1.96 * df$sigw)

df <- round(df, 1)
df$dt <- sapply(rownames(df), fancy.date)

write.csv(df, paste0(fpath.fig, "sd-thresholds.csv"), quote = F)

# clear difference across subpanels in behaviour in black image
# check in most recent acquisition
{
    zz <- array(unlist(lapply(lapply(alply(pw.sd, 4, 
                                           apply, 3, subpanels, .dims = T),
                                     array, dim = c(128, 1024, 32, 3)),
                              apply, 3:4, sd, na.rm = T)), 
                dim = c(32, 3, 12), dimnames = list(NULL, dimnames(pw.sd)[[3]], dimnames(pw.sd)[[4]]))

    pdf(paste0(fpath.fig, "sd-sigma-per-subpanel.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            matplot(zz[1:16,,i], type = "l", ylim = c(0,80), lty = 1)
            matplot(zz[17:32,,i], type = "l", add = T, lty = 2)
            text(14, 75, dimnames(zz)[[3]][i])
        }
        dev.off()
    }

}


####################################################################################################

# CLASSIFICATION USING EXISTING THRESHOLDS                                                      ####

th <- apply(pw.sd, 3, apply, 3, function(x) median(x) + 6 * sd(x))

# are same px identified as noisy in grey & white images? Or is there a difference?

# identify all noisy pixels in each image
npx <- lapply(apply(pw.sd, 4, apply, 3, function(x) which(x > (median(x) + 6 * sd(x)), arr.ind = T)),
              function(x) merge(merge(data.frame(x$black, src.b = T),
                                      data.frame(x$grey, src.g = T), by = c(1:2), all = T),
                                data.frame(x$white, src.w = T), by = c(1:2), all = T))

# list all pixels identified as noisy
{
    noisy <- rbind.fill(lapply(npx, "[", 1:2))
    noisy <- noisy[!dupliCated(noisy),]
    
    n.summ <- lapply(npx, function(x) count(x[,3:5]))
    
    zz <- rmerge.df.list(noisy, npx, by = c(1:2), all = T)
    
    qq <- rmerge.df.list(n.summ[[1]], n.summ[2:12], by = c(1:3), all = T)
    qq <- qq[rev(order(rowMeans(qq, na.rm = T))),]
    
    qq <- setNames(data.frame(gsub("NA", "-", apply(cbind(rep("b", 7)[qq$src.b],
                                                          rep("g", 7)[qq$src.g], 
                                                          rep("w", 7)[qq$src.w]), 1, paste, collapse = "")),
                              qq[,4:15]),
                   c("Found", sapply(names(df.list), fancy.date)))
    
    qq[is.na(qq)] <- "-"
    
    write.csv(qq, paste0(fpath.fig, "tbl-noisy-px-by-col.csv"), quote = F, row.names = F)
}

# not a lot of crossover between noisy pixels in black & white/grey images; also not masses between white & grey.

# merge noisy pixels with bad pixel maps

bpx <- list()
for (dt in dimnames(pw.sd)[[4]]) {
    tmp <- rbind(data.frame(npx[[dt]][is.na(npx[[dt]]$src.b),1:2], noisy = "gw"),
                data.frame(npx[[dt]][,1:2], noisy = "b"))
    tmp <- tmp[!dupliCated(tmp[,1:2]),]
    
    bpx[[dt]] <- merge(bp[[dt]], tmp, by = c(1:2), all = T)
    levels(bpx[[dt]]$type) <- c(levels(bpx[[dt]]$type), "noisy")
    bpx[[dt]]$type[is.na(bpx[[dt]]$type)] <- "noisy"
}

saveRDS(bpx, "./Notes/Standard-deviations/bad-px-maps-with-noise.rds")

# state transitions including noise
{
    tr <- list()
    
    for (i in 1:(length(bpx) - 1)) {
        tr[[i]] <- table("From" = ordered(longCat[c(bpx2im(bpx[[i]])) + 1], levels = longCat),
                         "To" = ordered(longCat[c(bpx2im(bpx[[i+1]])) + 1], levels = longCat))
        names(tr)[[i]] <- paste(names(bpx)[c(i,(i+1))], collapse = "-")
    }
}

# mean transition rates
tr <- array(unlist(tr), dim = c(17, 17, 11), dimnames = list(longCat, c("Normal", fancyCat, "Noisy"), names(tr)))

write.csv(prep.csv(apply(tr, 1:2, mean), dp = 0)[-11, -11],
          paste0(fpath.fig, "transition-px-mean.csv"), quote = F)

write.csv(prep.csv(apply(tr, 1:2, sd), dp = 0)[-11, -11],
          paste0(fpath.fig, "transition-px-sd.csv"), quote = F)

write.csv(prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, mean, na.rm = T) * 100)[-11, -11],
          paste0(fpath.fig, "transition-prop-mean.csv"), quote = F)

write.csv(prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, sd, na.rm = T) * 100)[-11, -11],
          paste0(fpath.fig, "transition-prop-sd.csv"), quote = F)

# pixels identified as each type
vv <- array(unlist(lapply(lapply(bpx, "[", 4:3), table, useNA = "ifany")), dim = c(3, 16, 12),
            dimnames = list(c("gw", "b", "-"), longCat[-1], names(bpx)))

write.csv(prep.csv(apply(vv, 1:2, mean), dp = 0)[, -10],
          paste0(fpath.fig, "mean-bad-px-types.csv"), quote = F)

write.csv(prep.csv(apply(vv, 1:2, sd), dp = 0)[, -10],
          paste0(fpath.fig, "sd-bad-px-types.csv"), quote = F)

####################################################################################################

# NOISE VS LOCAL BRIGHTNESS                                                                     ####



####################################################################################################

# STABILITY OF NOISY PX IN UA IMAGES                                                            ####

# import pixelwise SDs
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
    
    # get pixelwise mean
    m <- apply(ims, c(1, 2), sd)
    
    # assign to new object
    assign(paste0("ua", lvl), m, envir = .GlobalEnv)
}

j.sd <- lapply(c(20, 40, 60, 80, 100), jay.load)
j.mean <- lapply(c(20, 40, 60, 80, 100), jay.load)

th <- data.frame(med = unlist(lapply(j, median)),
                 sd = unlist(lapply(j, sd)),
                 th = unlist(lapply(j, function(x) median(x) + 6 * sd(x))))

npx <- lapply(j.sd, function(x) which(x > (median(x) + 6 * sd(x)), arr.ind = T))

zz <- rbind.fill(lapply(npx, as.data.frame))
zz <- zz[!dupliCated(zz),]

hh <- merge(merge(merge(merge(merge(zz, data.frame(npx[[1]], ua20 = T), by = c(1,2), all = T),
            data.frame(npx[[2]], ua40 = T), by = c(1,2), all = T),
            data.frame(npx[[3]], ua60 = T), by = c(1,2), all = T),
            data.frame(npx[[4]], ua80 = T), by = c(1,2), all = T),
            data.frame(npx[[5]], ua100 = T), by = c(1,2), all = T)
hh[is.na(hh)] <- F

table(hh[,3:4]); table(hh[,4:5]); table(hh[,5:6]); table(hh[,6:7]); table(hh[,c(7, 3)])

length(which(apply(hh[3:7], 1, all)))

count(hh[3:7])[rev(order(count(hh[3:7])[,6])),]

# get bad pixel map & compare: which of those px are already identified?
{
    # temporary functions to decouple from pw.m
    {
        spots <- function(im, smooth.span = 1/5, min.diam = 5, edge.width = 10, auto.threshold = T) {
            
            sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
            
            # apply lowess smoothing
            smoo <- lowess.per.column(im, span = smooth.span)
            res <- im - smoo
            
            # flatten further by setting brighter pixels to mean value 
            res[res > mean(res)] <- mean(res)
            
            # dilate resulting image
            dilated <- dilate(res, sk)
            
            # erode resulting image (complete morphological closing)
            eroded <- erode(dilated, sk)
            
            # use k-means thresholding to identify spots
            # use 1-thresholded value to assign 1 to spots, 0 to background
            if (auto.threshold) {
                dim <- array(1, dim = c(1996, 1996)) - threshold(eroded, method = "kmeans")
            } else {
                dim <- array(1, dim = c(1996, 1996)) - threshold(eroded, mean(eroded) - 3*sd(eroded))
            }
            
            #------------------------------------------------------------------------------
            # convert to raster & clump values to identify individual spots
            # (values need to be reordered to maintain image orientation)
            blobs <- clump(raster(t(dim[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5), dir = 4)
            
            # summarise & remove any clusters that are too close to an edge
            # check loCation of cluster (looking for edge clusters)
            sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                                   id = getValues(blobs)[!is.na(getValues(blobs))]),
                        .(id), summarise, 
                        xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y),
                        xm = round(mean(x),0), ym = round(mean(y),0))
            
            sc$n.id <- sc$id
            sc$n.id[(sc$ym >= 1996 - edge.width) | (sc$ym <= edge.width) | 
                        (sc$xm >= 1996 - edge.width) | (sc$xm <= edge.width)] <- NA
            
            zz <- subs(blobs, sc[,c("id", "n.id")])
            
            if (all(is.na(getValues(zz)))) {
                # return empty data frame in order to not 
                return(data.frame("row" = double(), "col" = double(), "type" = character()))
            } else {
                # extract coordinates of cells covered by spots identified, return as data frame
                return(setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                                c("row", "col", "type")))
            }
        }  
    }
    
    md <- med.diffs(j.mean[[5]])
    mdb <- med.diffs(j.mean[[1]])
    bp <- rbind(data.frame(edge.px(j.mean[[1]]), type = ordered("edge", levels = c(Cat, "noisy"))),
                data.frame(which(j.mean[[1]] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(j.mean[[5]] == 0, arr.ind = T), type = "dead"),
                spots(j.mean[[5]]),
                get.dim.bright.px(j.mean[[1]]),
                get.dim.bright.px(j.mean[[5]]),
                data.frame(which(find.lines(j.mean[[5]]) > 0, arr.ind = T), type = "line.b"),
                data.frame(which(threshold(md, level = mad(j.mean[[5]]) * 2) > 0, arr.ind = T), type = "l.bright"),
                data.frame(which(threshold(md, level = mad(j.mean[[5]]) * -2) == 0, arr.ind = T), type = "l.dim"),
                data.frame(which(threshold(mdb, level = mad(j.mean[[1]]) * 2) > 0, arr.ind = T), type = "l.bright"),
                data.frame(which(threshold(mdb, level = mad(j.mean[[1]]) * -2) == 0, arr.ind = T), type = "l.dim"))
    
    bp <- bp[order(bp$type),]
    bp <- bp[!dupliCated(bp[,1:2]),]
}

hh <- merge(hh, bp, by = c(1:2), all.x = T)
table(hh[,7:8], useNA = "ifany")
table(hh[,c(3, 8)], useNA = "ifany")

count(hh[3:7])[rev(order(count(hh[3:7])[,6])),]

hh$noise <- apply(hh[,3:7] * 1, 1, paste, collapse = ".")
hh$noise.n <- apply(hh[,3:7] * 1, 1, sum)

table(hh$noise.n)           
table(hh$noise, hh$type, useNA = "ifany")
table(hh$noise.n, hh$type, useNA = "ifany")

tt <- table(hh$noise.n, hh$type, useNA = "ifany")
tt <- cbind(tt, "TOTAL" = rowSums(tt))
tt <- tt[, colSums(tt) > 0]
tt[tt == 0] <- "-"
colnames(tt) <- gsub("[ ]", "", gsub("[.]", "", colnames(tt)))

write.csv(tt, paste0(fpath.fig, "tbl-noisy-by-type.csv"), quote = F)

# N  v.bright  bright  line.b  l.bright  s.bright  edge  v.dim  dim  l.dim  s.dim  <NA>   | TOTAL
# 1      51      41       2       204         5      6     0     0      0      0    87    |  396
# 2      46      28       2        85         1      4     0     0      0      0    23    |  189
# 3      26       9       1        34         1      2     0     0      0      0    10    |  83
# 4      23       7       0        25         0      0     0     0      0      0     4    |  59
# 5      12       9       0        21         1      0     0     0      0      0     7    |  50

# v. small number of pixels identified as noisy in all 5 images (50), only 7 of which aren't picked up already.

par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))
for (i in 1:5) {
    hist(j.sd[[i]][as.matrix(hh[,1:2])], breaks = c(0:1500) * 10, xlab = "", ylab = "", xlim = c(0,4000), ylim = c(0,70), 
         main = paste0("ua", c(20, 40, 60, 80, 100)[i], ": ", nrow(npx[[i]]), "px found"), col = "black")
    abline(v = th$th[i], col = "red")
    hist(j.sd[[i]][as.matrix(hh[is.na(hh$type),1:2])], breaks = c(0:1500) * 10, add = T, col = "cyan3", border = "cyan3")
}
par(mfrow = c(1,1))

for (i in 1:5) {
    png(paste0(fpath, "power-sd-plots-ua", c(20, 40, 60, 80, 100)[i] , ".png"), width = 600, height = 600, pointsize = 28); {
        par(mar = c(2, 2, 0.5, 1))
        
        hist(j.sd[[i]], breaks = c(0:7500) * 2, xlab = "", ylab = "", xlim = c(0,700), ylim = c(0, 70), col = "black",
             main = "")
        #main = paste0("ua", c(20, 40, 60, 80, 100)[i], ": ", nrow(npx[[i]]), "px, ", length(which(hh[is.na(hh$type),i + 2])), " unclassified"))
        hist(j.sd[[i]][as.matrix(hh[,1:2])], breaks = c(0:7500) * 2, add = T, col = "cyan3", border = "cyan3")
        hist(j.sd[[i]][as.matrix(hh[is.na(hh$type),1:2])], breaks = c(0:7500) * 2, add = T, col = "gold", border = "gold")
        
        abline(v = th$th[i], col = "red", lwd = 2)
    }
    dev.off()
}

#==============================================================================

# check timings of UA images

zz <- summarise.profiles(fpath = "./Other-data/Other-images/line_investigation/")

zz[order(zz$uA),c("uA", "kV", "Power", "ImageNotes")]



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

####################################################################################################

# NORMALISE BLACK IMAGE & RE-THRESHOLD                                                          ####

# black SDs (unlike white/grey) are not normally distributed.
hist.cols <- adjustcolor(c("gold", "skyblue", "cyan3", "green3", "orange", "purple",
                           "violet", "blue", "magenta3", "red", "grey", "blue"), alpha = 0.2)
hist(pw.sd[,,"black", 1], breaks = "fd", xlim = c(0,100), col = hist.cols[1], border = hist.cols[1])
for (i in 2:12) {
    hist(pw.sd[,,"black", i], breaks = "fd", add = T, col = hist.cols[i], border = hist.cols[i])
}

# only look at distribution of pixels not already known to be erratic/unreliable
sdn.b <- pw.sd[ , , "black", "160430"]
sdn.b[as.matrix(bp$"160430"[,1:2])] <- NA

sdn.g <- pw.sd[ , , "grey", "160430"]
sdn.g[as.matrix(bp$"160430"[,1:2])] <- NA


hist(pw.sd[ , , "black", "160430"], breaks = "fd", xlim = c(0,50), prob = T)
lines(c(0:500)/10, dJohnson(c(0:500)/10, JohnsonFit(sdn.b[!is.na(sdn.b)])), col = "orange", lwd = 2)

hist(pw.sd[ , , "grey", "160430"], breaks = "fd", xlim = c(0,350), prob = T)
lines(c(0:500), dnorm(c(0:500), mean(pw.sd[ , , "grey", "160430"]), sd(pw.sd[ , , "grey", "160430"])), col = "skyblue", lwd = 2)
lines(c(0:500), dJohnson(c(0:500), JohnsonFit(sdn.g[!is.na(sdn.g)])), col = "orange", lwd = 2)


# apply Johnson transformation
