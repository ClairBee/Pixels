
# misbehaving subpanel on 16-06-20

zz <- t(readTIFF("./Other-data/Other-images/misbehaving_panel.tif", as.is = T)[1996:1,,drop = F])

####################################################################################################

# THINGS TO CHECK WHEN FULL CALIBRATION SET ARRIVES

    #   per-panel gradients in all three images
    #   per-panel boxplots in all three images
    #   difference since last image at all power levels
    #   edge behaviour
    #   any damage/odd behaviour close to sensor?

####################################################################################################

# PLOTS ETC                                                                                     ####

pixel.image(zz)
draw.panels()

sp <- subpanels(zz)

pixel.image(sp[,,"L4"])

o.plot(zz[,256], xlim = c(255, 639), ylim = c(34000, 44000))
o.plot(zz[,512], xlim = c(255, 639), add = T, col = adjustcolor("blue", alpha = 0.4))
o.plot(zz[,768], xlim = c(255, 639), add = T, col = adjustcolor("green3", alpha = 0.4))
o.plot(zz[,992], xlim = c(255, 639), add = T, col = adjustcolor("orange", alpha = 0.4))
o.plot(zz[,1], xlim = c(255, 639), add = T, col = adjustcolor("purple", alpha = 0.4))

pdf("./Notes/Damaged-subpanel/fig/Panel-transect.pdf"); {
    o.plot(zz[256,], xlim = c(1, 992), col = adjustcolor("red", alpha = 0.4))
    o.plot(zz[512,], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(zz[384,], add = T)
    legend("topleft", pch = 20, col = c("red", "black", "orange"), legend = c("column 256", "column 384", "column 512"))
    dev.off()
}

abline(v = panel.edges()$x - 0.5, col = "red")

pdf("./Notes/Damaged-subpanel/fig/Panel-hist.pdf")
{
    hist(sp[,,"L4"], breaks = c(0:65535), xlim = c(30000, 42000), main = "Histogram of damaged subpanel")
    #hist(sp[,,"L3"], breaks = c(0:65535), xlim = c(30000, 42000), main = "Histogram of neighbouring subpanel")
    
    lines(c(38912, 40959, NA, 36864, 38911, NA, 36863, 34816), 
          c(19000, 19000, NA, 16000, 16000, NA, 4000, 4000), col = "red", lty = 2)
    
    lines(c(38911, 37888, NA, 36864, 37887, NA, 33792, 34815), 
          c(7000, 7000, NA, 5000, 5000, NA, 2000, 2000), col = "cyan3", lty = 2)
    
    legend("topleft", lty = 2, col = c("red", "cyan3"), legend = c("2048", "1024"), bty = "n")
    dev.off()
}

range(sp[,,"L4"], na.rm = T)

pdf("./Notes/Damaged-subpanel/fig/Panel-image.pdf"); {
    image(1:384, 1:992, zz[255:638, 1:992], breaks = 32768 + (-2:12) * 1024, asp = T,
          col = c("black", "midnightblue", "blue", "dodgerblue3", "cyan3", "greenyellow", "yellow", "gold", "orange",
                  "red3", "violetred", "magenta3", "purple", "orchid"))
    dev.off()
}

pixel.image()

range(zz[255:638, 1:992])
638-255
31597-32768
    49622

length(32768 + (-2:12) * 1024)
length(sd.colours())
length(sd.levels(zz))

#============================================================================

qq <- hist(sp[,,"L4"], breaks = c(0:65535), xlim = c(30000, 42000),
           main = "Histogram of damaged subpanel")

qq <- data.frame(mids = qq$mids, count = qq$counts)
qq <- qq[rev(order(qq$count)),]

qq.h <- qq[1:10,]
qq.h$mids[1:09] - qq.h$mids[2:10]

head(qq, 20)
qq <- qq[qq$count > 0,]
qq$diff <- qq$mids - min(qq$mids)

qq$l2 <- log2(qq$diff)

qq$mod <- qq$diff %% 16

pixel.image(pw.m[,,"grey", "141009"])
pixel.image(pw.m[,,"grey", "141118"])

load.pixel.means()

mean(zz); mean(pw.m[,,"grey", "160430"]); mean(pw.m[,,"white", "160430"])

# what has changed since last images were taken?
ch <- zz - pw.m[,,"grey", "160430"]

pixel.image(ch)

####################################################################################################

# SUBPANELS                                                                                     ####

sp <- subpanels(zz)

# per-subpanel boxplots
dimnames(sp) <- list(NULL, NULL, apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                             rep(c(rep("0", 9), rep("", 7)), 2),
                                             rep(c(1:16), 2)), 1, paste, collapse = ""))

boxplot(value ~ X3, melt(sp[,,1:16]), pch = 20, ylim = c(25000, 50000))
boxplot(value ~ X3, melt(sp[,,17:32]), pch = 20, ylim = c(25000, 50000))

# upper half of detector less variable with more outliers - because spot is offset?

# per-panel gradients?
zz.panels <- panel.lm(zz)
prev.panels.grey <- panel.lm(pw.m[,,"grey", "160430"])
prev.panels.white <- panel.lm(pw.m[,,"white", "160430"])

pixel.image(zz.panels$fitted.values)

####################################################################################################

# COLUMN 427                                                                                    ####

o.plot(zz[427,], xlim = c(993, 1996))
lines(zz[426,], col = adjustcolor("cyan3", alpha = 0.4))
lines(zz[428,], col = adjustcolor("green3", alpha = 0.4))



o.plot(zz[427,] - zz[426,], xlim = c(993, 1996))
lines(pw.m[427,,"grey", "160430"] - pw.m[426,,"grey", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
lines(pw.m[427,,"black", "160314"] - pw.m[426,,"black", "160314"], col = adjustcolor("orange", alpha = 0.4))


focal.plot(pw.m[,,"grey", "160430"], centre = c(427, 1200), bad.px = bp$"160430", bpx.cex = 5)
focal.plot(zz, centre = c(427, 1200), bad.px = bp$"160430", bpx.cex = 5)

# bright lines in the new image?
ll <- find.lines(zz)

xy <- which(ll > 0, arr.ind = T)
xy <- data.frame(xy, ll = ll[xy])

lll <- ddply(xy, .(ll, x = row), summarise,
             length = length(row), ymin = min(col), ymax = max(col))

lll <- lll[rev(order(lll$length)),]
