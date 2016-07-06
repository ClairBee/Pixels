
# ALTERNATIVE THRESHOLDING APPROACHES

library("IO.Pixels"); library("CB.Misc")

fpath <- "./Other-data/"

acq <- readRDS("./02_Objects/images/pwm-160430.rds")
md <- readRDS("./02_Objects/med-diffs/md-160430.rds")

####################################################################################################

# higher threshold on local px gives greater stability, but only up to a point
# include plot showing movement of screen spots & why they must be tracked each time

####################################################################################################

# COMPARE EXTREME-VALUED PIXELS IN GREY & BLACK                                                 ####

th <- apply(acq[,,c("black", "grey")], 3, function(im) {
    med <- median(im, na.rm = T)
    c(v.dim = med * 0.5, dim = med * 0.75,
      bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
})

{
    chk.vb <- length(which(acq[,,"black"] > median(acq[,,"black"], na.rm = T) + (65535 - median(acq[,,"black"], na.rm = T)) / 2 & 
                               acq[,,"grey"] > median(acq[,,"grey"], na.rm = T) + (65535 - median(acq[,,"grey"], na.rm = T)) / 2)) / 
    length(which(acq[,,"black"] > median(acq[,,"black"], na.rm = T) + (65535 - median(acq[,,"black"], na.rm = T)) / 2 | 
                     acq[,,"grey"] > median(acq[,,"grey"], na.rm = T) + (65535 - median(acq[,,"grey"], na.rm = T)) / 2))
    
    chk.b <- length(which(acq[,,"black"] > median(acq[,,"black"], na.rm = T) + (65535 - median(acq[,,"black"], na.rm = T)) / 4 & 
                              acq[,,"grey"] > median(acq[,,"grey"], na.rm = T) + (65535 - median(acq[,,"grey"], na.rm = T)) / 4)) / 
        length(which(acq[,,"black"] > median(acq[,,"black"], na.rm = T) + (65535 - median(acq[,,"black"], na.rm = T)) / 4 | 
                         acq[,,"grey"] > median(acq[,,"grey"], na.rm = T) + (65535 - median(acq[,,"grey"], na.rm = T)) / 4))
}

# countour plot of grey threshold vs black threshold vs % same
# black threshold from 5134 to 65535; grey from 17506
{
    zz.arr <- array(dim = c(13, 11), dimnames = list(c(1:13) * 5000, c(3:13) * 5000))
    zz <- melt(zz.arr, varnames = c("b", "g"))[,1:2]
    zz$p <- apply(zz, 1, 
                  function(x) length(which(acq[,,"black"] > x["b"] & acq[,,"grey"] > x["g"])) / 
                      length(which(acq[,,"black"] > x["b"] | acq[,,"grey"] > x["g"])))
    
    zz.arr <- array(zz$p, dim = dim(zz.arr), dimnames = dimnames(zz.arr))
    
    contour(c(1:13) * 5000, c(3:13) * 5000, zz.arr, xlab = "black", ylab = "grey")
    image(c(1:13) * 5000, c(3:13) * 5000, zz.arr, xlab = "black", ylab = "grey")
    
    abline(v = median(acq[,,"black"], na.rm = T), col = "red")
    abline(v = (3 * median(acq[,,"black"], na.rm = T) + 65535) / 4, col = "red", lty = 2)
    abline(v = (median(acq[,,"black"], na.rm = T) + 65535) / 2, col = "red", lty = 3)
    
    abline(h = median(acq[,,"grey"], na.rm = T), col = "red")
    abline(h = (3 * median(acq[,,"grey"], na.rm = T) + 65535) / 4, col = "red", lty = 2)
    abline(h = (median(acq[,,"grey"], na.rm = T) + 65535) / 2, col = "red", lty = 3)
    
    rect((3 * median(acq[,,"black"], na.rm = T) + 65535) / 4, (3 * median(acq[,,"grey"], na.rm = T) + 65535) / 4,
         65535, 65535, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect((median(acq[,,"black"], na.rm = T) + 65535) / 2, (median(acq[,,"grey"], na.rm = T) + 65535) / 2,
         65535, 65535, col = adjustcolor("red", alpha = 0.3), border = NA)
}

# check degree of correlation between high-valued px in grey vs black images
{
    plot(acq[,,"black"][which(acq[,,"black"] > 20000 | acq[,,"grey"] > 20000, arr.ind = T)],
         acq[,,"grey"][which(acq[,,"black"] > 20000 | acq[,,"grey"] > 20000, arr.ind = T)], pch = 20,
         xlab = "Value in black images", ylab = "Value in grey images")
}

smoothScatter(acq[,,"black"], acq[,,"grey"], nrpoints = Inf,
              xlab = "Value in black images", ylab = "Value in grey images")
abline(line(acq[,,"black"][which((acq[,,"black"] > 10000 | acq[,,"grey"] > 20000) & acq[,,"grey"] < 65535, arr.ind = T)],
            acq[,,"grey"][which((acq[,,"black"] > 10000 | acq[,,"grey"] > 20000) & acq[,,"grey"] < 65535, arr.ind = T)]),
       col = adjustcolor("darkred", alpha = 0.4))

smoothScatter(acq[,,"black"], acq[,,"white"], nrpoints = Inf,
              xlab = "Value in black images", ylab = "Value in white images")
abline(line(acq[,,"black"][which((acq[,,"black"] > 10000 | acq[,,"white"] > 51000) & acq[,,"white"] < 65535, arr.ind = T)],
            acq[,,"white"][which((acq[,,"black"] > 10000 | acq[,,"white"] > 51000) & acq[,,"white"] < 65535, arr.ind = T)]),
       col = adjustcolor("darkred", alpha = 0.4))

smoothScatter(acq[,,"grey"], acq[,,"white"], nrpoints = Inf,
              xlab = "Value in grey images", ylab = "Value in white images")
abline(line(acq[,,"grey"][which((acq[,,"grey"] > 23000 | acq[,,"white"] > 51000) & acq[,,"white"] < 65535, arr.ind = T)],
            acq[,,"white"][which((acq[,,"grey"] > 23000 | acq[,,"white"] > 51000) & acq[,,"white"] < 65535, arr.ind = T)]),
       col = adjustcolor("darkred", alpha = 0.4))

abline(v = 20000)

####################################################################################################

# BASIC THRESHOLDING - BLACK & GREY                                                             ####

org.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "skyblue", "green3")

# get thresholds for extreme-valued pixels
th <- apply(acq[,,c("black", "grey")], 3, function(im) {
    med <- median(im, na.rm = T)
    c(v.dim = med * 0.5, dim = med * 0.75,
      bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
})

px.org <- rbind(data.frame(edge.px(acq, edge.width = 40), type = "edge"),
            data.frame(no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]), type = "no.resp"),
            data.frame(which(acq[,,"black"] == 65535, arr.ind = T), type = "hot"),
            data.frame(which(acq[,,"grey"] == 0, arr.ind = T), type = "dead"),
            data.frame(screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40), type = "screen.spot"),
            data.frame(which(find.lines(acq[, , "black"], midline = 1024.5) > 0, arr.ind = T), type = "line.b"),
            data.frame(which(acq[, , "black"] > th["v.bright", "black"], arr.ind = T), type = "v.bright"),
            data.frame(which(acq[, , "black"] > th["bright", "black"], arr.ind = T), type = "bright"),
            data.frame(which(acq[, , "black"] < th["v.dim", "black"], arr.ind = T), type = "v.dim"),
            data.frame(which(acq[, , "black"] < th["dim", "black"], arr.ind = T), type = "dim"),
            data.frame(which(acq[, , "grey"] > th["v.bright", "grey"], arr.ind = T), type = "v.bright"),
            data.frame(which(acq[, , "grey"] > th["bright", "grey"], arr.ind = T), type = "bright"),
            data.frame(which(acq[, , "grey"] < th["v.dim", "grey"], arr.ind = T), type = "v.dim"),
            data.frame(which(acq[, , "grey"] < th["dim", "grey"], arr.ind = T), type = "dim"),
            data.frame(which(md[, , "black"] > 1200, arr.ind = T), type = "l.bright"),
            data.frame(which(md[, , "grey"] > 1200, arr.ind = T), type = "l.bright"),
            data.frame(which(md[, , "black"] < -1200, arr.ind = T), type = "l.dim"),
            data.frame(which(md[, , "grey"] < -1200, arr.ind = T), type = "l.dim"))

Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "edge", "l.bright", "line.d", "screen.spot", "v.dim", "dim", "l.dim")
px.org$type <- ordered(px.org$type, levels = Cat)

px.org <- px.org[order(px.org$type),]
px.org <- px.org[!duplicated(px.org[,1:2]),]
px.org <- px.org[!px.org$type %in% c("screen.spot"),]

sc.pad <- array(dim = dim(acq[,,1]))
sc.pad[3:1998, 33:2028] <- sc[,,12]

smoothScatter(acq[,,"grey"], sc.pad, nrpoints = Inf, xlim = c(0,65535),
              ylab = "Shading-corrected", xlab = "Pixelwise mean in grey images")

points(acq[,,"grey"][as.matrix(px[,1:2])], sc.pad[as.matrix(px[,1:2])],
       pch = ".", cex = 2,  col = Cat.cols[px$type])

# now remove identified bad pixels & check shading correction again
sc.corr <- sc.pad
sc.corr[as.matrix(px[,1:2])] <- NA

smoothScatter(acq[,,"grey"], sc.corr, nrpoints = Inf, 
              ylab = "Shading-corrected", xlab = "Pixelwise mean in grey images")

zz <- which(sc.corr < 12000, arr.ind = T)

focal.plot(acq[,,"black"], c(116, 48) + c(2, 32))       # locally bright in black (+ 8900)
focal.plot(acq[,,"black"], c(1496, 193) + c(2, 32))     # locally bright in black (+ 6800)
focal.plot(acq[,,"black"], c(572, 604) + c(2, 32))      # locally bright in black (+ 7600)
focal.plot(acq[,,"black"], c(163, 1017) + c(2, 32))     # locally bright in black (+ 7900)

focal.plot(acq[,,"black"], c(119, 1511) + c(2, 32))     # locally bright in black (+ 3700)
focal.plot(acq[,,"black"], c(90, 1532) + c(2, 32))      # locally bright in black (+ 9942)

md[90+2, 1532+32,]

####################################################################################################

# QUALITATIVE THRESHOLDING                                                                      ####

th <- apply(acq, 3, function(im) {
    med <- median(im, na.rm = T)
    c(v.dim = med * 0.5, dim = med * 0.75,
      bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
})

# identify points that are too bright in black image
bright <- which(acq[,,"black"] > th["bright", "black"], arr.ind = T)
plot(bright, pch = 20)

# then points that are too bright in grey image
warm <- which(acq[,,"grey"] > th["bright", "grey"] & acq[,,"black"] <= th["bright", "black"], arr.ind = T)
points(warm, col = "red")

# then points with non-linear response (high residual) in white image - excluding edges
{
    # data frame of all variables for active region of image
    df <- setNames(data.frame(melt(acq[,,"black"]), 
                              melt(acq[,,"grey"]),
                              melt(acq[,,"white"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    df <- df[findInterval(df$x, c(40.5, 2008.5)) == 1 & findInterval(df$y, c(40.5, 2008.5)) == 1,]
    df <- df[!is.na(df$b),]
    
    
    # fit linear model
    w.lm <- rlm(w ~ b * g, data = df)
    
    df$fv <- w.lm$fitted.values
    df$res <- w.lm$residuals
    
    smoothScatter(df$fv, df$res, ylim = c(-3000, 3000), colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
}
nonlinear <- setNames(df[abs(df$res) > 1000, 1:2], nm = c("row", "col"))
points(nonlinear, col = "cyan3")

# then points that are non-responsive (same in black & grey)
stuck <- which(acq[,,"grey"] - acq[,,"black"] < 500 & acq[,,"grey"] > acq[,,"black"], arr.ind = T)
{
    diffs <- acq[,,"grey"] - acq[,,"black"]
    hist(diffs, breaks = "fd", ylim = c(0,30))
    abline(v = 500, col = "red", lty = 2)
    summary(c(diffs))
    

    smoothScatter(acq[,,"black"][px], acq[,,"grey"][px], nrpoints = Inf)
    smoothScatter(acq[,,"black"][px], acq[,,"white"][px], nrpoints = Inf)
    
    abline(0,1, col = adjustcolor("darkred", alpha = 0.4), lty = 3)
    
    sum(acq[,,"grey"][px] > 50000 & acq[,,"black"][px] > 60000)
}
points(stuck, pch = 15, col = "purple")

# other points that are dim in black images
dim <- which(acq[,,"black"] < th["dim", "black"], arr.ind = T)
cool <- which(acq[,,"grey"] < th["dim", "grey"], arr.ind = T)

# bright / dim lines as usual - don't bother reworking

# gather all bad pixels
px.qual <- rbind(data.frame(bright, type = "bright"),
                 data.frame(warm, type = "warm"),
                 data.frame(nonlinear, type = "nonlinear"),
                 #data.frame(cool, type = "cool"),
                 data.frame(stuck, type = "stuck"),
                 data.frame(dim, type = "dim"))


px.qual$type <- ordered(px.qual$type, levels = c("bright", "warm", "stuck", "nonlinear", "dim", "cool"))
px.qual <- px.qual[order(px.qual$type),]
px.qual <- px.qual[!duplicated(px.qual[,1:2]),]

qual.cols <- c("magenta3", "orange",  "blue", "cyan3", "forestgreen", "green")
plot(px.qual[,1:2], pch = 20, col = qual.cols[px.qual$type])

saveRDS(px.qual, paste0(fpath, "bpm-qualitative.rds"))

bp <- merge(px.org, px.qual, by = c(1,2), all = T, suffix = c(".org", ".qual"))

table(bp$type.org, bp$type.qual, useNA = "ifany")

sc.qual <- sc.org <- sc <- 60000 * (acq[,,"grey"] - acq[,,"black"]) / (acq[,,"white"] - acq[,,"black"])
sc.qual[is.infinite(sc.qual)] <- 0; sc.org[is.infinite(sc.org)] <- 0

# ignore defective pixels
sc.qual[as.matrix(px.qual[,1:2])] <- NA
sc.org[as.matrix(px.org[,1:2])] <- NA

# ignore edges
smoothScatter(acq[,,"black"][41:2008,41:2008], sc.qual[41:2008,41:2008], asp = T, nrpoints = 0, xlim = c(4000,9000),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
points(acq[,,"black"][as.matrix(px.qual[,1:2])], sc[as.matrix(px.qual[,1:2])],
       pch = ".", cex = 2, col = qual.cols[px.qual$type])

smoothScatter(acq[,,"black"][41:2008,41:2008], sc.org[41:2008,41:2008], asp = T, nrpoints = 0, xlim = c(4000,9000),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
points(acq[,,"black"][as.matrix(px.org[,1:2])], sc[as.matrix(px.org[,1:2])],
       pch = ".", cex = 2, col = org.cols[px.org$type])
points(acq[,,"black"][as.matrix(px.qual[,1:2])],
       sc[as.matrix(px.qual[,1:2])],
       cex = 0.3, col = adjustcolor("green3", alpha = 0.3))

####################################################################################################

# DARK-ADJUSTED ERRORS                                                                          ####

# identify errors in black images as usual

# then subtract black image from grey and look for errors in the resulting image

####################################################################################################

# SCREEN SPOTS - GREY VS WHITE                                                                  ####

# white image gives clearer view in this case - dips stand out against noise
{
    ss.g <- screen.spots(acq[,,"grey"], enlarge = T, ignore.edges = 40)
    ss.w <- screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40)
    
    plot(ss.g, pch = 15, asp = T, cex = 0.1, col = adjustcolor("cyan3", alpha = 0.2),
         xlim = c(0,2048), ylim = c(0,2048), xlab = "", ylab = "")
    points(ss.w, pch = 15, asp = T, cex = 0.1, col = adjustcolor("gold", alpha = 0.2))
    draw.panels.2048(col = "grey")
    
    plot(ss.w, xlim = c(960, 1030), ylim = c(540, 600), pch = 15, asp = T, cex = 0.6)
    pixel.image(acq[,,"white"], xlim = c(960, 1030), ylim = c(540, 600))
    points(ss.w, pch = 1)
    
    o.plot(acq[,570,"grey"], xlim = c(900, 1100))
    points(acq[,570,"grey"], col = c(NA, "red")[c(1:1200) %in% ss.w[ss.w[,2] == 570,1] + 1])
    
    o.plot(acq[,570,"white"] - mean(acq[900:1100,570,"white"]), xlim = c(900, 1100))
    points(acq[,570,"white"] - mean(acq[900:1100,570,"white"]),
           col = c(NA, "red")[c(1:1200) %in% ss.w[ss.w[,2] == 570,1] + 1])

}

####################################################################################################

# EXAMINE BRIGHT LINES                                                                          ####

summarise.lines(find.lines(acq[, , "black"], midline = 1024.5, threshold.at = 4000), midline = 1024.5)
# nothing gained by reducing threshold - only noise

# what about pattern of line segments?
{
    ll.summ <- summarise.lines(find.lines(acq[, , "black"], midline = 1024.5), midline = 1024.5)
    unlist(lapply(unique(ll.summ$col),
                  function(cc) max(ll.summ$ymin[ll.summ$col == cc]) - 
                      min(ll.summ$ymax[ll.summ$col == cc])))
    
    # both segment breaks are approximately 200px in length
    
    # plots aligned so that central panel is at LHS
    o.plot(acq[429,,"black"] - acq[428,,"black"], xlim = c(1025, 2048), ylim = c(-1000,2000))
    o.plot(acq[811,,"black"] - acq[810,,"black"], xlim = c(1024, 0), ylim = c(-1000,2000))
    
    plot(acq[429,1025:2048,"black"] - acq[428,1025:2048,"black"], ylim = c(-1000,2000), type = "l", col = "blue",
         xlab = "Distance from midline ->", ylab = "Difference from neighbour grey value")
    lines(acq[811,1024:1,"black"] - acq[810,1024:1,"black"], col = "red")
}
# not only is length between line segments very similar - also very similar brightness in 
# 'tail' segment closest to midline


# brightness vs neighbouring columns at each power setting?
{
    o.plot(acq[429,,"black"] - acq[428,,"black"], xlim = c(1025, 2048), ylim = c(-1000,2000))
    lines(acq[429,,"grey"] - acq[428,,"grey"], col = adjustcolor("cyan3", alpha = 0.5))
    lines(acq[429,,"white"] - acq[428,,"white"], col = adjustcolor("gold", alpha = 0.5))
    
    o.plot(acq[811,,"black"] - acq[812,,"black"], xlim = c(0,1024), ylim = c(-1000,2000))
    lines(acq[811,,"grey"] - acq[812,,"grey"], col = adjustcolor("cyan3", alpha = 0.5))
    lines(acq[811,,"white"] - acq[812,,"white"], col = adjustcolor("gold", alpha = 0.5))
    
    # normal column, for comparison
    o.plot(acq[429,,"black"] - acq[428,,"black"], xlim = c(0,1024), ylim = c(-1000,2000))
    lines(acq[429,,"grey"] - acq[428,,"grey"], col = adjustcolor("cyan3", alpha = 0.5))
    lines(acq[429,,"white"] - acq[428,,"white"], col = adjustcolor("gold", alpha = 0.5))
}
# both columns show extremely similar behaviour at each power setting.
# suggests additional electrical charge 'leaking' along column?

# check for any other types of line
{
    # no dim lines found
    tt <- find.lines(acq[, , "black"], midline = 1024.5, dim.lines = T, horizontal = F)
    pixel.image(tt)
    summarise.lines(tt, midline = 1024.5)
    
    # somehow finding verical lines - problem with smoothing perhaps? Check.
    tt <- find.lines(acq[, , "black"], midline = 1024.5, dim.lines = T, horizontal = T)
    pixel.image(tt)
    summarise.lines(tt, midline = 1024.5)
    
    # somehow finding verical lines - problem with smoothing perhaps? Check.
    tt <- find.lines(acq[, , "black"], midline = 1024.5, dim.lines = F, horizontal = T)
    pixel.image(tt)
    summarise.lines(tt, midline = 1024.5)
}

# progression over time
{
    # compare also to single image from defect panel investigation
    zz <- t(readTIFF("./Other-data/Other-images/misbehaving_panel.tif", as.is = T)[1996:1,,drop = F])
    acq.160314 <- readRDS("./02_Objects/images/pwm-160314.rds")
    
    o.plot(acq.160314[429,,"black"] - acq.160314[428,,"black"], 
           xlim = c(1025, 2048), ylim = c(-1000,2000))
    lines(acq[429,,"black"] - acq[428,,"black"], col = "blue")
    lines(zz[427,] - zz[426,], col = adjustcolor("purple", alpha = 0.5))
    
    focal.plot(acq[,,"black"], c(x = 429, y = 1200))
}

####################################################################################################

# LOCALLY BRIGHT/DIM PX                                                                         ####

# better to use SD or use MAD?
{
    s.hist(md[,,"black"])
    abline(v = 2 * mad(acq[,,"black"], na.rm = T), col = "red")
    
    s.hist(md[,,"grey"])
    abline(v = 2 * mad(acq[,,"grey"], na.rm = T), col = "red")
    
    # significantly fewer points when using SD rather than MAD.
    lb.b.sd <- which(md[,,"black"] > 2 * sd(acq[,,"black"], na.rm = T), arr.ind = T)
    lb.b.mad <- which(md[,,"black"] > 2 * mad(acq[,,"black"], na.rm = T), arr.ind = T)
    lb.b.sd3 <- which(md[,,"black"] > 3 * sd(acq[,,"black"], na.rm = T), arr.ind = T)
    lb.b.sd4 <- which(md[,,"black"] > 4 * sd(acq[,,"black"], na.rm = T), arr.ind = T)
    
    
    lb.g.sd <- which(md[,,"grey"] > 2 * sd(acq[,,"grey"], na.rm = T), arr.ind = T)
    lb.g.mad <- which(md[,,"grey"] > 2 * mad(acq[,,"grey"], na.rm = T), arr.ind = T)
    lb.g.sd3 <- which(md[,,"grey"] > 3 * sd(acq[,,"grey"], na.rm = T), arr.ind = T)
    lb.g.sd4 <- which(md[,,"grey"] > 4 * sd(acq[,,"grey"], na.rm = T), arr.ind = T)
    
    
    nrow(lb.b.sd); nrow(lb.b.mad); nrow(lb.g.sd); nrow(lb.g.mad)
    # much greater similarity between grey & black lists when SD is used: 86% in common
    
    tt <- merge(data.frame(lb.b.sd, b = 1), data.frame(lb.g.sd, g = 1), by = c(1, 2), all = T)
    tt[is.na(tt)] <- 0
    tt$ttl <- tt$b + tt$g
    round(table(grey = tt$g, black = tt$b) / nrow(tt) * 100, 2)

    mm <- merge(data.frame(lb.b.mad, b = 1), data.frame(lb.g.mad, g = 1), by = c(1, 2), all = T)
    mm[is.na(mm)] <- 0
    mm$mml <- mm$b + mm$g
    round(table(grey = mm$g, black = mm$b) / nrow(mm) * 100, 2)
    
    tt3 <- merge(data.frame(lb.b.sd3, b = 1), data.frame(lb.g.sd3, g = 1), by = c(1, 2), all = T)
    tt3[is.na(tt3)] <- 0
    tt3$tt3l <- tt3$b + tt3$g
    round(table(grey = tt3$g, black = tt3$b) / nrow(tt3) * 100, 2)
    
    tt4 <- merge(data.frame(lb.b.sd4, b = 1), data.frame(lb.g.sd4, g = 1), by = c(1, 2), all = T)
    tt4[is.na(tt4)] <- 0
    tt4$tt4l <- tt4$b + tt4$g
    round(table(grey = tt4$g, black = tt4$b) / nrow(tt4) * 100, 2)
}
# SD better. Plot # px identified as l. bright vs threshold (with SD multiples added afterwards)

hh <- lapply(c(3:15) * 100, 
             function(lim) length(which(md[,,"black"] > lim & md[,,"grey"] > lim)) / 
                 length(which(md[,,"black"] > lim | md[,,"grey"] > lim)))

plot(c(3:15) * 100, unlist(hh) * 100, pch = 20, ylab = "% points marked in both images", xlab = "threshold")
abline(v = c(mad(acq[,,"black"], na.rm = T), mad(acq[,,"grey"], na.rm = T)), col = "orange", lty = c(1,2))
abline(v = c(sd(acq[,,"black"], na.rm = T), sd(acq[,,"grey"], na.rm = T)), col = "green3", lty = c(1,2))
abline(v = 2 * c(sd(acq[,,"black"], na.rm = T), sd(acq[,,"grey"], na.rm = T)), col = "blue", lty = c(1,2))

points(c(sd(acq[,,"black"], na.rm = T), sd(acq[,,"grey"], na.rm = T)), rep(length(which(md[,,"black"] > sd(acq[,,"black"], na.rm = T) & md[,,"grey"] >  sd(acq[,,"grey"], na.rm = T))) / 
           length(which(md[,,"black"] >  sd(acq[,,"black"], na.rm = T) | md[,,"grey"] > sd(acq[,,"grey"], na.rm = T))) * 100, 2), col = "green3", type = "o", lty = 4)

points(c(mad(acq[,,"black"], na.rm = T), mad(acq[,,"grey"], na.rm = T)), rep(length(which(md[,,"black"] > mad(acq[,,"black"], na.rm = T) & md[,,"grey"] >  mad(acq[,,"grey"], na.rm = T))) / 
           length(which(md[,,"black"] >  mad(acq[,,"black"], na.rm = T) | md[,,"grey"] > mad(acq[,,"grey"], na.rm = T))) * 100, 2), col = "orange", type = "o", lty = 4)

legend("bottomright", col = c("orange", "green3", "blue"), lty = 2, 
       legend = c("1 x MAD", "1 x SD", "2 x SD"))

# spatial plots of points identified
{
    plot(which(md[,,"black"] > 2 * mad(acq[,,"black"], na.rm = T), arr.ind = T), pch = 16, cex = 0.2,
         col = "cyan3")
    points(which(md[,,"black"] > 2 * sd(acq[,,"black"], na.rm = T), arr.ind = T), pch = 16, cex = 0.2,
           col = "red")
}

# plots of threshold vs % persistent points in all acquisitions
{
    load.pixel.means()
    md.b <- readRDS("./Other-data/Median-diffs-black.rds")
    md.g <- readRDS("./Other-data/Median-diffs-grey.rds")
    
    hh <- do.call("rbind", lapply(names(md.b),
                 function(dt) unlist(lapply(c(3:15) * 100, 
                                     function(lim) length(which(md.b[[dt]] > lim & md.g[[dt]] > lim)) / 
                                         length(which(md.b[[dt]] > lim | md.g[[dt]] > lim))))))
    
    # proportion of bright pixels matched in grey/black images for each acquisition
    pdf("./Image-plots/Misc/Local-threshold.pdf"); {
        par(mar = c(3, 3, 1, 1), mfrow = c(4, 3))
        for (i in (1:dim(pw.m)[[4]])) {
            plot(c(3:15) * 100, hh[i,] * 100, pch = 20, xlab = "Threshold", ylab = "% bright in both images", 
                 ylim = c(0,100), cex.axis = 0.6, cex = 0.8)
            abline(v = sd(pw.m[,,"black", i], na.rm = T) * c(1, 2), col = "blue", lty = c(1,2))
            abline(v = sd(pw.m[,,"grey", i], na.rm = T) * c(1, 2), col = "green3", lty = c(1,2))
            abline(h = c(0:10) * 10, col = adjustcolor("grey", alpha = 0.5), lty = 3)
            
            legend("bottomright", col = c("blue", "green3", "blue", "green3"), lty = c(1, 1, 2, 2), 
                   title = fancy.date(dimnames(pw.m)[[4]][i]), bg = "white", cex = 0.7,
                   legend = c("Black SD", "Grey SD", "2 x black SD", "2 x grey SD"))
        }
        dev.off()
    }
}

####################################################################################################

# PREDICT WHITE VALUE FROM GREY/BLACK IMAGES                                                    ####
fpath <- "./Notes/WV-prediction/fig/"

{
    # data frame of all variables for active region of image
    df <- merge(setNames(data.frame(melt(acq[,,"black"]), 
                              melt(acq[,,"grey"]),
                              melt(acq[,,"white"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w")),
                px, by = c(1:2), all.x = T)
    df <- df[!is.na(df$b),]
    
    # fit linear model to healthy px only, check line through fitted/actual
    rlm <- rlm(w ~ b * g, data = df[is.na(df$type),])       # 231.160       0.995
    hlm <- lm(w ~ b * g, data = df[is.na(df$type),])        # 167.6664      0.9964
    alm <- lm(w ~ b * g, data = df)                         # 236.2960      0.9949
    
    hlm.p <- predict(hlm, interval = "conf")
    
    df$h.fv[is.na(df$type)] <- hlm$fitted.values
    df$h.res[is.na(df$type)] <- hlm$residuals
    
    df$a.fv <- alm$fitted.values
    df$a.res <- alm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(alm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(alm)$sigma, 2)),
          paste0(fpath, "fitted-wv-all.txt"))
    write(paste0("Adj. $r^2$ ", round(summary(hlm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(hlm)$sigma, 2)),
          paste0(fpath, "fitted-wv-healthy.txt"))
    write(paste0(sum(abs(hlm$residuals) > 2 * sd(df$h.fv, na.rm = T)), " px > ", round(2 * sd(df$h.fv, na.rm = T), 0), " res"),
          paste0(fpath, "fitted-wv-res.txt"))
    
    pdf(paste0(fpath, "fitted-wv-all-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$w, df$a.fv, xlim = c(0,65535), ylim = c(0,65535),
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Observed white value", ylab = "Fitted white value")
        abline(line(df$w, df$a.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
        dev.off()
    }
   
    pdf(paste0(fpath, "fitted-wv-healthy-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$w, df$h.fv, xlim = c(0,65535), ylim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed white value", ylab = "Fitted white value")
        abline(line(df$w, df$h.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
        dev.off()
    }
    
    pdf(paste0(fpath, "fitted-wv-healthy-res.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$h.fv, df$h.res, xlim = c(0,65535), ylim = c(-4000,4000),
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Fitted white value", ylab = "Residual")
        abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
        abline(h = 1200 * c(-1,1), col = adjustcolor("darkblue", alpha = 0.4), lty = 2)  
        dev.off()
    }
}

####################################################################################################

# PREDICT WHITE VALUE IN OLD DATA - 131122                                                      ####

pw.m <- readRDS("./Other-data/Old-data/Pixelwise-means.rds")
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")$"131122"
bp <- bp[!bp$type %in% c("s.bright", "l.bright", "line.b", "l.dim", "s.dim"),]

# create data frame
{
    dfo <- setNames(data.frame(melt(pw.m[,,"black", "131122"]), 
                               melt(pw.m[,,"grey", "131122"]),
                               melt(pw.m[,,"white", "131122"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                    nm = c("x", "y", "b", "g", "w"))
    
    dfo <- merge(dfo, bp, by = c(1:2), all.x = T)
    dfo <- dfo[!is.na(dfo$b),]
}

# fit models
{
    hlm <- lm(w ~ b * g, data = dfo[is.na(dfo$type),])        # 1814.3848     0.9588
    alm <- lm(w ~ b * g, data = dfo)                         # 236.2960      0.9949
    
    dfo$h.fv[is.na(dfo$type)] <- hlm$fitted.values
    dfo$h.res[is.na(dfo$type)] <- hlm$residuals
    
    dfo$a.fv <- alm$fitted.values
    dfo$a.res <- alm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(alm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(alm)$sigma, 2)),
          paste0(fpath, "fitted-wv-all-old.txt"))
    write(paste0("Adj. $r^2$ ", round(summary(hlm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(hlm)$sigma, 2)),
          paste0(fpath, "fitted-wv-healthy-old.txt"))
    write(paste0(sum(abs(hlm$residuals) > 2 * sd(df$h.fv, na.rm = T)), " px > ", round(2 * sd(df$h.fv, na.rm = T), 0), " res"),
          paste0(fpath, "fitted-wv-res-old.txt"))
}

pdf(paste0(fpath, "fitted-wv-all-px-old.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$w, dfo$a.fv, xlim = c(0,65535), ylim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed white value", ylab = "Fitted white value")
    abline(line(dfo$w, dfo$a.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
}

pdf(paste0(fpath, "fitted-wv-healthy-px-old.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$w, dfo$h.fv, xlim = c(0,65535), ylim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed white value", ylab = "Fitted white value")
    abline(line(dfo$w, dfo$h.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    dev.off()
}

pdf(paste0(fpath, "fitted-wv-healthy-res-old.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$h.fv, dfo$h.res, xlim = c(0,65535), ylim = c(-4000,4000),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Fitted white value", ylab = "Residual")
    abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    abline(h = 1200 * c(-1,1), lty = 2, col = adjustcolor("darkblue", alpha = 0.2))
    dev.off()
}

####################################################################################################

# PREDICT WHITE VALUE IN OLD DATA - 140128                                                      ####

pw.m <- readRDS("./Other-data/Old-data/Pixelwise-means.rds")
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")$"140128"
bp <- bp[!bp$type %in% c("s.bright", "l.bright", "line.b", "l.dim", "s.dim"),]

# create data frame
{
    dfo <- setNames(data.frame(melt(pw.m[,,"black", "140128"]), 
                               melt(pw.m[,,"grey", "140128"]),
                               melt(pw.m[,,"white", "140128"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                    nm = c("x", "y", "b", "g", "w"))
    
    dfo <- merge(dfo, bp, by = c(1:2), all.x = T)
    dfo <- dfo[!is.na(dfo$b),]
}

# fit models
{
    hlm <- lm(w ~ b * g, data = dfo[is.na(dfo$type),])        # 1814.3848     0.9588
    alm <- lm(w ~ b * g, data = dfo)                         # 236.2960      0.9949
    
    dfo$h.fv[is.na(dfo$type)] <- hlm$fitted.values
    dfo$h.res[is.na(dfo$type)] <- hlm$residuals
    
    dfo$a.fv <- alm$fitted.values
    dfo$a.res <- alm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(alm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(alm)$sigma, 2)),
          paste0(fpath, "fitted-wv-all-refurb.txt"))
    write(paste0("Adj. $r^2$ ", round(summary(hlm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(hlm)$sigma, 2)),
          paste0(fpath, "fitted-wv-healthy-refurb.txt"))
    write(paste0(sum(abs(hlm$residuals) > 2 * sd(df$h.fv, na.rm = T)), " px > ", round(2 * sd(df$h.fv, na.rm = T), 0), " res"),
          paste0(fpath, "fitted-wv-res-refurb.txt"))
}

pdf(paste0(fpath, "fitted-wv-all-px-refurb.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$w, dfo$a.fv, xlim = c(0,65535), ylim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed white value", ylab = "Fitted white value")
    abline(line(dfo$w, dfo$a.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
}

pdf(paste0(fpath, "fitted-wv-healthy-px-refurb.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$w, dfo$h.fv, xlim = c(0,65535), ylim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed white value", ylab = "Fitted white value")
    abline(line(dfo$w, dfo$h.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    dev.off()
}

pdf(paste0(fpath, "fitted-wv-healthy-res-refurb.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(dfo$h.fv, dfo$h.res, xlim = c(0,65535), ylim = c(-4000,4000),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Fitted white value", ylab = "Residual")
    abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    abline(h = 1200 * c(-1,1), lty = 2, col = adjustcolor("darkblue", alpha = 0.2))
    dev.off()
}

####################################################################################################

# PREDICT SC FROM GREY/BLACK IMAGES                                                             ####

sc <- array(dim = c(2048, 2048))
sc[3:1998, 33:2028] <- readRDS("./Other-data/Shading-corrections.rds")[,,"160430"]
sc[is.infinite(sc) | is.na(sc)] <- 0
# where pixel is hot or dead in more than one image, set to 0

# data frame of all variables for active region of image
{
    df <- setNames(data.frame(melt(acq[,,"black"]), 
                                    melt(acq[,,"grey"]),
                                    melt(acq[,,"white"]),
                                    melt(sc))[,c("X1", "X2", "value", "value.1", "value.2", "value.3")],
                         nm = c("x", "y", "b", "g", "w", "sc"))
    df <- merge(df, px, by = c(1,2), all.x = T)
    df <- df[!is.na(df$b),]
    
    # fit linear model to healthy px only, check line through fitted/actual
    hlm <- lm(sc ~ b * g, data = df[is.na(df$type),])        # 167.6664      0.9964
    alm <- lm(sc ~ b * g, data = df)                         # 236.2960      0.9949
    
    df$h.fv[is.na(df$type)] <- hlm$fitted.values
    df$h.res[is.na(df$type)] <- hlm$residuals
    
    df$a.fv <- alm$fitted.values
    df$a.res <- alm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(alm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(alm)$sigma, 2)),
          paste0(fpath, "fitted-sc-all.txt"))
    write(paste0("Adj. $r^2$ ", round(summary(hlm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(hlm)$sigma, 2)),
          paste0(fpath, "fitted-sc-healthy.txt"))
    write(paste0(sum(abs(hlm$residuals) > 2 * sd(df$h.fv, na.rm = T)), " px > ", round(2 * sd(df$h.fv, na.rm = T), 0), " res"),
          paste0(fpath, "fitted-sc-res.txt"))
}
    
pdf(paste0(fpath, "fitted-sc-all-px.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(df$sc, df$a.fv, 
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed value", ylab = "Fitted value")
    abline(line(df$sc, df$a.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
}

pdf(paste0(fpath, "fitted-sc-healthy-px.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(df$sc, df$h.fv, asp = T,
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed value", ylab = "Fitted value")
    abline(line(df$sc, df$h.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    dev.off()
}

pdf(paste0(fpath, "fitted-sc-healthy-res.pdf")); {
    par(mar = c(4, 4, 1, 1))
    smoothScatter(df$h.fv, df$h.res,
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Fitted value", ylab = "Residual")
    abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
    abline(h = 500 * c(-1,1), col = adjustcolor("darkblue", alpha = 0.4), lty = 2)  
    dev.off()
}

####################################################################################################

# PREDICT BLACK FROM WHITE/GREY IMAGES                                                          ####

# does process work in reverse? Suspect not
fpath <- "./Notes/WV-prediction/fig/"

{
    # data frame of all variables for active region of image
    df <- merge(setNames(data.frame(melt(acq[,,"black"]), 
                                    melt(acq[,,"grey"]),
                                    melt(acq[,,"white"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                         nm = c("x", "y", "b", "g", "w")),
                px, by = c(1:2), all.x = T)
    df <- df[!is.na(df$b),]
    
    # fit linear model to healthy px only, check line through fitted/actual
    hlm <- lm(b ~ w * g, data = df[is.na(df$type),])        # 5458      -0.0001347
    alm <- lm(b ~ w * g, data = df)                         # 5458      -0.0001347
    
    df$h.fv[is.na(df$type)] <- hlm$fitted.values
    df$h.res[is.na(df$type)] <- hlm$residuals
    
    df$a.fv <- alm$fitted.values
    df$a.res <- alm$residuals
    
    write(paste0("Adj. $r^2$ ", round(summary(alm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(alm)$sigma, 2)),
          paste0(fpath, "fitted-bv-all.txt"))
    write(paste0("Adj. $r^2$ ", round(summary(hlm)$adj.r.squared, 3), "; ",
                 "residual SD ", round(summary(hlm)$sigma, 2)),
          paste0(fpath, "fitted-bv-healthy.txt"))
    write(paste0(sum(abs(hlm$residuals) > 2 * sd(df$h.fv, na.rm = T)), " px > ", round(2 * sd(df$h.fv, na.rm = T), 0), " res"),
          paste0(fpath, "fitted-bv-res.txt"))
    
    pdf(paste0(fpath, "fitted-bv-all-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$b, df$a.fv, xlim = c(0,65535), ylim = c(0,65535),
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Observed black value", ylab = "Fitted black value")
        abline(line(df$b, df$a.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
        dev.off()
    }
    
    pdf(paste0(fpath, "fitted-bv-healthy-px.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$b, df$h.fv, xlim = c(0,65535), ylim = c(0,65535),
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Observed black value", ylab = "Fitted black value")
        abline(line(df$b, df$h.fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
        dev.off()
    }
    
    pdf(paste0(fpath, "fitted-bv-healthy-res.pdf")); {
        par(mar = c(4, 4, 1, 1))
        smoothScatter(df$h.fv, df$h.res, 
                      colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                      xlab = "Fitted value", ylab = "Residual")
        abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)  
        abline(h = 1200 * c(-1,1), col = adjustcolor("darkblue", alpha = 0.4), lty = 2)  
        dev.off()
    }
    
    df.xt <- df[which(abs(df$h.res) > 2 * sd(df$h.fv, na.rm = T)), ]
    nrow(df.xt)
    # only 39px
    
    plot(df.xt[,1:2], pch = 20, col = c("red", "black")[(df.xt$h.res > 0) + 1])
    
    df.xt$md.b <- md[,,"black"][as.matrix(df.xt[,1:2])]
    df.xt$md.g <- md[,,"grey"][as.matrix(df.xt[,1:2])]
}



