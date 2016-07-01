
# ALTERNATIVE THRESHOLDING APPROACHES

library("IO.Pixels"); library("CB.Misc")

acq <- readRDS("./02_Objects/images/pwm-160430.rds")
md <- readRDS("./02_Objects/med-diffs/md-160430.rds")

####################################################################################################

# SCREEN SPOTS - GREY VS WHITE                                                                  ####

# white image gives clearer view in this case - dips stand out against noise

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

####################################################################################################

# BASIC THRESHOLDING - BLACK & GREY                                                             ####

Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "skyblue", "green3")

px <- rbind(data.frame(edge.px(acq, edge.width = 40), type = "edge"),
            data.frame(no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]), type = "no.resp"),
            data.frame(which(acq[,,"black"] == 65535, arr.ind = T), type = "hot"),
            data.frame(which(acq[,,"grey"] == 0, arr.ind = T), type = "dead"),
            data.frame(screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40), type = "screen.spot"),
            data.frame(which(find.lines(acq[, , "black"], midline = 1024.5) > 0, arr.ind = T), type = "line.b"),
            data.frame(which()))


rbind(data.frame(which(md.b[[dt]] > 2 * mad(pw.m[,,"black", dt]), arr.ind = T),
                 type = "l.bright"),
      data.frame(which(md.b[[dt]] < -2 * mad(pw.m[,,"black", dt]), arr.ind = T),
                 type = "l.dim"),
      data.frame(which(md.g[[dt]] > 2 * mad(pw.m[,,"grey", dt]), arr.ind = T),
                 type = "l.bright"),
      data.frame(which(md.g[[dt]] < -2 * mad(pw.m[,,"grey", dt]), arr.ind = T),
                 type = "l.dim")))

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

legend("bottomright", col = c("orange", "green3"), lty = 2, 
       legend = c("1 x MAD", "1 x SD"))
# plots
{
    plot(which(md[,,"black"] > 2 * mad(acq[,,"black"], na.rm = T), arr.ind = T), pch = 16, cex = 0.2,
         col = "cyan3")
    points(which(md[,,"black"] > 2 * sd(acq[,,"black"], na.rm = T), arr.ind = T), pch = 16, cex = 0.2,
           col = "red")
}

