
library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Line-detection/Column-427/"

load.pixel.means()

# last 2px at 'inner' end are artefact of smoothing by 5-square kernel & should be removed

focal.plot <- function(im, centre, surround = 5, dp = 1, scale.by = 1000, pt.cex = 0.7, bad.px, bpx.cex = 2.5, ...) {
    
    ff <- get.focus(data.frame(x = centre[1], y = centre[2]), surround = surround)
    
    pixel.image(im, xlim = range(ff[,1]), ylim = range(ff[,2]), ...)
    text(ff, labels = round(im[ff]/scale.by, dp), cex = pt.cex)
    
    if (!missing(bad.px)) {
        points(bad.px[,1:2], pch = 0, cex = bpx.cex)
    }
}

####################################################################################################

# 427 - INITIAL PLOTS                                                                           ####

# line transect vs neighbours
{
    o.plot(pw.m[427,,"black", "160430"], xlim = c(992, 1996), ylim = c(4500,8000),
           main = "Bright column 427 & neighbours", xlab = "", ylab = "")
    lines(pw.m[426,,"black", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[428,,"black", "160430"], col = adjustcolor("orange", alpha = 0.4))
    abline(v = 992.5, col = "blue")
    legend("topleft", col = c("cyan3", "black", "orange"), lty = 1, legend = c(426:428), cex = 0.7)
}

{
    o.plot(pw.m[427,,"black", "160430"], xlim = c(992, 1996), ylim = c(4500,8000),
           main = "Bright column 427 & neighbours", xlab = "", ylab = "")
    lines(pw.m[426,,"black", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[428,,"black", "160430"], col = adjustcolor("orange", alpha = 0.4))
    abline(v = 992.5, col = "blue")
    legend("topleft", col = c("cyan3", "black", "orange"), lty = 1, legend = c(426:428), cex = 0.7)
}
    

# track difference between neighours all the way along column
{
    plot(pw.m[427, , "black", "160314"] - pw.m[426, , "black", "160314"], type = "l",
          xlim = c(1180, 1996), ylim = c(-200,1000))
    
    diff1 <- findInterval(pw.m[427, , "black", "160314"] - pw.m[426, , "black", "160314"], c(300))
    diff2 <- findInterval(pw.m[427, , "black", "160314"] - pw.m[428, , "black", "160314"], c(300))
    diff <- (diff1 + diff2 > 0) + 1
    
    points(pw.m[427, , "black", "160314"] - pw.m[426, , "black", "160314"],
           col = c("red", "gold")[diff], pch = 20)
    
    
    abline(h = 300, col = "red")
    points(pw.m[427, , "black", "160314"] - pw.m[426, , "black", "160314"], pch = 20,
           col = c("black", "red")[(pw.m[427, 993:1996, "black", "160314"] - pw.m[426, 993:1996, "black", "160314"] > 300)+1])
    
    which(pw.m[427, 993:1996, "black", "160430"] - pw.m[426, 993:1996, "black", "160430"] < 300)
}
# line transect: development over 4 most recent acquisitions
{
    o.plot(pw.m[427,,"black", "160430"], xlim = c(992, 1996), ylim = c(4500,8000),
           main = "Development of bright column 427", xlab = "", ylab = "")
    lines(pw.m[427,,"black", "160314"], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[427,,"black", "151015"], col = adjustcolor("orange", alpha = 0.4))
    lines(pw.m[427,,"black", "150828"], col = adjustcolor("gold", alpha = 0.4))
    legend("topleft", col = c("black", "cyan3", "orange", "gold"), lty = 1, 
           legend = sapply(dimnames(pw.m)[[4]][12:9], fancy.date), cex = 0.7)
}

# plot area at end of bright line (row 1200, with 5px radius)
pdf(paste0(fpath, "col-427-end-plots.pdf")); {
    par(mar = c(2, 3, 2, 2), mfrow = c(4, 3))
    
    lapply(dimnames(pw.m)[[4]],
           function(y) lapply(dimnames(pw.m)[[3]],
                              function(x) focal.plot(pw.m[ , , x, y], 
                                                     centre = c(427, 1200),
                                                     bad.px = bp[[y]],
                                                     pt.cex = 0.5, cex.main = 0.8,
                                                     title = paste0(fancy.date(y), ", ", x))))
    dev.off()
}

focal.plot(pw.m[ , , "black", "160430"], centre = c(427, 1200), bad.px = bp$"160430", pt.cex = 0.5, cex.main = 0.8)

####################################################################################################

# 809 - INITIAL PLOTS                                                                           ####

o.plot(pw.m[809,,"black", "160430"], xlim = c(0,992), ylim = c(4000,6000))
lines(pw.m[808,,"black", "160430"], col = "cyan3")
lines(pw.m[810,,"black", "160430"], col = "orange")

# plot area at end of bright line (row 1200, with 5px radius)
pdf(paste0(fpath, "col-809-end-plots.pdf")); {
    par(mar = c(2, 3, 2, 2), mfrow = c(4, 3))
    
    lapply(dimnames(pw.m)[[4]],
           function(y) lapply(dimnames(pw.m)[[3]],
                              function(x) focal.plot(pw.m[ , , x, y], 
                                                     centre = c(809, 177),
                                                     bad.px = bp[[y]],
                                                     pt.cex = 0.5, cex.main = 0.8,
                                                     title = paste0(fancy.date(y), ", ", x))))
    dev.off()
}

####################################################################################################

# UPDATE FUNCTION TO RETAIN LINE SEGMENTS                                                       ####

# convolve & threshold
conv <- threshold(convolve.lines(pw.m[,,"black", "160430"], k.size = 5), level = 5500)
image(1:1996, 1:1996, conv, col = c(NA, "red"), xlim = c(400, 900))

# smooth & repair lines
sm <- smooth.lines(conv, sm.size = 11, min.length = 6)
image(1:1996, 1:1996, sm, col = c(NA, "cyan3"), xlim = c(400, 900))


#        ---~~~===###===~~~--- FINE UP TO THIS POINT ---~~~===###===~~~---


# NOT FINE.
{
    # filter out short line segments and return image of long segments identified
    lines <- filter.lines(sm, filter.at = list(cover = 0.5, filled = 20), midline = 992.5)
    image(1:1996, 1:1996, lines, col = c(NA, "green3", "green3"), xlim = c(400, 900))
}
# need to adjust filtering algorithm so that line segments aren't automatically merged
# problem is in line summary: everything collapses to single column


# adjust line summary algorithm to retain individual segments
{
    # convert to raster & clump to get line IDs
    lines <- clump(m2r(sm), dir = 4)
    
    xy <- data.frame(xyFromCell(lines, which(!is.na(getValues(lines)))),
                     id = getValues(lines)[!is.na(getValues(lines))])
    
    # summarise line segments
    ll <- rbind(ddply(xy[xy$y > 992.5,], .(id, col = x), summarise, panel = "U",
                      ymin = min(y), ymax = max(y), length = length(x)),
                ddply(xy[xy$y < 992.5,], .(id, col = x), summarise, panel = "L",
                      ymin = min(y), ymax = max(y), length = length(x)))
    
    ### filtering ###
    
    # remove line segments within 10 columns of panel edge
    ll <- ll[ll$col > 10 & ll$col < 1996-10,]
    
    # remove line segments less than 10 in length
    # relate to kernel parameters - what is likely result of convolution & smoothing?
    ll <- ll[ll$length > 10,]
    
    # adjust line ends to panel midline/edge if within 10px (truncated by smoothing)
    # relate this to smoothing parameters - what is distance likely to be?
    ll$ymin[ll$ymin %in% (ceiling(992.5) + c(0:9))] <- ceiling(992.5)
    ll$ymin[ll$ymin %in% (1 + c(0:9))] <- 1
    
    ll$ymax[ll$ymax %in% (floor(992.5) - c(0:9))] <- floor(992.5)
    ll$ymax[ll$ymax %in% (1996 - c(0:9))] <- 1996
    
    # recalculate lengths according to new line ends
    ll$length <- ll$ymax - ll$ymin + 1
    
    # remove any segments less than 20px in length - UNLESS associated with longer segment
    long <- ll[ll$length > 20,]
    
    long <- ll[ll$col %in% long$col,]
}

# transfer line segments to image array
{
    new.im <- array(0, dim = c(1996, 1996))
    
    if (nrow(long) > 0) {
        for (i in 1:nrow(long)) {
            new.im[long$col[i], long$ymin[i]:long$ymax[i]] <- i
        }
    }
    
}

# now adapted into 'official' functions
ll <- summarise.lines(sm)
ff <- filter.lines(sm)
image(1:1996, 1:1996, ff, col = c(NA, "red", "red", "blue", "blue"), xlim = c(400, 900))
draw.panels()
table(c(ff))

qq <- find.lines(pw.m[,,"black", "160314"])
image(1:1996, 1:1996, qq, col = c(NA, "red", "red", "blue", "blue"), xlim = c(400, 900))

o.plot(pw.m[427,,"black", "160314"], xlim = c(993, 1996), ylim = c(4000,7000))
lines(pw.m[426,,"black", "160314"], col = adjustcolor("skyblue", alpha = 0.4))
lines(pw.m[428,,"black", "160314"], col = adjustcolor("orange", alpha = 0.4))
lines((qq[427,] > 0) * 4500, col = "red")

o.plot(pw.m[427,,"black", "160430"], xlim = c(993, 1996), ylim = c(4000,7000))
lines(pw.m[426,,"black", "160430"], col = adjustcolor("skyblue", alpha = 0.4))
lines(pw.m[428,,"black", "160430"], col = adjustcolor("orange", alpha = 0.4))
lines((qq[427,] > 0) * 4500, col = "red")

focal.plot(pw.m[,,"black", "160314"], centre = c(427, 1820), surround = 20)
focal.plot(pw.m[,,"black", "160430"], centre = c(427, 1820), surround = 20)

plot((qq[427,] > 0) * 4500, col = "red", pch = 16, cex = 0.4, xlim = c(1800, 1850))

# some breaks in upper line, but only when first identified.
# (due to neighbouring brighter pixels)

