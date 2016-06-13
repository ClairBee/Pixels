
library("IO.Pixels"); library("CB.Misc")

load.pixel.means()

focal.plot <- function(im, centre, surround = 5, dp = 1, scale.by = 1000, pt.cex = 0.7, bad.px, bpx.cex = 2.5, ...) {
    
    ff <- get.focus(data.frame(x = centre[1], y = centre[2]), surround = surround)
    
    pixel.image(im, xlim = range(ff[,1]), ylim = range(ff[,2]), ...)
    text(ff, labels = round(im[ff]/scale.by, dp), cex = pt.cex)
    
    if (!missing(bad.px)) {
        points(bad.px[,1:2], pch = 0, cex = bpx.cex)
    }
}
####################################################################################################

# line transect vs neighbours
{
    o.plot(pw.m[427,,"black", "160430"], xlim = c(992, 1996), ylim = c(4500,8000),
           main = "Bright column 427 & neighbours", xlab = "", ylab = "")
    lines(pw.m[426,,"black", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[428,,"black", "160430"], col = adjustcolor("orange", alpha = 0.4))
    abline(v = 992.5, col = "blue")
    legend("topleft", col = c("cyan3", "black", "orange"), lty = 1, legend = c(426:428), cex = 0.7)
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

