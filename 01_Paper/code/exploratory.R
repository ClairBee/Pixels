
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

load.pixel.means.2048()

####################################################################################################

# SCREEN SPOTS                                                                                  ####

clump.centres <- function(px) {
    
    # clump adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = c(2048, 2048))), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    df <- ddply(xy, .(id), summarise,
                xm = mean(x), ym = mean(y),
                r = ceiling(max(max(x) - min(x), max(y) - min(y)) / 2))
    df
}

# get centres of screen spots
spots <- screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40)
ss <- clump.centres(spots)

# white image with screen spots manually marked
jpeg(paste0(fpath, "pwm-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(pw.m[,,"white", "141009"])
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# shading-corrected image with screen spots manually marked
jpeg(paste0(fpath, "sc-image-141009-white.jpg")); {
    par(mar = c(2,2,1,1))
    pixel.image(shading.corrected(pw.m[,,,"141009"]))
    symbols(ss$x, ss$y, circles = 2 * ss$r, add = T, inches = F)
    dev.off()
}

# position of screen spots in successive acquisitions overplotted
jpeg(paste0(fpath, "spots-overplotted.jpg")); {
    par(mar = c(2,2,1,1))
    plot(screen.spots(pw.m[,,"white", "141009"], enlarge = T, ignore.edges = 40),
         col = "cyan3", pch = ".", xlab = "", ylab = "")
    points(screen.spots(pw.m[,,"white", "141118"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("gold", alpha = 0.4), pch = ".")
    points(screen.spots(pw.m[,,"white", "141217"], enlarge = T, ignore.edges = 40),
           col = adjustcolor("magenta3", alpha = 0.4), pch = ".")
    legend("topright", pch = 15, col = c("cyan3", "gold", "magenta3"), bty = "n",
           legend = sapply(dimnames(pw.m)[[4]][1:3], fancy.date))
    dev.off()
}

####################################################################################################
