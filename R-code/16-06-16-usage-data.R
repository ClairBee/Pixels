
xtek <- read.csv("./Other-data/xtek-usage.csv", row.names = 1, stringsAsFactors = F)
bp <- readRDS("./Notes/Spatial/fig/bad-px.rds")

# plots to output:
    # usage vs pixel damage (cumulative)
    # moving average of usage vs new damaged pixels?
    
    # usage vs screen spots (non-cumulative: label lines with # cells covered by spots)

# extract & rename only columns of interest
{
    xt <- data.frame(acq.date = strptime(lapply(sapply(lapply(sapply(xtek$filenm, strsplit, "[[]"), 
                                                              "[", 2), strsplit, "[]]"), "[", 1),
                                         "%Y-%m-%d %H.%M.%S"),
                     kV = xtek$XraySettings.Settings.kV,
                     uA = xtek$XraySettings.Settings.uA,
                     exp.time = xtek$ImagingSettings..attrs.exposure,
                     n.proj = xtek$Projections,
                     frames.per.proj = xtek$FramesPerProjection,
                     filter.mat = xtek$XrayFilterMaterial,
                     filter.thk = xtek$XrayFilterThickness,
                     filenm = xtek$filenm,
                     stringsAsFactors = F)
    
    xt <- xt[order(xt$acq.date),]
    xt <- xt[231:540,]
    
    xt$usage <- (xt$kV * xt$uA * xt$exp.time / 1000 * xt$n.proj * xt$frames.per.proj) / 1000000
}

fpath <- "./Notes/Usage/fig/"

# function to scrape acquisition dates in correct format
acq.datetimes <- function() {
    ff <- list.dirs("./Image-data", recursive = F, full.names = T)
    ff <- ff[ff != "./Image-data/150702"]     # corrupted files
    
    ff <- unlist(lapply(ff, paste, "/black", sep = ""))
    
    as.POSIXct(unlist(lapply(ff, 
                             function (x) xmlToList(xmlParse(list.files(x, pattern = "\\.xml", full.names = T)[1]))$TimeAcquired)))
}

# add IO acquisitions to data frame & resort
{
    xt <- rbind(xt,
                data.frame(acq.date = acq.datetimes(),
                           kV = 85,
                           uA = 100,
                           exp.time = 1000,
                           n.proj = 20,
                           frames.per.proj = 1,
                           filter.mat = "None",
                           filter.thk = 0,
                           filenm = names(bp),
                           usage = (85*100*20) / 1000000))
    
    xt <- xt[order(xt$acq.date),]
}
xt$cum.usage <- cumsum(xt$usage)
           
####################################################################################################

# USAGE OVER TIME                                                                               ####

usage.lm <- lm(cum.usage ~ acq.date, xt)
usage.qm <- lm(cum.usage ~ poly(acq.date, 2), xt)

pdf(paste0(fpath, "usage-vs-time.pdf"), height = 4, width = 7); {
    par(mar = c(2, 3, 1, 1))
    plot(xt$acq.date, xt$usage, type = "o", pch = 20,
         xlab = "", ylab = "Usage (arbitrary units)")
    abline(v = acq.datetimes(), col = "red", lty = 3)
    
    lines(xt$acq.date, lm(usage ~ poly(acq.date, 2), xt)$fitted.values, col = "cyan3")
    lines(xt$acq.date, lm(usage ~ acq.date, xt)$fitted.values, col = "gold")
    dev.off()
}

pdf(paste0(fpath, "cum-usage-vs-time.pdf"), height = 4, width = 7); {
    par(mar = c(2, 3, 1, 1))
    plot(xt$acq.date, cumsum(xt$usage), type = "o", pch = 20,
         xlab = "", ylab = "Usage (arbitrary units)")
    abline(v = acq.datetimes(), col = "red", lty = 3)
    
    lines(xt$acq.date, usage.qm$fitted.values, col = "cyan3")
    lines(xt$acq.date, usage.lm$fitted.values, col = "gold")
    
    dev.off()
}

xt$ma.usage <- filter(xt$usage, filter = rep(0.1, 10), sides = 1)

pdf(paste0(fpath, "ma-usage-vs-time.pdf"), height = 4, width = 7); {
    par(mar = c(2, 3, 1, 1))
    plot(xt$acq.date, xt$ma.usage, type = "o", pch = 20,
         xlab = "", ylab = "Usage (arbitrary units)")
    abline(v = acq.datetimes(), col = "red", lty = 3)
    
    dev.off()
}

####################################################################################################

# USAGE VS DAMAGE                                                                               ####

plot(acq.datetimes(), unlist(lapply(bp, nrow)) / 1996^2 * 100, type = "o", pch = 20,
     xlab = "", ylab = "Proportion of pixels defective")

pdf(paste0(fpath, "cum-usage-vs-bp-prop.pdf"), height = 4, width = 7); {
    par(mar = c(2,4,1,1))
    plot(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"],
         unlist(lapply(bp, nrow)) / 1996^2 * 100,
         type = "o", pch = 20, xlab = "Cumulative usage", ylab = "Percentage of pixels defective")
    lines(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"],
          lm(unlist(lapply(bp, nrow)) / 1996^2 * 100 ~ xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], xt)$fitted.values,
          col = "gold")
    lines(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"],
          lm(unlist(lapply(bp, nrow)) / 1996^2 * 100 ~ poly(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 2), xt)$fitted.values,
          col = "cyan3")
    dev.off()
}


# plot each type of defective pixel in turn
{
    pdf(paste0(fpath, "cum-usage-vs-bp-prop-by-type.pdf"), height = 4, width = 7); {
        par(mar = c(2,4,1,1))
        plot(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"],
             unlist(lapply(bp, function(x) sum(x$type == "l.bright"))) / 1996^2 * 100, 
             col = "gold", ylim = c(0, 0.5), type = "o", pch = 20,
             xlab = "Cumulative usage", ylab = "Percentage of pixels defective")
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "line.b"))) / 1996^2 * 100, 
               col = "orange", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "hot"))) / 1996^2 * 100, 
               col = "magenta3", ylim = c(0, 0.01), type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "v.bright"))) / 1996^2 * 100, 
               col = "red", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "bright"))) / 1996^2 * 100, 
               col = "orange", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "l.dim"))) / 1996^2 * 100, 
               col = "green3", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "no.resp"))) / 1996^2 * 100, 
               col = "purple", type = "o", pch = 20)
        dev.off()
    }
    
    pdf(paste0(fpath, "cum-usage-vs-bp-prop-by-type-zoom.pdf"), height = 4, width = 7); {
        par(mar = c(2,4,1,1))
        
        plot(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "hot"))) / 1996^2 * 100, 
               col = "magenta3", ylim = c(0, 0.01), type = "o", pch = 20,
             xlab = "Cumulative usage", ylab = "Percentage of pixels defective")
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "v.bright"))) / 1996^2 * 100, 
               col = "red", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "bright"))) / 1996^2 * 100, 
               col = "orange", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "l.dim"))) / 1996^2 * 100, 
               col = "green3", type = "o", pch = 20)
        
        points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], 
               unlist(lapply(bp, function(x) sum(x$type == "no.resp"))) / 1996^2 * 100, 
               col = "purple", type = "o", pch = 20)
        dev.off()
    }
}

# anything unusual in image mean/variance relationship?
{
    plot(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.m[,,"black",], 3, mean), type = "o", pch = 20, ylim = c(0,65535))
    points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.m[,,"grey",], 3, mean), type = "o", col = "grey", pch = 20)
    points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.m[,,"white",], 3, mean), type = "o", col = "gold", pch = 20)
    
    plot(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.sd[,,"white",], 3, mean), type = "o", col = "gold", pch = 20, ylim = c(0,350))
    points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.sd[,,"grey",], 3, mean), type = "o", col = "grey", pch = 20)
    points(xt[xt$acq.date %in% acq.datetimes(), "cum.usage"], apply(pw.sd[,,"black",], 3, mean), type = "o", pch = 20)
    
    plot(apply(pw.m[,,"black",], 3, mean), apply(pw.sd[,,"black",], 3, mean), pch = 21,
         bg = c(rep("black", 10), "red", "black"))
    plot(apply(pw.m[,,"grey",], 3, mean), apply(pw.sd[,,"grey",], 3, mean), pch = 21,
         bg = c(rep("black", 10), "red", "black"))
    plot(apply(pw.m[,,"white",], 3, mean), apply(pw.sd[,,"white",], 3, mean), pch = 21,
         bg = c(rep("black", 10), "red", "black"))
}


df <- xt[xt$acq.date %in% acq.datetimes(),]
df$n.pix <- unlist(lapply(bp, nrow))
df$pix.pc <- df$n.pix / 1996^2 * 100

pdf(paste0(fpath, "usage-vs-defect-prop.pdf"), height = 4, width = 7); {
    par(mar = c(2,4,1,4))
    plot(xt$acq.date, xt$usage, type = "o", pch = 20, ylim = c(0,1600),
         xlab = "", ylab = "Usage (arbitrary units)")
    
    abline(v = acq.datetimes(), col = "red", lty = 3)
    points(df$acq.date, df$pix.pc * 3000,
           pch = 21, bg = "skyblue", col = "blue")
    
    axis(4, at = c(0:5)*300, labels = c(0:5)*300 / 3000, col = "blue")
    mtext(side = 4, line = 3, "Percentage of detector covered", col = "blue")
    dev.off()
}

pdf(paste0(fpath, "defects-vs-usage.pdf"), height = 4, width = 7); {
    par(mar = c(2,4,1,4))
    plot(df$cum.usage, df$pix.pc, type = "o", pch = 20,
         xlab = "", ylab = "Percentage of pxiels defective")
    abline(coef(lm(pix.pc ~ cum.usage, df)), col = "gold")
    lines(df$cum.usage, lm(pix.pc ~ poly(cum.usage, 2), df)$fitted.values, col = "cyan3")
    dev.off()
}

pdf(paste0(fpath, "defects-vs-time.pdf"), height = 4, width = 7); {
    par(mar = c(2,4,1,4))
    plot(df$acq.date, df$pix.pc, type = "o", pch = 20,
         xlab = "", ylab = "Percentage of pxiels defective")
    abline(coef(lm(pix.pc ~ acq.date, df)), col = "gold")
    lines(df$acq.date, lm(pix.pc ~ poly(acq.date, 2), df)$fitted.values, col = "cyan3")
    dev.off()
}

usage.qm <- lm(pix.pc ~ poly(cum.usage, 2), df)


# what about fitting line to first 10 images only?
{
    plot(df$cum.usage, df$pix.pc,
         type = "o", pch = 20, xlab = "Cumulative usage", ylab = "Percentage of pixels defective")
    abline(coef(lm(pix.pc ~ cum.usage, df[1:10,])), col = "gold", lty = 2)
    abline(coef(lm(pix.pc ~ cum.usage, df)), col = "orange", lty = 2)
    
}

####################################################################################################

# SCREEN SPOTS                                                                                  ####

# also repeat for fitting of screen spots

sp <- readRDS("./Notes/Final-classifications/fig/bad-px-screenspots.rds")
pc <- unlist(c(summarise.bpx(sp)))
pc[is.na(pc)] <- 0
pc <- pc / 1996^2 * 100


pdf(paste0(fpath, "usage-vs-screen-spot-prop.pdf"), height = 4, width = 7); {
    par(mar = c(2,4,1,4))
    plot(xt$acq.date, xt$usage, type = "o", pch = 20, ylim = c(0,1600),
         xlab = "", ylab = "Usage (arbitrary units)")
    
    abline(v = acq.datetimes(), col = "red", lty = 3)
    points(xt[xt$acq.date %in% acq.datetimes(), "acq.date"], pc * 3000,
           pch = 21, bg = "skyblue", col = "blue")
    
    axis(4, at = c(0:5)*300, labels = c(0:5)*300 / 3000, col = "blue")
    mtext(side = 4, line = 3, "Percentage of detector covered", col = "blue")
    dev.off()
}

# relationships between explantory variables?
xt$filter.mat[is.na(xt$filter.mat)] <- "None"
xt$filter.mat <- factor(xt$filter.mat)

plot(xt[,2:6], pch = 20, col = c("red", "black", "blue")[xt$filter.mat])


plot(xt$acq.date, xt$usage, pch = 20, 
     col = c(NA, "black", NA)[xt$filter.mat], xlab = "No filter", ylab = "Usage (arbitrary units)")
lines(xt$acq.date, xt$usage)

abline(v = acq.datetimes(), col = "red", lty = 3)
points(xt[xt$acq.date %in% acq.datetimes(), "acq.date"],
       (unlist(lapply(sp, function (x) if (is.null(nrow(x))) {0} else {nrow(x)})) / 1996^2 * 100) * 3000,
       pch = 21, bg = "gold", col = "orange")


plot(sp[[1]][,1:2], pch = 15, cex = 0.4, asp = T, xlim = c(0,1996), ylim = c(0,1996))
points(sp[[2]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("blue", alpha = 0.3))
points(sp[[3]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("green3", alpha = 0.3))
points(sp[[4]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("cyan3", alpha = 0.3))
points(sp[[5]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("skyblue", alpha = 0.3))
points(sp[[6]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("purple", alpha = 0.3))

plot(sp[[7]][,1:2], pch = 15, cex = 0.4, asp = T, xlim = c(0,1996), ylim = c(0,1996), col = adjustcolor("gold", alpha = 0.3))
points(sp[[8]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("gold", alpha = 0.3))
points(sp[[9]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("gold", alpha = 0.3))
points(sp[[10]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("gold", alpha = 0.3))
points(sp[[11]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("orange", alpha = 0.3))

points(sp[[12]][,1:2], pch = 15, cex = 0.4, col = adjustcolor("red", alpha = 0.3))

plot(xt$acq.date, xt$kV, pch = 20, type = "o", xlab = "", ylab = "Usage (arbitrary units)", ylim = c(0,350))
abline(v = acq.datetimes(), col = "red", lty = 3)
points(xt[xt$acq.date %in% acq.datetimes(), "acq.date"],
       (unlist(lapply(sp, function (x) if (is.null(nrow(x))) {0} else {nrow(x)})) / 1996^2 * 100) * 600,
       pch = 21, bg = "gold", col = "orange")

plot(xt$acq.date, xt$uA, pch = 20, type = "o", xlab = "", ylab = "Usage (arbitrary units)", ylim = c(0,650))
abline(v = acq.datetimes(), col = "red", lty = 3)
points(xt[xt$acq.date %in% acq.datetimes(), "acq.date"],
       (unlist(lapply(sp, function (x) if (is.null(nrow(x))) {0} else {nrow(x)})) / 1996^2 * 100) * 1200,
       pch = 21, bg = "gold", col = "orange")

plot(xt$acq.date, xt$exp.time, pch = 20, type = "o", xlab = "", ylab = "Usage (arbitrary units)", ylim = c(0,4000))
abline(v = acq.datetimes(), col = "red", lty = 3)
points(xt[xt$acq.date %in% acq.datetimes(), "acq.date"],
       (unlist(lapply(sp, function (x) if (is.null(nrow(x))) {0} else {nrow(x)})) / 1996^2 * 100) * 8000,
       pch = 21, bg = "gold", col = "orange")

####################################################################################################
