
# START FROM RAW IMAGES, CLASSIFY EVERYTHING FROM SCRATCH

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Med-diff-classification/fig/"

cat <- c("no response", "dead", "hot", 
         "v.bright", "bright", "line.b", "l.bright", "s.bright", 
         "screen spot", "line.d", "edge",
         "v.dim", "dim", "l.dim", "s.dim")

cat.cols <- c("purple", "black", "magenta3", 
              "red", "orange", "gold", "gold", NA,
              "grey", "violet", NA, 
              "green3", "green", "lightskyblue", NA)

####################################################################################################

# revised process to omit parametric model: instead, use difference from local median

####################################################################################################

# then repeat the process over the old images
# clean up functions: should be able to apply to old panel structure too!

# should also be able to specify degree of edge cropping

####################################################################################################

# FIND GLOBALLY EXTREME VALUES                                                                  ####
load.pixel.means()
dt <- "160430"

# quick function to get thresholds
qq <- get.dim.bright.px(pw.m[,,"white", dt])

# hot, v.bright, v.dim, dead, no response
bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                         rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                ordered("edge", levels = cat)),
                     c("row", "col", "type")),
            data.frame(no.response(dt), type = "no response"),
            data.frame(which(pw.m[, , "black", dt] == 65535, arr.ind = T), type = "hot"),
            data.frame(which(pw.m[, , "white", dt] == 0, arr.ind = T), type = "dead"),
            screen.spots.xy(dt),
            get.dim.bright.px(pw.m[,,"white", dt]),
            get.dim.bright.px(pw.m[,,"grey", dt]),
            get.dim.bright.px(pw.m[,,"black", dt]))

bp <- bp[order(bp$type),]
bp <- bp[!duplicated(bp[,1:2]),]
bp <- bp[!bp$type %in% c("s.bright", "s.dim"),]     # retain only more severe categories

table(bp$type)

plot.bad.px(bp)

# compare classification in each of these images
{
    db <- lapply(list(black = get.dim.bright.px(pw.m[,,"black", dt]),
                      grey = get.dim.bright.px(pw.m[,,"grey", dt]),
                      white = get.dim.bright.px(pw.m[,,"white", dt])),
                 function (x) x[!x$type %in% c("s.bright", "s.dim"),])
    
    db2 <- merge(merge(db$black, db$grey, by = c("row", "col"), suffix = c(".black", ".grey"), all = T),
                 db$white, by = c("row", "col"), suffix = c("", ".white"), all = T)
    
    table("black" = db2$type.black, "grey" = db2$type.grey, useNA = "ifany")
    {
        #                                   grey
        #    black      v.bright bright s.bright v.dim dim s.dim <NA>
        #    v.bright       194      0        0     0   0     0    0
        #    bright          19     67        0     0   0     0    2
        #    s.bright         0      0        0     0   0     0    0
        #    v.dim            0      0        0     5   0     0    0
        #    dim              0      0        0     0   0     0    1
        #    s.dim            0      0        0     0   0     0    0
        #    <NA>             0     49        0   181 467     0  652
    }

    table("black" = db2$type.black, "white" = db2$type, useNA = "ifany")
    {
           #                                  white
           #    black      v.bright bright s.bright v.dim dim s.dim <NA>
           #    v.bright       193      1        0     0   0     0    0
           #    bright          80      7        0     0   0     0    1
           #    s.bright         0      0        0     0   0     0    0
           #    v.dim            0      0        0     4   1     0    0
           #    dim              0      0        0     0   0     0    1
           #    s.dim            0      0        0     0   0     0    0
           #    <NA>            51    276        0   395 626     0    1
       }

    table("grey" = db2$type.grey, "white" = db2$type, useNA = "ifany")
    {
        #                                    white
        #   grey       v.bright bright s.bright v.dim dim s.dim <NA>
        #   v.bright       212      1        0     0   0     0    0
        #   bright          95     19        0     0   0     0    2
        #   s.bright         0      0        0     0   0     0    0
        #   v.dim            0      0        0   185   1     0    0
        #   dim              0      0        0   213 254     0    0
        #   s.dim            0      0        0     0   0     0    0
        #   <NA>            17    264        0     1 372     0    1
    }    
}

# plots of thresholds
{
    lim <- list(b = bad.px.limits(pw.m[,,"black", dt]),
                g = bad.px.limits(pw.m[,,"grey", dt]),
                w = bad.px.limits(pw.m[,,"white", dt]))
    pdf(paste0(fpath, "threshold-histograms.pdf"), height = 4, width = 7); {
        par(mar = c(2, 2, 1, 1))
        hist(pw.m[,,"black", dt], breaks = "fd", col = "darkgreen", border = "green3", ylim = c(0,250), xlab = "", ylab = "", main = ""); {
            rect(0, 0, lim$b$dv, 600, col = adjustcolor("skyblue", alpha = 0.5), border = NA)
            rect(lim$b$dv, 0, lim$b$dm, 600, col = adjustcolor("cyan3", alpha = 0.5), border = NA)
            rect(lim$b$bm, 0, lim$b$bv, 600, col = adjustcolor("gold", alpha = 0.5), border = NA)
            rect(lim$b$bv, 0, 65535, 600, col = adjustcolor("orange", alpha = 0.5), border = NA)
            lines(c(0,0), c(0,length(which(pw.m[,,"white", dt] == 0))), col = "blue", lwd = 4)
            lines(c(65535,65535), c(0,length(which(pw.m[,,"black", dt] == 65535))), col = "red", lwd = 4)
            hist(pw.m[11:1985,11:1985,"black", dt], breaks = "fd", add = T, col = "black")
        }
        hist(pw.m[,,"grey", dt], breaks = "fd", col = "darkgreen", border = "green3", ylim = c(0,250), xlab = "", ylab = "", main = ""); {
            rect(0, 0, lim$g$dv, 600, col = adjustcolor("skyblue", alpha = 0.5), border = NA)
            rect(lim$g$dv, 0, lim$g$dm, 600, col = adjustcolor("cyan3", alpha = 0.5), border = NA)
            rect(lim$g$bm, 0, lim$g$bv, 600, col = adjustcolor("gold", alpha = 0.5), border = NA)
            rect(lim$g$bv, 0, 65535, 600, col = adjustcolor("orange", alpha = 0.5), border = NA)
            rect(qJohnson(0.001, JohnsonFit(pw.m[,,"black", dt])), 0, 
                 qJohnson(0.999, JohnsonFit(pw.m[,,"black", dt])), 600, col = adjustcolor("slateblue1", alpha = 0.5), border = NA)
            lines(c(0,0), c(0,length(which(pw.m[,,"white", dt] == 0))), col = "blue", lwd = 4)
            lines(c(65535,65535), c(0,length(which(pw.m[,,"black", dt] == 65535))), col = "red", lwd = 4)
            hist(pw.m[11:1985,11:1985,"grey", dt], breaks = "fd", add = T, col = "black")
        }
        hist(pw.m[,,"white", dt], breaks = "fd", col = "darkgreen", border = "green3", ylim = c(0,250), xlab = "", ylab = "", main = ""); {
            rect(0, 0, lim$w$dv, 600, col = adjustcolor("skyblue", alpha = 0.5), border = NA)
            rect(lim$w$dv, 0, lim$w$dm, 600, col = adjustcolor("cyan3", alpha = 0.5), border = NA)
            rect(lim$w$bm, 0, lim$w$bv, 600, col = adjustcolor("gold", alpha = 0.5), border = NA)
            rect(lim$w$bv, 0, 65535, 600, col = adjustcolor("orange", alpha = 0.5), border = NA)
            rect(qJohnson(0.001, JohnsonFit(pw.m[,,"black", dt])), 0, 
                 qJohnson(0.999, JohnsonFit(pw.m[,,"black", dt])), 600, col = adjustcolor("slateblue1", alpha = 0.5), border = NA)
            lines(c(0,0), c(0,length(which(pw.m[,,"white", dt] == 0))), col = "blue", lwd = 4)
            lines(c(65535,65535), c(0,length(which(pw.m[,,"black", dt] == 65535))), col = "red", lwd = 4)
            hist(pw.m[11:1985,11:1985,"white", dt], breaks = "fd", add = T, col = "black")
        }
        dev.off()
    }
    pdf(paste0(fpath, "threshold-histograms-legend1.pdf")); {
        plot.new()
        legend("center", pch = 15, cex = 0.8, col = adjustcolor(c("blue", "skyblue", "cyan3", "gold", "orange"), alpha = 0.4),
               legend = c("No response", "Very dim", "Dim", "Bright", "Very Bright"), horiz = T, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "threshold-histograms-legend1.pdf"))
    }
    pdf(paste0(fpath, "threshold-histograms-legend2.pdf")); {
    plot.new()
    legend("center", cex = 0.8, lty = 1, lwd = 2, col = c("blue", "red", "green3", "black"),
           legend = c("Dead", "Hot", "Hist incl. edge px", "Hist excl. edge px"), horiz = T, bty = "n")
    dev.off()
    crop.pdf(paste0(fpath, "threshold-histograms-legend2.pdf"))
    }
}

####################################################################################################

# IDENTIFY BRIGHT/DIM COLUMNS & ROWS                                                            ####

bp <- rbind(bp,
            data.frame(which(find.lines(pw.m[,,"black", dt]) > 0, arr.ind = T), type = "line.b"))
    
    newlines <- find.lines(pw.m[,,"black", dt], dim.lines = T)      # strip of 35px in bottom corner
    newlines <- find.lines(pw.m[,,"black", dt], horizontal = T)     # nothing
    newlines <- find.lines(pw.m[,,"black", dt], horizontal = T, dim.lines = T)
    table(newlines)

    bp <- bp[order(bp$type),]
    bp <- bp[!duplicated(bp[,1:2]),]
    table(bp$type)
    
    
}

plot(bp[,1:2], pch = 15, col = cat.cols[bp$type], asp = T)


####################################################################################################

# FIND LOCALLY EXTREME VALUES                                                                   ####

# load median-filtered differences
{
    md.b <- readRDS("./Other-data/Median-diffs-black.rds")
    md.g <- readRDS("./Other-data/Median-diffs-grey.rds")
    md.w <- readRDS("./Other-data/Median-diffs-white.rds")
}

s.hist(md.b[[dt]]); abline(v = mad(md.b[[dt]], na.rm = T) * c(-3, 3), col = "red")
mad(md.b[[dt]], na.rm = T); mad(md.g[[dt]], na.rm = T); mad(md.w[[dt]], na.rm = T)


image(1:1996, 1:1996, threshold(md.b[[dt]], level = mad(md.b[[dt]], na.rm = T) * 3), col = c(NA, "red"))
image(1:1996, 1:1996, threshold(md.g[[dt]], level = mad(md.g[[dt]], na.rm = T) * 3), col = c(NA, "blue"))
image(1:1996, 1:1996, threshold(md.w[[dt]], level = mad(md.w[[dt]], na.rm = T) * 3), col = c(NA, "green3"))

bp <- rbind(bp,
            data.frame(which(threshold(md.b[[dt]], level = mad(md.b[[dt]], na.rm = T) * 3) > 0, arr.ind = T),
                       type = "l.bright"),
            data.frame(which(threshold(md.b[[dt]], level = mad(md.b[[dt]], na.rm = T) *-3) == 0, arr.ind = T),
                       type = "l.dim"))

pixel.image(md.b[[dt]])

bp <- bp[order(bp$type),]
bp <- bp[!duplicated(bp[,1:2]),]

table(bp$type)

# histogram with & without locally bright/dim px, for comparison
tmp <- pw.m[,,"black", dt]
tmp[as.matrix(bp[,1:2])] <- NA

hist(pw.m[,,"black", dt], col = "black", breaks = "fd", ylim = c(0,500), xlim = c(0,15000))
hist(tmp, add = T, breaks = "fd", col = adjustcolor("cyan3", alpha = 0.4), border = adjustcolor("cyan3", alpha = 0.1))

hist(pw.m[,,"black", dt][as.matrix(which(threshold(md.b[[dt]], level = mad(md.b[[dt]], na.rm = T) * 7) > 0, arr.ind = T))],
     breaks = "fd", border = "gold", col = "gold", add = T)

qq <- unlist(lapply(c(3:15), function(x) length(which(threshold(md.g[[dt]], level = mad(md.g[[dt]], na.rm = T) * x) > 0))))
plot(qq)

qJohnson(0.999, JohnsonFit(md.g[[dt]][!is.na(md.g[[dt]])]))

####################################################################################################

# FIND SINGLETON BRIGHT/DIM PIXELS USING SINGLE-POINT KERNEL                                    ####

####################################################################################################

# CLASSIFY BAD PIXELS & SUPERCLUSTERS                                                           ####

# actually, should probably do this before the 
####################################################################################################

# MERGE INTO STATE SPACE                                                                        ####

####################################################################################################

# TRACK CHANGES IN STATE SPACE                                                                  ####

####################################################################################################
