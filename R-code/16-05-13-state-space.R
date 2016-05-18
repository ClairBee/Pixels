
# START FROM RAW IMAGES, CLASSIFY EVERYTHING FROM SCRATCH

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Med-diff-classification/fig/"


####################################################################################################
# GET ALL BAD PIXEL MAPS (WITH LINES)
# TRANSFER MAPS TO IMAGE ARRAYS WITH NUMBERED CATEGORIES
# MARK SUPERCLUSTERS
# CREATE STATE SPACE/TRANSITION DIAGRAMS
# OVERPLOT 'LOCALLY BRIGHT' POINTS (BY MED DIFF OR BY POINT CONVOLUTION)
# LOOK AT STABILITY OF 'LOCALLY BRIGHT' POINTS & IDENTIFY USEFUL THRESHOLD
####################################################################################################

# revised process to omit parametric model: instead, use difference from local median

####################################################################################################

# then repeat the process over the old images
# clean up functions: should be able to apply to old panel structure too!

# should also be able to specify degree of edge cropping

####################################################################################################

# SETUP                                                                                         ####

bp <- readRDS(paste0(fpath, "bad-px-maps.rds"))
load.pixel.means()

# load median-filtered differences
{
    md.b <- readRDS("./Other-data/Median-diffs-black.rds")
    md.g <- readRDS("./Other-data/Median-diffs-grey.rds")
    md.w <- readRDS("./Other-data/Median-diffs-white.rds")
}

cat <- c("no response", "dead", "hot", 
         "v.bright", "bright", "line.b", "l.bright", "s.bright", 
         "screen spot", "line.d", "edge",
         "v.dim", "dim", "l.dim", "s.dim")

cat.cols <- c("purple", "black", "magenta3", 
              "red", "orange", "gold", "gold", NA,
              "grey", "violet", NA, 
              "green3", "green", "lightskyblue", NA)

headerCat <- c("normal", gsub("[ ]", "", gsub("[.]", "", cat)))
fancyCat <- c("Normal", "No response", "Dead", "Hot", "V. bright", "Bright", "Bright line", "Locally bright", "Slighly bright", "Screen spot", "Dim line", "Edge", "V. dim", "Dim", "Locally dim", "Slightly dim")

####################################################################################################

# FIND GLOBALLY EXTREME VALUES                                                                  ####

bp <- list()

# hot, v.bright, v.dim, dead, no response, edge, dim spot
for (dt in dimnames(pw.m)[[4]]) {
    bp[[dt]] <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
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
    
    bp[[dt]] <- bp[[dt]][order(bp[[dt]]$type),]
    bp[[dt]] <- bp[[dt]][!duplicated(bp[[dt]][,1:2]),]
    bp[[dt]] <- bp[[dt]][!bp[[dt]]$type %in% c("s.bright", "s.dim"),]     # retain only more severe categories
}


table(bp$"160430"$type)

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

for (dt in dimnames(pw.m)[[4]]) {
    bp[[dt]] <- rbind(bp[[dt]],
            data.frame(which(find.lines(pw.m[,,"black", dt]) > 0, arr.ind = T), type = "line.b"))
    
    bp[[dt]] <- bp[[dt]][order(bp[[dt]]$type),]
    bp[[dt]] <- bp[[dt]][!duplicated(bp[[dt]][,1:2]),]
    bp[[dt]] <- bp[[dt]][!bp[[dt]]$type %in% c("s.bright", "s.dim"),]     # retain only more severe categories
}

# check for other line types (none were found, so not included in main algorithm for sake of time)
{   
    newlines <- find.lines(pw.m[,,"black", dt], dim.lines = T)      # strip of 35px in bottom corner
    newlines <- find.lines(pw.m[,,"black", dt], horizontal = T)     # nothing
    newlines <- find.lines(pw.m[,,"black", dt], horizontal = T, dim.lines = T)
    table(newlines)
}

plot(bp[,1:2], pch = 15, col = cat.cols[bp$type], asp = T)

####################################################################################################

# STATE SPACE (EXCLUDING SUPERCLUSTERS)                                                         ####

# turn bad pixel map into list of image files
bp.im <- lapply(bp, bpx2im)

bp.chng <- list()


for (i in 1:11) {
    tmp <- table(bp.im[[i]], bp.im[[i+1]])
    colnames(tmp) <- fancyCat[as.integer(colnames(tmp)) + 1]
    rownames(tmp) <- fancyCat[as.integer(rownames(tmp)) + 1]
    
    tmp[tmp == 0] <- "-"
    write.csv(tmp, paste0(fpath, "bpx-change-", names(bp)[i], "-", names(bp)[i+1], ".csv"), quote = F)
}

for (i in 1:11) {
    tmp <- table(bp.im[[i]], bp.im[[i+1]])
    tmp[tmp == 0] <- NA
    tmp <- round((tmp / rowSums(tmp, na.rm = T)) * 100, 1)
    colnames(tmp) <- fancyCat[as.integer(colnames(tmp)) + 1]
    rownames(tmp) <- fancyCat[as.integer(rownames(tmp)) + 1]
    tmp[is.na(tmp)] <- "-"
    
    write.csv(tmp, paste0(fpath, "bpx-change-", names(bp)[i], "-", names(bp)[i+1], "-prop.csv"), quote = F)
}

bp.change <- array(dim = c(16, 16, 11),
                   dimnames = list(c("normal", cat),
                                   c("normal", cat),
                                   apply(cbind(names(bp)[1:11], names(bp)[2:12]), 1, paste, collapse = ".")))

for (i in 1:11) {
    tmp <- table(bp.im[[i]], bp.im[[i+1]])
    bp.change[as.integer(rownames(tmp)) + 1, as.integer(colnames(tmp)) + 1, i] <- tmp
}

saveRDS(bp.change, paste0(fpath, "bp-change.rds"))

zz <- array(apply(bp.change, 3, function(x) 100 * x / rowSums(x, na.rm = T)), 
            dim = dim(bp.change), dimnames = dimnames(bp.change))

mean.change <- apply(zz, 1:2, mean, na.rm = T)
mean.change[mean.change == 0] <- NA
mean.change <- round(mean.change[-c(8,9, 11, 15, 16), -c(8,9, 11, 15, 16)], 2)
mean.change[is.na(mean.change)] <- "-"
write.csv(mean.change, paste0(fpath, "mean-change.csv"), quote = F)

mean.px <- apply(bp.change, 1:2, mean, na.rm = T)
mean.px[mean.px == 0] <- NA
mean.px <- round(mean.px[-c(8,9, 11, 15, 16), -c(8,9, 11, 15, 16)], 0)
mean.px[is.na(mean.px)] <- "-"
write.csv(mean.px, paste0(fpath, "mean-px.csv"), quote = F)

####################################################################################################

# STATE SPACE (INCLUDING SUPERCLUSTERS)                                                         ####

# get superclusters from bad pixel map
{
    bp <- lapply(bp, superclusters)
    saveRDS(bp, paste0(fpath, "bad-px-maps.rds"))
}

sc <- lapply(lapply(lapply(bp, "[", , c("x", "y", "sc.type")), setNames, c("x", "y", "type")), bpx2im)

sc.change <- array(dim = c(5, 5, 11),
                   dimnames = list(c("normal", levels(bp[[1]]$sc.type)),
                                   c("normal", levels(bp[[1]]$sc.type)),
                                   apply(cbind(names(bp)[1:11], names(bp)[2:12]), 1, paste, collapse = ".")))

for (i in 1:11) {
    sc.change[ , , i] <- table(sc[[i]], sc[[i+1]])
}

sc.rate <- array(apply(sc.change, 3, function(x) 100 * x / rowSums(x)), dim = dim(sc.change), dimnames = dimnames(sc.change))

write.csv(round(apply(sc.change, 1:2, mean), 0), paste0(fpath, "sc-mean-change.csv"), quote = F)
write.csv(round(apply(sc.rate, 1:2, mean), 1), paste0(fpath, "sc-mean-rate.csv"), quote = F)

# classify states including superclusters

####################################################################################################

# LOOK AT CORRELATION BETWEEN NOISE PATTERNS IN ALL IMAGES                                      ####

# also test for complete spatial randomness

# FIND LOCALLY EXTREME VALUES                                                                   ####

s.hist(md.b[["160430"]], ylim = c(0,1000))
abline(v = mad(pw.m[,,"black", "160430"], na.rm = T) * c(-2, 2), col = "red")
mad(md.b[[dt]], na.rm = T); mad(md.g[[dt]], na.rm = T); mad(md.w[[dt]], na.rm = T)


image(1:1996, 1:1996, threshold(md.b[["160430"]], level = mad(pw.m[,,"black", "160430"]) * 2), col = c(NA, "red"))
image(1:1996, 1:1996, threshold(md.g[[dt]], level = mad(md.g[[dt]], na.rm = T) * 3), col = c(NA, "blue"))
image(1:1996, 1:1996, threshold(md.w[[dt]], level = mad(md.w[[dt]], na.rm = T) * 3), col = c(NA, "green3"))

bpx <- list()
for (dt in names(bp)) {
    bpx[[dt]] <- rbind(bp[[dt]][,1:3],
                       setNames(data.frame(which(threshold(md.b[[dt]], level = mad(pw.m[,,"black", dt]) * 2) > 0, arr.ind = T),
                                           type = "l.bright"), c("x", "y", "type")),
                       setNames(data.frame(which(threshold(md.b[[dt]], level = mad(pw.m[,,"black", dt]) * -2) == 0, arr.ind = T),
                                           type = "l.dim"), c("x", "y", "type")))
    bpx[[dt]] <- bpx[[dt]][order(bpx[[dt]]$type),]
    bpx[[dt]] <- bpx[[dt]][!duplicated(bpx[[dt]][,1:2]),]
}

table(bpx[["160314"]]$type)

bp.im <- lapply(bpx, bpx2im)

bp.change <- array(dim = c(16, 16, 11),
                   dimnames = list(c("normal", cat),
                                   c("normal", cat),
                                   apply(cbind(names(bp)[1:11], names(bp)[2:12]), 1, paste, collapse = ".")))

for (i in 1:11) {
    tmp <- table(bp.im[[i]], bp.im[[i+1]])
    bp.change[as.integer(rownames(tmp)) + 1, as.integer(colnames(tmp)) + 1, i] <- tmp
}

zz <- array(apply(bp.change, 3, function(x) 100 * x / rowSums(x, na.rm = T)), 
            dim = dim(bp.change), dimnames = dimnames(bp.change))

mean.change <- apply(zz, 1:2, mean, na.rm = T)
mean.change[mean.change == 0] <- NA
mean.change <- round(mean.change[-c(9, 11, 16), -c(9, 11, 16)], 2)
mean.change[is.na(mean.change)] <- "-"
write.csv(mean.change, paste0(fpath, "mean-change-local-2-mad.csv"), quote = F)

mean.px <- apply(bp.change, 1:2, mean, na.rm = T)
mean.px[mean.px == 0] <- NA
mean.px <- round(mean.px[-c(9, 11, 16), -c(9, 11, 16)], 0)
mean.px[is.na(mean.px)] <- "-"
write.csv(mean.px, paste0(fpath, "mean-px-local-2-mad.csv"), quote = F)

# remove bad pixels from image and check remaining image.

im <- pw.m[,,"white", "160430"]
im[as.matrix(bpx$"160430"[,1:2])] <- NA

s.hist(im)

pdf(paste0(fpath, "tmp.pdf"))
pixel.image(im, breaks = sd.levels())
dev.off()


#---------------------------------------------------------------------------------------------
# scratch
{
    
# check stability of categories

md.b.th.600 <- array(unlist(lapply(md.b, threshold, level = 600)), dim = c(1996, 1996, 12))
md.b.th.sum <- apply(md.b.th, 1:2, sum)

md.b.th.600.str <- apply(md.b.th.600, 1:2, paste, collapse = "")

table(md.b.th.600.str)


unlist(lapply(md.b.th, sum, na.rm = T))

{
    pdf(paste0(fpath, "tmp.pdf"))
    image(1:1996, 1:1996, md.b.th$"160430", col = c(NA, "gold"))
    image(1:1996, 1:1996, md.b.th$"160314", col = c(NA, "red"), add = T)
    dev.off()
}



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
}

####################################################################################################

# FIND SINGLETON BRIGHT/DIM PIXELS USING SINGLE-POINT KERNEL                                    ####
# check stability of categories
# better or worse than median differencing?

k3 <- matrix(c(rep(-1, 4), 8, rep(-1,4)), ncol = 3)
k5 <- matrix(c(rep(-1, 12), 24, rep(-1,12)), ncol = 5)

# 3-square kernel
{
    hh.3 <- r2m(focal(m2r(pw.m[,,"black", "160430"]), k3))
    
    r2m(focal(m2r(matrix(c(rep(0, 40), ceiling(2 * mad(pw.m[,,"black", "160430"])),
                           rep(0, 40)), ncol = 9)), k3))
    
    bp.k3 <- rbind(setNames(data.frame(which(threshold(hh.3, 
                                                       level = floor(mad(pw.m[,,"black", "160430"]) * 2)) > 0,
                                             arr.ind = T),
                                       type = "l.bright"), c("x", "y", "type")),
                   setNames(data.frame(which(threshold(hh.3, 
                                                       level = floor(mad(pw.m[,,"black", "160430"]) * -2)) == 0,
                                             arr.ind = T),
                                       type = "l.dim"), c("x", "y", "type")))
    
    # repeat for grey & white images
    {
        bp.k3.g <- rbind(setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"grey", "160430"]), k3)), 
                                                             level = floor(mad(pw.m[,,"grey", "160430"]) * 2)) > 0,
                                                   arr.ind = T),
                                             type = "l.bright"), c("x", "y", "type")),
                         setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"grey", "160430"]), k3)), 
                                                             level = floor(mad(pw.m[,,"grey", "160430"]) * -2)) == 0,
                                                   arr.ind = T),
                                             type = "l.dim"), c("x", "y", "type")))
        
        bp.k3.w <- rbind(setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"white", "160430"]), k3)), 
                                                             level = floor(mad(pw.m[,,"white", "160430"]) * 2)) > 0,
                                                   arr.ind = T),
                                             type = "l.bright"), c("x", "y", "type")),
                         setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"white", "160430"]), k3)), 
                                                             level = floor(mad(pw.m[,,"white", "160430"]) * -2)) == 0,
                                                   arr.ind = T),
                                             type = "l.dim"), c("x", "y", "type")))
    }
    
    zz <- data.frame(zz,
                     l.b.3 = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k3)) + 1], levels = c("normal", "l.bright", "l.dim")),
                     l.g.3 = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k3.g)) + 1], levels = c("normal", "l.bright", "l.dim")),
                     l.w.3 = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k3.w)) + 1], levels = c("normal", "l.bright", "l.dim")))
    
    table("med.diff" = zz$lw, "3-kernel" = zz$l.w.3)
}


# 5-square kernel
{
    bp.k5.b <- rbind(setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"black", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"black", "160430"]) * 24)) > 0,
                                               arr.ind = T),
                                         type = "l.bright"), c("x", "y", "type")),
                     setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"black", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"black", "160430"]) * -24)) == 0,
                                               arr.ind = T),
                                         type = "l.dim"), c("x", "y", "type")))
    
    bp.k5.g <- rbind(setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"grey", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"grey", "160430"]) * 24)) > 0,
                                               arr.ind = T),
                                         type = "l.bright"), c("x", "y", "type")),
                     setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"grey", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"grey", "160430"]) * -24)) == 0,
                                               arr.ind = T),
                                         type = "l.dim"), c("x", "y", "type")))
    
    bp.k5.w <- rbind(setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"white", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"white", "160430"]) * 24)) > 0,
                                               arr.ind = T),
                                         type = "l.bright"), c("x", "y", "type")),
                     setNames(data.frame(which(threshold(r2m(focal(m2r(pw.m[,,"white", "160430"]), k5)), 
                                                         level = floor(mad(pw.m[,,"white", "160430"]) * -24)) == 0,
                                               arr.ind = T),
                                         type = "l.dim"), c("x", "y", "type")))
    
    zz$l.b.5 <- ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k5.b)) + 1], levels = c("normal", "l.bright", "l.dim"))
    zz$l.g.5 <- ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k5.g)) + 1], levels = c("normal", "l.bright", "l.dim"))
    zz$l.w.5 <- ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.k5.w)) + 1], levels = c("normal", "l.bright", "l.dim"))
    
    
}

    cbind(table("med.diff" = zz$local, "3-kernel" = zz$l.w.3),
          table("med.diff" = zz$local, "5-kernel" = zz$l.w.5))

# tables of results: median differencing vs kernel
# all thresholded using 2 * image MAD as cutoff (scaled where necessary)
{
    #    BLACK IMAGES
    #                               3-kernel                           5-kernel
    #    med.diff       normal l.bright   l.dim            normal l.bright   l.dim
    #    normal        2066822   863514 1036222           3909977    34219   22362
    #    l.bright           11    16906     250               177    16831     159
    #    l.dim               0        0     291                 5        0     286
    
    #    GREY IMAGES
    #                   normal l.bright  l.dim             normal l.bright l.dim
    #    normal        2153540   838861 974157            3925408    25469 15681
    #    l.bright           25    16916    226                365    16673   129
    #    l.dim               0        0    291                 24        3   264
    
    #    WHITE IMAGES
    #                   normal l.bright  l.dim             normal l.bright l.dim
    #    normal        3619041   197087 150430            3964889     1201   468
    #    l.bright         1080    15948    139              14092     3073     2
    #    l.dim              29        2    260                265       11    15
}
# far larger number of pixels identified as bad using essentially same thresholds.
# Stick with median-differencing as a more direct & intuitive approach.

####################################################################################################

# PIXELWISE SD OF EACH CLASS                                                                    ####

# are noisy points already picked up by these classifications?

####################################################################################################

# FIND BRIGHT/DIM LINES IN OLD DATA                                                             ####

{
    old.b <- list("130613" = readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapBlack.tif", as.is = T),
                  "131002" = readTIFF("~/Documents/Pixels/Other-data/Old-data/131002-bad-pixel-map/BadPixelMapBlack.tif", as.is = T),
                  "140128" = readTIFF("~/Documents/Pixels/Other-data/Old-data/140128-bad-pixel-map/BadPixelMapBlack.tif", as.is = T),
                  "150126" = readTIFF("~/Documents/Pixels/Other-data/Old-data/150126-bad-pixel-map/BadPixelMapBlack.tif", as.is = T))
    old.b <- lapply(old.b, function(x) t(x[nrow(x):1, , drop = F]))
}

old.bc <- lapply(old.b, find.lines)
old.dc <- lapply(old.b, find.lines, dim.lines = T)

table(old.dc$"130613")      # column 988 (lower)
table(old.dc$"131002")
table(old.dc$"140128")      # column 745 (lower)
table(old.dc$"150126")      # column 809 (lower)


# transect plots
{
    # bright column 988 (lower) in 13-06-13
    pdf("./Other-data/Old-data/plots/bright-130613.pdf", height = 4, width = 7); {
        par(mar = c(2,2,3,1))
        o.plot(old.b$"130613"[988, ], xlim = c(0,992), ylim = c(5000, 6000), main = "bright column 988 (lower) in 13-06-13")
        o.plot(old.b$"130613"[987, ], add = T, col = adjustcolor("cyan3", alpha = 0.4))
        o.plot(old.b$"130613"[989, ], add = T, col = adjustcolor("orange", alpha = 0.4))
        dev.off()
    }

# bright column 745 (lower) in 14-01-28
    pdf("./Other-data/Old-data/plots/bright-140128.pdf", height = 4, width = 7); {
        par(mar = c(2,2,3,1))
        o.plot(old.b$"140128"[745, ], xlim = c(1, 992), main = "bright column 745 (lower) in 14-01-28")
        o.plot(old.b$"140128"[744, ], add = T, col = adjustcolor("cyan3", alpha = 0.4))
        o.plot(old.b$"140128"[746, ], add = T, col = adjustcolor("orange", alpha = 0.4))
        dev.off()
    }

# bright column 809 (lower) in 15-01-26
    pdf("./Other-data/Old-data/plots/bright-150126.pdf", height = 4, width = 7); {
        par(mar = c(2,2,3,1))
        
    o.plot(old.b$"150126"[809, ], xlim = c(0, 992), ylim = c(5000,6500), main = "bright column 809 (lower) in 15-01-26")
    o.plot(old.b$"150126"[808, ], add = T, col = adjustcolor("cyan3", alpha = 0.4))
    o.plot(old.b$"150126"[810, ], add = T, col = adjustcolor("orange", alpha = 0.4))
    dev.off()
}

zz <- which(old.dc$"130613" > 0, arr.ind = T)

# dim column 208 (upper) in 13-06-13
    pdf("./Other-data/Old-data/plots/dim-130613-1.pdf", height = 4, width = 7); {
        par(mar = c(2,2,3,1))
    o.plot(old.b$"130613"[208, ], xlim = c(993, 1996), ylim = c(4000, 5600), main = "dim column 208 (upper) in 13-06-13")
    o.plot(old.b$"130613"[207, ], add = T, col = adjustcolor("cyan3", alpha = 0.4))
    o.plot(old.b$"130613"[206, ], add = T, col = adjustcolor("orange", alpha = 0.4))
    dev.off()
}

# dim column 1372 (upper) in 13-06-13
    pdf("./Other-data/Old-data/plots/dim-130613-2.pdf", height = 4, width = 7); {
        par(mar = c(2,2,3,1))
    o.plot(old.b$"130613"[1372, ], xlim = c(993, 1996), ylim = c(4000, 5600), main = "dim column 1372 (upper) in 13-06-13")
    o.plot(old.b$"130613"[1371, ], add = T, col = adjustcolor("cyan3", alpha = 0.4))
    o.plot(old.b$"130613"[1373, ], add = T, col = adjustcolor("orange", alpha = 0.4))
    dev.off()
    }
}

o.plot(old.b$"130613"[1372, ],add = T, col = "red")

####################################################################################################

# DOES MEDIAN-FILTERED DIFF PICK UP ALL 'GLOBALLY' BAD PIXELS?                                  ####

bp.local <- rbind(setNames(data.frame(which(threshold(md.b[["160430"]], level = mad(pw.m[,,"black", "160430"]) * 2) > 0, arr.ind = T),
                                       type = "l.bright"), c("x", "y", "type")),
                   setNames(data.frame(which(threshold(md.b[["160430"]], level = mad(pw.m[,,"black", "160430"]) * -2) == 0, arr.ind = T),
                                       type = "l.dim"), c("x", "y", "type")))

bp.global <- bp$"160430"

hh <- data.frame(global = ordered(c("normal", cat)[c(bpx2im(bp.global)) + 1], levels = c("normal", cat)),
                 local = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.local)) + 1], levels = c("normal", "l.bright", "l.dim")))

# repeat for grey & white images, to check composition
{
    bp.l.g <- rbind(setNames(data.frame(which(threshold(md.g[["160430"]], level = mad(pw.m[,,"grey", "160430"]) * 2) > 0, arr.ind = T),
                                        type = "l.bright"), c("x", "y", "type")),
                    setNames(data.frame(which(threshold(md.g[["160430"]], level = mad(pw.m[,,"grey", "160430"]) * -2) == 0, arr.ind = T),
                                        type = "l.dim"), c("x", "y", "type")))
    
    bp.l.w <- rbind(setNames(data.frame(which(threshold(md.w[["160430"]], level = mad(pw.m[,,"white", "160430"]) * 2) > 0, arr.ind = T),
                                        type = "l.bright"), c("x", "y", "type")),
                    setNames(data.frame(which(threshold(md.w[["160430"]], level = mad(pw.m[,,"white", "160430"]) * -2) == 0, arr.ind = T),
                                        type = "l.dim"), c("x", "y", "type")))
    
    zz <- data.frame(zz,
                     lg = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.l.g)) + 1], levels = c("normal", "l.bright", "l.dim")),
                     lw = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(bp.l.w)) + 1], levels = c("normal", "l.bright", "l.dim")))
}

qq <- cbind(table(zz$global, zz$local),
            table(zz$global, zz$lg),
            table(zz$global, zz$lw))

rownames(qq) <- fancyCat
qq <- qq[rowSums(qq) > 0,]
qq[qq == 0] <- "-"

write.csv(qq, paste0(fpath, "global-vs-local-cat.csv"), quote = F)

# compare px identified in black/white/grey images
{
    tbl <- cbind(table(zz[zz$local == "normal", c("lg", "lw")]),
                 table(zz[zz$local == "l.bright", c("lg", "lw")]),
                 table(zz[zz$local == "l.dim", c("lg", "lw")]))
    
    rownames(tbl) <- c("Uniform", "Locally bright", "Locally dim")
    tbl[tbl == 0] <- "-"
    
    write.csv(tbl, paste0(fpath, "local-classes.csv"), quote = F)
}

####################################################################################################

# CHECK LOCALLY EXTREME PX IN ALL IMAGES, ALL POWERS                                            ####

# does including the white image actually add anything?

md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")
md.w <- readRDS("./Other-data/Median-diffs-white.rds")

{
    zz <- lapply(dimnames(pw.m)[[4]], 
                 function(x) data.frame("b" = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(rbind(setNames(data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T),
                                                                                                                          type = "l.bright"), c("x", "y", "type")),
                                                                                                      setNames(data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * -2) == 0, arr.ind = T),
                                                                                                                          type = "l.dim"), c("x", "y", "type"))))) + 1],
                                                      levels = c("normal", "l.bright", "l.dim")),
                                        "g" = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(rbind(setNames(data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * 2) > 0, arr.ind = T),
                                                                                                                          type = "l.bright"), c("x", "y", "type")),
                                                                                                      setNames(data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * -2) == 0, arr.ind = T),
                                                                                                                          type = "l.dim"), c("x", "y", "type"))))) + 1],
                                                      levels = c("normal", "l.bright", "l.dim")),
                                        "w" = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(rbind(setNames(data.frame(which(threshold(md.w[[x]], level = mad(pw.m[,,"white", x]) * 2) > 0, arr.ind = T),
                                                                                                                          type = "l.bright"), c("x", "y", "type")),
                                                                                                      setNames(data.frame(which(threshold(md.w[[x]], level = mad(pw.m[,,"white", x]) * -2) == 0, arr.ind = T),
                                                                                                                          type = "l.dim"), c("x", "y", "type"))))) + 1],
                                                      levels = c("normal", "l.bright", "l.dim"))))
}

        
names(zz) <- dimnames(pw.m)[[4]]
                               
rbind.fill(lapply(lapply(lapply(lapply(zz, function(x) table(x[x$b == "normal" & x$g == "normal", c("w")])), as.matrix), t), as.data.frame))
