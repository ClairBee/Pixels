
# DETECTION OF SLIGHTLY BRIGHT OR NON-RESPONSIVE LINES
# (creating plots for ./Notes/Line-detection/Line-detection.tex)

library("IO.Pixels"); library("CB.Misc")
if (getwd() != "/home/clair/Documents/Pixels") {setwd("/home/clair/Documents/Pixels")}
fpath <- "./Notes/Line-detection/fig/"

####################################################################################################

# SETUP & CONVOLUTION                                                                           ####

load.pixel.means()

# perform convolution with 5x5 kernel
k5 <- matrix(c(rep(-1, 10), rep(4, 5), rep(-1, 10)), ncol = 5)
conv.sq5 <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k5), r2m)


qq <- alply(pw.m[,,"black", ], 3, convolve.lines, k = 5, .dims = T)
            


# define regions of interest
{
    im1 <- list(xrng = c(359:550), yrng = c(1100:1300), col = 427, row = 1199)
    im2 <- list(xrng = c(766:896), yrng = c(50:200), col = 809)
    
    line.px <- rbind(cbind(x = 427, y = c(1201:1996)),
                     cbind(x = 809, y = c(1:177)))
    
    line1 <- cbind(x = 427, y = c(1201:1996))
    line2 <- cbind(x = 809, y = c(1:177))
    
    edge.px <- matrix(c(rep(panel.edges()$x[1:16], 992), rep(panel.edges()$x[2:17]-1, 1004), 
                        sort(rep(1:992, 16)), sort(rep(993:1996, 16))), ncol = 2)
}


####################################################################################################

# RUN OVER VARIOUS THRESHOLDS                                                                   ####

th5 <- list("3500" = threshold(conv.sq5$"160314", level = 3500),
            "4000" = threshold(conv.sq5$"160314", level = 4000),
            "4500" = threshold(conv.sq5$"160314", level = 4500),
            "5000" = threshold(conv.sq5$"160314", level = 5000),
            "5500" = threshold(conv.sq5$"160314", level = 5500),
            "6000" = threshold(conv.sq5$"160314", level = 6000),
            "6500" = threshold(conv.sq5$"160314", level = 6500))

# smooth all thresholded images with 11-kernel
bright.lines <- lapply(lapply(lapply(lapply(th5, m2r), focal, sm), r2m), threshold, level = 5.5)

# overplot each threshold
{
    line.cols <- c("gold", "green2", "cyan3", "dodgerblue", "purple", "magenta3", "red")
    
    pdf(paste0(fpath, "bright-lines-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, bright.lines$"3500", col = c(NA, line.cols[1]), xlim = range(im1$xrng), ylim = range(im1$yrng))
        for (i in 2:length(th5)) {
            image(1:1996, 1:1996, bright.lines[[i]], col = c(NA, line.cols[i]), add = T)
        }
    }; dev.off()
    pdf(paste0(fpath, "bright-lines-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, bright.lines$"3500", col = c(NA, line.cols[1]), xlim = range(im2$xrng), ylim = range(im2$yrng))
        for (i in 2:length(th5)) {
            image(1:1996, 1:1996, bright.lines[[i]], col = c(NA, line.cols[i]), add = T)
        }
    }; dev.off()
    pdf(paste0(fpath, "bright-lines-thresholds.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, bright.lines$"3500", col = c(NA, line.cols[1]))
        for (i in 2:length(th5)) {
            image(1:1996, 1:1996, bright.lines[[i]], col = c(NA, line.cols[i]), add = T)
        }
    }; dev.off()
    
    pdf(paste0(fpath, "bright-lines-legend.pdf")); {
        plot.new()
        legend("center", col = line.cols, pch = 16,
               legend = names(th5), title = "  Thresholded at  ")
    }; dev.off()
    crop.pdf(paste0(fpath, "bright-lines-legend.pdf"))
    
}

####################################################################################################

# CLASSIFICATION RATES                                                                          ####

# create mask for pos/neg expected result
l.mask <- array(0, dim = c(1996, 1996))
l.mask[line.px] <- 1

rates <- rbind.fill(lapply(lapply(lapply(lapply(bright.lines, table, l.mask), c), t), data.frame))
{
    colnames(rates) <- c("TN", "FP", "FN", "TP"); rownames(rates) <- names(th5)
    rates$true.pos.rate <- round((rates$TP / (rates$TP + rates$FN)) * 100, 1)
#    rates$true.neg.rate <- round((rates$TN / (rates$TN + rates$FP)) * 100, 1)
    rates$false.disc.rate <- round((rates$FP / (rates$TP + rates$FP)) * 100, 1)
#    rates$false.neg <- round((rates$FN / (rates$TP + rates$FN)) * 100, 1)
#    rates$precision <- round((rates$TP / (rates$TP + rates$FP)) * 100, 1)
#    rates$accuracy <- round(((rates$TP + rates$TN) / 1996^2) * 100, 1)
    rates$F1.score <- round(2 * rates$TP / (2 * rates$TP + rates$FP + rates$FN), 2)
}
write.csv(rates, paste0(fpath, "sq5-classification-rates.csv"), quote = F)

pdf(paste0(fpath, "tmp.pdf"))
image(1:1996, 1:1996, bright.lines$"5500", col = c(NA, "magenta3"))
dev.off()

# try cutting at 6 instead: effect?
{
    bright.lines2 <- lapply(lapply(lapply(lapply(th5, m2r), focal, sm), r2m), threshold, level = 6.5)
    rates2 <- rbind.fill(lapply(lapply(lapply(lapply(bright.lines2, table, l.mask), c), t), data.frame))
    colnames(rates2) <- c("TN", "FP", "FN", "TP"); rownames(rates2) <- names(th5)
    rates2$true.pos.rate <- round((rates2$TP / (rates2$TP + rates2$FN)) * 100, 1)
    rates2$false.disc.rate <- round((rates2$FP / (rates2$TP + rates2$FP)) * 100, 1)
    rates2$F1.score <- round(2 * rates2$TP / (2 * rates2$TP + rates2$FP + rates2$FN), 2)
}
####################################################################################################

# FINAL CLEANING & CLASSIFICATION                                                               ####

# for each column, check max & min pixel identified as part of a column,
# and proportion in between.

bl <- bright.lines$"5500"

cc <- clump(m2r(bl))


# get x-y coordinates of each point
xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                 id = getValues(cc)[!is.na(getValues(cc))])

mids <- ddply(xy, .(id), summarise, xm = mean(x), ym = mean(y), ymax = max(y), ymin = min(y))
    

# retain only columns where more than 1% of pixels were identified as bright
bl[as.integer(names(which(table(xy$x[xy$y < 993]) <= 10))), 1:992] <- NA        # clean lower panel
bl[as.integer(names(which(table(xy$x[xy$y >= 993]) <= 10))), 993:1996] <- NA    # clean upper panel

image(1:1996, 1:1996, bl, col = c(NA, "blue"))

cc.bl <- clump(m2r(bl))

xy.bl <- data.frame(xyFromCell(cc.bl, which(!is.na(getValues(cc.bl)))),
                 id = getValues(cc.bl)[!is.na(getValues(cc.bl))])

# get distance to nearest (horizontal) edge
# upper
df <- rbind(ddply(xy.bl[xy.bl$y > 992,], .(x), summarise, edge.dist = 1996 - min(y), 
                  freq = length(x), edge.prop = freq / edge.dist),
            ddply(xy.bl[xy.bl$y <= 992,], .(x), summarise, edge.dist = max(y), 
                  freq = length(x), edge.prop = freq / edge.dist))

df <- ddply(xy.bl, .(x), summarise, edge.dist = 1996 - min((y.min > 992) * 1996, u.dist = (y.min > 992) * y.min,
            l.dist = (y.min <= 992) * y.max)

table(xy$id, xy$x)

####################################################################################################

# DIM & HORIZONTAL LINES                                                                        ####

# MATCH SHORT SEGMENTS TO BAD PIXELS                                                            ####

mids <- ddply(xy, .(id), summarise, xm = mean(x), ym = mean(y), ymax = max(y), ymin = min(y))

bp <- readRDS("./Models/Simple-parametric/Bad-px-new-thresholds.rds")$"160314"

bpx <- bp[bp$type %in% c("v.bright", "bright", "s.bright", "hot"), 1:2]

image(1:1996, 1:1996, bright.lines$"3500", col = c(NA, "orange"), asp = T, xlab = "", ylab = "",
      xlim = c(0,100), ylim = c(0,100))
points(bpx, pch = 15)

####################################################################################################

# POINT KERNEL                                                                                  ####

 kp <- matrix(c(0,0,0,0,0,   0,-1,-1,-1,0,  0,)
              
####################################################################################################

# HORIZONTAL CONVOLUTION                                                                        ####

conv.h <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, t(k5)), r2m)
c.h.th <- threshold(conv.h$"160314", level = 5500)

image(1:1996, 1:1996, c.h.th, col = c(NA, "magenta3"))
pixel.image(conv.h$"160314", xlim = range(im1$xrng), ylim = range(im1$yrng))

conv.mixed <- conv.sq5$"160314" - conv.h$"160314"

pixel.image(conv.mixed, xlim = range(im1$xrng), ylim = range(im1$yrng))

o.plot(conv.mixed[427,1150:1300], ylim = c(-10000,10000))
abline(h = 5500, col = "red")

th.mixed <- threshold(conv.mixed, level = 5500)
image(1:1996, 1:1996, th.mixed, col = c(NA, "magenta3"),
      xlim = range(im1$xrng), ylim = range(im1$yrng))

th.smoothed <- r2m(focal(m2r(th.mixed), matrix(rep(1, 11), ncol = 1)))
image(1:1996, 1:1996, threshold(th.smoothed, level = 5.5), col = c(NA, "cyan3"),
      xlim = range(im1$xrng), ylim = range(im1$yrng))
image(1:1996, 1:1996, threshold(th.smoothed, level = 5.5), col = c(NA, "cyan3"),
      xlim = range(im2$xrng), ylim = range(im2$yrng))
pdf(paste0(fpath, "tmp.pdf"))
image(1:1996, 1:1996, threshold(th.smoothed, level = 5.5), col = c(NA, "cyan3"))
dev.off()

####################################################################################################

# WHAT IS CAUSING THE SHORT LINE SEGMENTS?                                                      ####

image(1:1996, 1:1996, bright.lines$"5500", col = c(NA, "orange"), asp = T, xlab = "", ylab = "",
      xlim = c(0,100), ylim = c(0,100))

chk.all <- c(3, 6, 8, 14, 20, 30, 34, 44, 56, 97)

pdf(paste0(fpath, "short-lines-after-thresholding.pdf")); {
    par(mar = c(2,2,1,1), mfrow = c(3, 4))
    for (chk.c in chk.all) {
        o.plot(pw.m[chk.c,1:100,"black", "160314"], ylim = c(5000, 20000))
        o.plot(bright.lines$"5500"[chk.c, 1:100] * max(pw.m[chk.c,1:100,"black", "160314"]), add = T, col = "red")
    }
    dev.off()
}

####################################################################################################

conv <- convolve.lines(pw.m[,,"black", "160314"], k.size = 5)
th <- threshold(conv, level = 5500)
sm <- r2m(focal(m2r(th), matrix(rep(1, 11), ncol = 1)))

c.col <- 6
o.plot(pw.m[c.col,1:100, "black", "160314"], ylim = c(-10000, 20000))
o.plot(conv[c.col,1:100], add = T, col = "red")
abline(h = 5500, col = "cyan3")
points((sm[c.col,1:100] > 5.5) * conv[c.col,1:100], pch = 20, col = "blue")


####################################################################################################

# POSSIBLE CAUSES
# may be bridging between 2 points
# surrounding a single point whose neighbours are also slightly bright


pts <- array(0, dim = c(15, 30))

test.pts <- matrix(c(rep(10, 9), c(12:20)), ncol = 2)

par(mfrow = c(3,3))
for (i in 1:nrow(test.pts)) {
    pts[8, test.pts[i,]] <- 1400
    th <- threshold(r2m(focal(m2r(threshold(r2m(focal(m2r(pts), k5)), level = 5500)), matrix(rep(1, 11), ncol = 1))), level = 5.5)
    image(1:dim(pts)[[1]], 1:dim(pts)[[2]], pts, asp = T, col = c(NA, "black"), xlab = "", ylab = "")
    image(1:dim(pts)[[1]], 1:dim(pts)[[2]], th, add = T, col = c(NA, adjustcolor("cyan3", alpha = 0.4)))
    title(paste0(test.pts[i,2] - test.pts[1,1] - 1, "px gap -> ", sum(th, na.rm = T), "px line"))
    pts <- array(0, dim = c(15, 30))
}
par(mfrow = c(1,1))


th <- threshold(r2m(focal(m2r(threshold(r2m(focal(m2r(pts), k5)), level = 5500)), matrix(rep(1, 11), ncol = 1))), level = 5.5)

image(1:dim(pts)[[1]], 1:dim(pts)[[2]], pts, asp = T, col = c(NA, "black"), xlab = "", ylab = "")
image(1:dim(pts)[[1]], 1:dim(pts)[[2]], th, add = T, col = c(NA, adjustcolor("cyan3", alpha = 0.4)))

####################################################################################################

# PLOT PROGRESSION (AFTER 5X5 CONVOLUTION) IN GREY/WHITE IMAGES                                 ####

sq5.grey <- lapply(lapply(apply(pw.m[,,"grey", ], 3, m2r), focal, k5), r2m)
sq5.white <- lapply(lapply(apply(pw.m[,,"white", ], 3, m2r), focal, k5), r2m)

pdf(paste0(fpath, "hist-5x5-grey-compare.pdf")); {
    par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
    for (i in 1:12) {
        hist(sq5.grey[[i]], breaks = "fd", xlim = c(-10000, 15000), ylim = c(0, 1500), col = "black", main = "", xlab = "", ylab = "")
        hist(sq5.grey[[i]][edge.px], breaks = "fd", add = T, border = "cyan3", col = "cyan3")
        hist(sq5.grey[[i]][line1], breaks = "fd", add = T, border = "red", col = "red")
        hist(sq5.grey[[i]][line2], breaks = "fd", add = T, border = "orange", col = "orange")
    }
    dev.off()
}

pdf(paste0(fpath, "hist-5x5-white-compare.pdf")); {
    par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
    for (i in 1:12) {
        hist(sq5.white[[i]], breaks = "fd", xlim = c(-20000, 20000), ylim = c(0, 1500), col = "black", main = "", xlab = "", ylab = "")
        hist(sq5.white[[i]][edge.px], breaks = "fd", add = T, border = "cyan3", col = "cyan3")
        hist(sq5.white[[i]][line1], breaks = "fd", add = T, border = "red", col = "red")
        hist(sq5.white[[i]][line2], breaks = "fd", add = T, border = "orange", col = "orange")
    }
    dev.off()
}

{
    pdf(paste0(fpath, "hist-line1.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(sq5.white[[i]][line1], breaks = "fd", xlim = c(-10000, 20000), col = "gold", border = "gold", main = "", xlab = "", ylab = "")
            hist(sq5.grey[[i]][line1], breaks = "fd", add = T,
                 col = adjustcolor("cyan3", alpha = 0.4), border = adjustcolor("cyan3", alpha = 0.4))
            hist(conv.sq5[[i]][line1], breaks = "fd", add = T, col = "black")
            abline(v = c(0, 5500), col = "red", lty = c(1, 2))
        }
        dev.off()
    }
    
    pdf(paste0(fpath, "hist-line2.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(sq5.white[[i]][line2], breaks = "fd", xlim = c(-10000, 20000), col = "gold", border = "gold", main = "", xlab = "", ylab = "")
            hist(sq5.grey[[i]][line2], breaks = "fd", add = T,
                 col = adjustcolor("cyan3", alpha = 0.4), border = adjustcolor("cyan3", alpha = 0.4))
            hist(conv.sq5[[i]][line2], breaks = "fd", add = T, col = "black")
            abline(v = c(0, 5500), col = "red", lty = c(1, 2))
        }
        dev.off()
    }
}