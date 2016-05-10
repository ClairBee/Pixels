
# DETECTION OF SLIGHTLY BRIGHT OR NON-RESPONSIVE LINES
# (creating plots for ./Notes/Line-detection/Line-detection.tex)

library("IO.Pixels"); library("CB.Misc")

####################################################################################################

# SETUP                                                                                         ####

if (getwd() != "/home/clair/Documents/Pixels") {setwd("/home/clair/Documents/Pixels")}
load.pixel.means()
fpath <- "./Notes/Line-detection/fig/"

th.cols <- c("white", "paleturquoise1",  "cyan3", "skyblue", "green3",
             "gold", "orange", "red", "magenta3", "blue", "black")

# identify regions to look at in detail
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

# specify kernels
{
    k3 <- matrix(c(rep(-1, 3), rep(2, 3), rep(-1, 3)), ncol = 3)
    k5 <- matrix(c(rep(-1, 10), rep(4, 5), rep(-1, 10)), ncol = 5)
}

# perform convolution
{
    conv.sq3 <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k3), r2m)
    conv.sq5 <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k5), r2m)
}

####################################################################################################

# DISTRIBUTION OF SUBPANEL EDGE PIXELS                                                          ####

hist.cols <- adjustcolor(c("orange", "gold", "chartreuse3", "green3", "cyan3", "dodgerblue2",
                           "blue", "slateblue1", "purple", "magenta", "red"), alpha = 0.3)

plot(c(1:12), rep(1, 12), pch = 16, col = hist.cols)
      
# histograms of all edge pixels after convolution
{
    pdf(paste0(fpath, "hist-3x3-edges.pdf"))
    par(mar = c(2,2,1,1))
    s.hist(conv.sq3[[1]][edge.px], border = "lemonchiffon", main = "", xlab = "", ylab = "")
    for (i in 2:12) {
        hist(conv.sq3[[i]][edge.px], add = T, breaks = "fd", col = hist.cols[i], border = hist.cols[i])
    }
    dev.off()
    
    pdf(paste0(fpath, "hist-5x5-edges.pdf"))
    par(mar = c(2,2,1,1))
    s.hist(conv.sq5[[1]][edge.px], border = "lemonchiffon", main = "", xlab = "", ylab = "")
    for (i in 12:2) {
        hist(conv.sq5[[i]][edge.px], add = T, breaks = "fd", col = hist.cols[i], border = hist.cols[i])
    }
    dev.off()
}

# histograms of edge px / image MAD after convolution
{
    pdf(paste0(fpath, "hist-3x3-edges-normalised.pdf"))
    par(mar = c(2,2,1,1))
        hist(conv.sq3[[1]][edge.px] / mad(conv.sq3[[1]], na.rm = T), breaks = "fd", xlim = c(-10, 15), border = "lemonchiffon", main = "", xlab = "", ylab = "")
        for (i in 2:12) {
            hist(conv.sq3[[i]][edge.px] / mad(conv.sq3[[i]], na.rm = T), add = T, breaks = "fd", col = hist.cols[i], border = hist.cols[i])
        }
    dev.off()
    
    pdf(paste0(fpath, "hist-5x5-edges-normalised.pdf"))
    par(mar = c(2,2,1,1))
        hist(conv.sq5[[1]][edge.px] / mad(conv.sq5[[1]], na.rm = T), breaks = "fd", xlim = c(-10, 15), border = "lemonchiffon", main = "", xlab = "", ylab = "")
        for (i in 2:12) {
            hist(conv.sq5[[i]][edge.px] / mad(conv.sq5[[i]], na.rm = T), add = T, breaks = "fd", col = hist.cols[i], border = hist.cols[i])
        }
    dev.off()
}

write.csv(round(rbind("3x3" = unlist(lapply(conv.sq3, mad, na.rm = T)),
                      "5x5" = unlist(lapply(conv.sq5, mad, na.rm = T))),0),
          paste0(fpath, "MAD-comparison.csv"), quote = F)


# compare panel edges vs bad lines in each image
{
    pdf(paste0(fpath, "hist-3x3-compare.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(conv.sq3[[i]][edge.px], breaks = "fd", xlim = c(-2000, 7000), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq3[[i]][line.px], breaks = "fd", col = "red", border = "red", add = T)
        }
        dev.off()
    }

    pdf(paste0(fpath, "hist-5x5-compare.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(conv.sq5[[i]][edge.px], breaks = "fd", xlim = c(-2000, 20000), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq5[[i]][line.px], breaks = "fd", col = "red", border = "red", add = T)
        }
        dev.off()
    }

    pdf(paste0(fpath, "hist-3x3-compare-normalised.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(conv.sq3[[i]][edge.px] / mad(conv.sq3[[i]], na.rm = T), breaks = "fd", xlim = c(-10, 30), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq3[[i]][line.px] / mad(conv.sq3[[i]], na.rm = T), breaks = "fd", col = "red", border = "red", add = T)
        }
        dev.off()
    }

    pdf(paste0(fpath, "hist-5x5-compare-normalised.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
        for (i in 1:12) {
            hist(conv.sq5[[i]][edge.px] / mad(conv.sq5[[i]], na.rm = T), breaks = "fd", xlim = c(-10, 30), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq5[[i]][line.px] / mad(conv.sq5[[i]], na.rm = T), breaks = "fd", col = "red", border = "red", add = T)
        }
        dev.off()
    }
}

# value enhancement by threshold (example data)
{
    m1 <- matrix(rep(c(rep(0, 3), 300, rep(0, 3)), 7), ncol = 7)

    r2m(focal(m2r(m1 + 100), k3))
    r2m(focal(m2r(m1 + 2000), k5))
    
    m2 <- matrix(rep(c(rep(200, 3), 300, rep(0, 3)), 7), ncol = 7)
    r2m(focal(m2r(m2), k3))
    r2m(focal(m2r(m2), k5))
    
        write.csv(m1, paste0(fpath, "M1.csv"), row.names = F, quote = F)
        write.csv(r2m(focal(m2r(m1), k3)), paste0(fpath, "M1-3x3.csv"), row.names = F, quote = F)
        write.csv(r2m(focal(m2r(m1), k5)), paste0(fpath, "M1-5x5.csv"), row.names = F, quote = F)
        
        write.csv(m2, paste0(fpath, "M2.csv"), row.names = F, quote = F)
        write.csv(r2m(focal(m2r(m2), k3)), paste0(fpath, "M2-3x3.csv"), row.names = F, quote = F)
        write.csv(r2m(focal(m2r(m2), k5)), paste0(fpath, "M2-5x5.csv"), row.names = F, quote = F)
}

# numerical summary of poss. thresholds
{
    summ.sq3 <- rbind("Mean" = unlist(lapply(conv.sq3, mean, na.rm = T)),
                      "Median" = unlist(lapply(conv.sq3, median, na.rm = T)),
                      "Edge mean" = unlist(lapply(lapply(conv.sq3, "[", edge.px), mean, na.rm = T)),
                      "Edge median" = unlist(lapply(lapply(conv.sq3, "[", edge.px), median, na.rm = T)),
                      "Edge Q.99" = unlist(lapply(lapply(lapply(lapply(conv.sq3, "[", edge.px), function(x) x[!is.na(x)]), 
                                                         JohnsonFit, moment = "quant"), function(x) qJohnson(0.99, x))),
                      "Line1 mean" = unlist(lapply(lapply(conv.sq3, "[", line1), mean, na.rm = T)),
                      "Line1 median" = unlist(lapply(lapply(conv.sq3, "[", line1), median, na.rm = T)),
                      "Line2 mean" = unlist(lapply(lapply(conv.sq3, "[", line2), mean, na.rm = T)),
                      "Line2 median" = unlist(lapply(lapply(conv.sq3, "[", line2), median, na.rm = T)),
                      "MAD" = unlist(lapply(conv.sq3, mad, na.rm = T)))
    summ.sq3 <- rbind(summ.sq3,
                      "Edge / MAD" = summ.sq3["Edge mean",] / summ.sq3["MAD",],
                      "Line1 / MAD" = summ.sq3["Line1 mean",] / summ.sq3["MAD",],
                      "Line2 / MAD" = summ.sq3["Line2 mean",] / summ.sq3["MAD",])
                
    write.csv(round(summ.sq3, 0), paste0(fpath, "3x3-summary.csv"), quote = F)
    
    summ.sq5 <- rbind("Mean" = unlist(lapply(conv.sq5, mean, na.rm = T)),
                      "Median" = unlist(lapply(conv.sq5, median, na.rm = T)),
                      "Edge mean" = unlist(lapply(lapply(conv.sq5, "[", edge.px), mean, na.rm = T)),
                      "Edge median" = unlist(lapply(lapply(conv.sq5, "[", edge.px), median, na.rm = T)),
                      "Edge Q.99" = unlist(lapply(lapply(lapply(lapply(conv.sq5, "[", edge.px), function(x) x[!is.na(x)]), 
                                                         JohnsonFit, moment = "quant"), function(x) qJohnson(0.99, x))),
                      "Line1 mean" = unlist(lapply(lapply(conv.sq5, "[", line1), mean, na.rm = T)),
                      "Line1 median" = unlist(lapply(lapply(conv.sq5, "[", line1), median, na.rm = T)),
                      "Line2 mean" = unlist(lapply(lapply(conv.sq5, "[", line2), mean, na.rm = T)),
                      "Line2 median" = unlist(lapply(lapply(conv.sq5, "[", line2), median, na.rm = T)),
                      "MAD" = unlist(lapply(conv.sq5, mad, na.rm = T)))
    summ.sq5 <- rbind(summ.sq5,
                      "Edge / MAD" = summ.sq5["Edge mean",] / summ.sq5["MAD",],
                      "Line1 / MAD" = summ.sq5["Line1 mean",] / summ.sq5["MAD",],
                      "Line2 / MAD" = summ.sq5["Line2 mean",] / summ.sq5["MAD",])
    write.csv(round(summ.sq5, 0), paste0(fpath, "5x5-summary.csv"), quote = F)
}

####################################################################################################

# THRESHOLD SETTING                                                                             ####

# quick histogram to show difference in lines
{
    pdf(paste0(fpath, "hist-3x3-compare.pdf")); {
        par(mfrow = c(4,3), mar = c(2,2,1,1))
        for (i in 1:12) {
            hist(conv.sq3[[i]], breaks = "fd", xlim = c(-3000, 7000), ylim = c(0, 1000), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq3[[i]][edge.px], breaks = "fd", add = T, border = "cyan3", col = "cyan3")
            hist(conv.sq3[[i]][line1], breaks = "fd", add = T, border = "red", col = "red")
            hist(conv.sq3[[i]][line2], breaks = "fd", add = T, border = "orange", col = "orange")
        }
        dev.off()
    }
    
    pdf(paste0(fpath, "hist-5x5-compare.pdf")); {
        par(mfrow = c(4,3), mar = c(2,2,1,1))
        for (i in 1:12) {
            hist(conv.sq5[[i]], breaks = "fd", xlim = c(-5000, 15000), ylim = c(0, 1000), col = "black", main = "", xlab = "", ylab = "")
            hist(conv.sq5[[i]][edge.px], breaks = "fd", add = T, border = "cyan3", col = "cyan3")
            hist(conv.sq5[[i]][line1], breaks = "fd", add = T, border = "red", col = "red")
            hist(conv.sq5[[i]][line2], breaks = "fd", add = T, border = "orange", col = "orange")
        }
        dev.off()
    }
}

# go back to running everything over 16-03-14 black images

# try a variety of thresholds
th3 <- list("1500" = threshold(conv.sq3$"160314", level = 1500),
            "1750" = threshold(conv.sq3$"160314", level = 1750),
            "2000" = threshold(conv.sq3$"160314", level = 2000),
            "2250" = threshold(conv.sq3$"160314", level = 2250),
            "2500" = threshold(conv.sq3$"160314", level = 2500))

th5 <- list("3500" = threshold(conv.sq5$"160314", level = 3500),
            "4000" = threshold(conv.sq5$"160314", level = 4000),
            "4500" = threshold(conv.sq5$"160314", level = 4500),
            "5000" = threshold(conv.sq5$"160314", level = 5000),
            "5500" = threshold(conv.sq5$"160314", level = 5500),
            "6000" = threshold(conv.sq5$"160314", level = 6000),
            "6500" = threshold(conv.sq5$"160314", level = 6500))
            
# plot results
{
    for (i in names(th3)) {
        pdf(paste0(fpath, "sq3-thresholds-", i, ".pdf"))
            par(mar = c(2,2,3,1))
            image(c(1:1996), c(1:1996), th3[[i]], asp = T, col = c(NA, "blue"), main = i, xlab = "", ylab = "")
        dev.off()
    }
    
    for (i in names(th5)) {
        pdf(paste0(fpath, "sq5-thresholds-", i, ".pdf"))
        par(mar = c(2,2,3,1))
        image(c(1:1996), c(1:1996), th5[[i]], asp = T, col = c(NA, "red"), main = i, xlab = "", ylab = "")
        dev.off()
    }
}

# classification success

# mask for pos/neg
l.mask <- array(0, dim = c(1996, 1996))
l.mask[line.px] <- 1

rate3 <- rbind.fill(lapply(lapply(lapply(lapply(th3, table, l.mask), c), t), data.frame))
{
    colnames(rate3) <- c("TN", "FP", "FN", "TP"); rownames(rate3) <- names(th3)
    rate3$true.pos.rate <- rate3$TP / (rate3$TP + rate3$FN)
    rate3$false.disc <- rate3$FP / (rate3$TP + rate3$FP)
    rate3$false.neg <- rate3$FN / (rate3$TP + rate3$FN)
    rate3$precision <- rate3$TP / (rate3$TP + rate3$FP)
    rate3$accuracy <- (rate3$TP + rate3$TN) / 1996^2
    rate3$F1 <- 2 * rate3$TP / (2 * rate3$TP + rate3$FP + rate3$FN)
}

rate5 <- rbind.fill(lapply(lapply(lapply(lapply(th5, table, l.mask), c), t), data.frame))
{
    colnames(rate5) <- c("TN", "FP", "FN", "TP"); rownames(rate5) <- names(th5)
    rate5$true.pos.rate <- rate5$TP / (rate5$TP + rate5$FN)
    rate5$false.disc <- rate5$FP / (rate5$TP + rate5$FP)
    rate5$false.neg <- rate5$FN / (rate5$TP + rate5$FN)
    rate5$precision <- rate5$TP / (rate5$TP + rate5$FP)
}

rate3; rate5

####################################################################################################

# CONVOLVE WITH SINGLE-POINT KERNEL                                                             ####

# can this filter out small clusters before convolution with linear kernel?

# CONVOLVE WITH RECTANGULAR KERNEL                                                              ####


####################################################################################################

# LINE SMOOTHING                                                                                ####

# run a vertical filter of length 2k+1 over the thresholded image, and retain only points with value k+1

ll <- matrix(rep(1, 11), ncol = 1)

th5.3500.sm <- r2m(focal(m2r(th5$"3500"), ll))
th5.6000.sm <- r2m(focal(m2r(th5$"6000"), ll))


l.cols <- c(NA, rep("gold", 5), rep("cyan3", 6))

# smoothing over 5x5 convolution with 11-vector
{
    pdf(paste0(fpath, "11-smoothed-3500-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5.3500.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im1$xrng), 
              ylim = range(im1$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "11-smoothed-3500-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5.3500.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im2$xrng), 
              ylim = range(im2$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "11-smoothed-6000-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5.6000.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im1$xrng), 
              ylim = range(im1$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "11-smoothed-6000-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5.6000.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im2$xrng), 
              ylim = range(im2$yrng), col = l.cols)
        dev.off()
    }
}

# smoothing over 3x3 convolution with 7-vector
ll <- matrix(rep(1, 7), ncol = 1)
th3.1500.sm <- r2m(focal(m2r(th3$"1500"), ll))
th3.2000.sm <- r2m(focal(m2r(th3$"2000"), ll))
l.cols <- c(NA, rep("gold", 3), rep("cyan3", 8))

{
    pdf(paste0(fpath, "7-smoothed-1500-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3.1500.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im1$xrng), 
              ylim = range(im1$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "7-smoothed-1500-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3.1500.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im2$xrng), 
              ylim = range(im2$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "7-smoothed-2000-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3.2000.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im1$xrng), 
              ylim = range(im1$yrng), col = l.cols)
        dev.off()
    }
    pdf(paste0(fpath, "7-smoothed-2000-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3.2000.sm, breaks = c(0:12)-0.5, asp = T, xlim = range(im2$xrng), 
              ylim = range(im2$yrng), col = l.cols)
        dev.off()
    }
}

####################################################################################################



######### SCRATCH ######### ####                   
# thresholds
{
    xt.lin.b <- list()
    for (n in names(conv.lin.b)) {
        im <- conv.lin.b[[n]]
        xt.lin.b[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # linear kernel, black
{
    xt.lin.g <- list()
    for (n in names(conv.lin.g)) {
        im <- conv.lin.g[[n]]
        xt.lin.g[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # linear kernel, grey
{
    xt.sq3.b <- list()
    for (n in names(conv.sq3.b)) {
        im <- conv.sq3.b[[n]]
        xt.sq3.b[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # 3x3 kernel, black
{
    xt.sq3.g <- list()
    for (n in names(conv.sq3.g)) {
        im <- conv.sq3.g[[n]]
        xt.sq3.g[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # 3x3 kernel, grey
{
    xt.sq5.b <- list()
    for (n in names(conv.sq5.b)) {
        im <- conv.sq5.b[[n]]
        xt.sq5.b[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # 5x5 kernel, black
{
    xt.sq5.g <- list()
    for (n in names(conv.sq5.g)) {
        im <- conv.sq5.g[[n]]
        xt.sq5.g[[n]] <- matrix(findInterval(im, median(im, na.rm = T) + c(0:9) * mad(im, na.rm = T)), ncol = 1996)
    }
}   # 5x5 kernel, grey

####################################################################################################

# DISTRIBUTION OF EDGE PIXELS IN ALL IMAGES
{
    # how bright are the panel edges?
    # expect RHS of panel join to be higher in value on lower, LHS on upper
    edge.px <- matrix(c(rep(panel.edges()$x[1:16], 992), rep(panel.edges()$x[2:17]-1, 1004), 
                        sort(rep(1:992, 16)), sort(rep(993:1996, 16))), ncol = 2)
    
    xt.lin.b.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.lin.b, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    xt.lin.g.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.lin.g, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    xt.sq3.b.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq3.b, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    xt.sq3.g.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq3.g, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    xt.sq5.b.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq5.b, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    xt.sq5.g.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq5.g, "[", edge.px), table), as.matrix), t), data.frame))) / length(edge.px)
    
    xt.lin.b.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.lin.b, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    xt.lin.g.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.lin.g, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    xt.sq3.b.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq3.b, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    xt.sq3.g.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq3.g, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    xt.sq5.b.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq5.b, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    xt.sq5.g.line.tab <- as.matrix(rbind.fill(lapply(lapply(lapply(lapply(lapply(xt.sq5.g, "[", line.px), table), as.matrix), t), data.frame))) / length(line.px)
    
    # plots of proportion in each set
    {
        pdf(paste0(fpath, "props-lin-black.pdf"))
        plot(c(1:11)-0.1, xt.lin.b.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.lin.b.tab) + 0.01),
             main = "Linear kernel, black", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.lin.b.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.lin.b.tab)) {
            points(c(1:11)-0.1, xt.lin.b.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.lin.b.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # linear kernel, black
    {
        pdf(paste0(fpath, "props-lin-grey.pdf"))
        plot(c(1:11)-0.1, xt.lin.g.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.lin.g.tab) + 0.01),
             main = "Linear kernel, grey", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.lin.g.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.lin.g.tab)) {
            points(c(1:11)-0.1, xt.lin.g.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.lin.g.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # linear kernel, grey
    {
        pdf(paste0(fpath, "props-3x3-black.pdf"))
        plot(c(1:11)-0.1, xt.sq3.b.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.sq3.b.tab) + 0.01),
             main = "3x3 kernel, black", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.sq3.b.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.sq3.b.tab)) {
            points(c(1:11)-0.1, xt.sq3.b.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.sq3.b.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # 3x3 kernel, black
    {
        pdf(paste0(fpath, "props-3x3-grey.pdf"))
        plot(c(1:11)-0.1, xt.sq3.g.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.sq3.g.tab, na.rm = T) + 0.01),
             main = "3x3 kernel, grey", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.sq3.g.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.sq3.g.tab)) {
            points(c(1:11)-0.1, xt.sq3.g.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.sq3.g.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # 3x3 kernel, grey
    {
        pdf(paste0(fpath, "props-5x5-black.pdf"))
        plot(c(1:11)-0.1, xt.sq5.b.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.sq5.b.tab) + 0.01),
             main = "5x5 kernel, black", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.sq5.b.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.sq5.b.tab)) {
            points(c(1:11)-0.1, xt.sq5.b.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.sq5.b.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # 5x5 kernel, black
    {
        pdf(paste0(fpath, "props-5x5-grey.pdf"))
        plot(c(1:11)-0.1, xt.sq5.g.tab[1,], pch = 21, bg = th.cols, ylim = c(0, max(xt.sq5.g.tab) + 0.01),
             main = "5x5 kernel, grey", ylab = "Proportion of points", xlab = "MAD above median")
        points(c(1:11)+0.1, xt.sq5.g.line.tab[1,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        for (i in 2:nrow(xt.sq5.g.tab)) {
            points(c(1:11)-0.1, xt.sq5.g.tab[i,], pch = 21, bg = th.cols)
            points(c(1:11)+0.1, xt.sq5.g.line.tab[i,], pch = 20, col = adjustcolor(th.cols, alpha = 0.4))
        }
        dev.off()
    }   # 5x5 kernel, grey
    
    # plots of proportion in each set over time          # UNCLEAR
    {
        matplot(xt.lin.b.tab, bg = th.cols, col = "black", pch = 21, type = "o")
        matplot(xt.lin.b.line.tab, add = T, col = adjustcolor(th.cols, alpha = 0.4), pch = 20, type = "o")
    }
    
    # maybe look at change in proportion of pixels at each level that fall within each type?
    
}

####################################################################################################

# CLASSIFICATION COMPARISON
{
    # treat as classification problem: table of T/F pos/neg

    # combinations to check:
    # 3x3       at 5, 6, 7, 8, 9                expect 5.5/6.5?
    # 5x5       at 5, 6, 7, 8, 9                expect 6.5/7.5
}

# CLUMPING VS FILTER TO CLASSIFY CONTINUOUS LINES

# can't simply merge all points and assume that they form a line.
# single points will be smeared by into line by a 3x3 or 5x5 kernel.
m3 <- matrix(c(rep(0, 60), 1375, rep(0, 60)), ncol = 11)
r2m(focal(m2r(m3), k3))
r2m(focal(m2r(m3), k5))
r2m(focal(m2r(m3), matrix(c(rep(-1, 21), rep(6, 7), rep(-1, 21)), ncol = 7)))

#--------------------------------------------------------------------------------------
# clumping
cc <- clump(m2r(th3$"1500"))

zz <- freq(cc)[rev(order(freq(cc)[,2])),]

head(zz)

# tend to end up with short segments - better to smooth out first

# perhaps filter by segment length before this point?

#--------------------------------------------------------------------------------------
# focal window to identify lines: try 3, 5, 7, 9
l <- 7

ll <- matrix(rep(1, l), ncol = 1)

qq <- r2m(focal(m2r(th5$"5000"), ll))

pdf(paste0(fpath, "tmp.pdf"))
image(1:1996, 1:1996, threshold(qq, level = l/2), asp = T, col = c(NA, "red"))
dev.off()

# now, clump after smoothing
cc <- clump(m2r(qq))

v.rng <- list(x = c(100,200), y = c(1500,1600))
pixel.image(pw.m[,,"black", "160314"], xlim = v.rng$x, ylim = v.rng$y, title = "Raw black image")
image(1:1996, 1:1996, conv.sq5$"160314", main = "5x5 convolution", xlim = v.rng$x, ylim = v.rng$y, 
      breaks = c(-220000, 0, 3500, 4500, 5500, 6500, 632000), col = c(sd.colours()[1:6]), asp = T)
image(1:1996, 1:1996, th5$"5000", xlim = c(0,100), ylim = c(0,100), asp = T, main = "5x5 cut at 5000", col = c(NA, "green3"))
image(1:1996, 1:1996, qq, xlim = c(0,100), ylim = c(0,100), asp = T, main = "linear filter, length 7",
      breaks = c(0:8) - 0.5, col = c(NA, sd.colours()[1:7]))
image(1:1996, 1:1996, threshold(qq, level = l/2), add = T, xlim = c(0,100), main = "cut at 3.5", ylim = c(0,100), asp = T, col = c(NA, "red"))
#--------------------------------------------------------------------------------------
# classify final FP, FN, TP, TN and use to assess performance

# USE POINT FILTER TO EXCLUDE SINGLE BRIGHT POINTS                                        

# try applying rotated kernel & comparing result
{
    v.rng <- list(x = c(100,200), y = c(1500,1600))
    image(1:1996, 1:1996, conv.sq5$"160314", main = "5x5 convolution", xlim = v.rng$x, ylim = v.rng$y, 
          breaks = c(-220000, 0, 3500, 4500, 5500, 6500, 632000), col = c(sd.colours()[1:6]), asp = T)
    
    c2 <- r2m(focal(m2r(pw.m[,,"black", "160314"]), t(k5)))
    
    image(1:1996, 1:1996, c2, main = "5x5 convolution rotated", xlim = v.rng$x, ylim = v.rng$y, 
          breaks = c(-220000, 0, 3500, 4500, 5500, 6500, 632000), col = c(sd.colours()[1:6]), asp = T)
    
    # and where there's an actual line...
    v.rng <- list(x = c(400,500), y = c(1150,1250))
    
    # so. Use the difference?
    conv.diff <- conv.sq5$"160314" - c2
    
    pixel.image(conv.diff)
    
    o.plot(conv.diff[427,1150:1250], ylim = c(-5000,10000))
    
    # threshold as before
    
    image(1:1996, 1:1996, threshold(conv.diff, level = 3500), asp = T, col = c(NA, "magenta3"))
}


# we know that a single point will be smeared to a 3-point line by either kernel 
# can we just filter these out after thresholding?

# clump thresholded data
cc <- clump(m2r(th5$"6000"), dir = 4)

zz <- data.frame(freq(cc, useNA = "no"))
zz$new <- zz$value
zz$new[zz$count <= 3] <- NA

cc <- reclassify(cc, rcl = zz[,c(1, 3)])


pdf(paste0(fpath, "tmp.pdf"))
plot(cc, asp = T, col = "red"); dev.off()

# pdf images after convolution & thresholding
{
    pdf(paste0(fpath, "th-sq5-3500-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5$"3500", main = "", xlim = range(im1$xrng), ylim = range(im1$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    pdf(paste0(fpath, "th-sq5-6000-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5$"6000", main = "", xlim = range(im1$xrng), ylim = range(im1$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    
    pdf(paste0(fpath, "th-sq5-3500-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5$"3500", main = "", xlim = range(im2$xrng), ylim = range(im2$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    pdf(paste0(fpath, "th-sq5-6000-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th5$"6000", main = "", xlim = range(im2$xrng), ylim = range(im2$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    
    pdf(paste0(fpath, "th-sq3-1500-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3$"1500", main = "", xlim = range(im1$xrng), ylim = range(im1$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    pdf(paste0(fpath, "th-sq3-2000-upper.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3$"2000", main = "", xlim = range(im1$xrng), ylim = range(im1$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    
    pdf(paste0(fpath, "th-sq3-1500-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3$"1500", main = "", xlim = range(im2$xrng), ylim = range(im2$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
    pdf(paste0(fpath, "th-sq3-2000-lower.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, th3$"2000", main = "", xlim = range(im2$xrng), ylim = range(im2$yrng), 
              col = c(NA, "magenta3"), asp = T)
        dev.off()
    }
}


####################################################################################################
