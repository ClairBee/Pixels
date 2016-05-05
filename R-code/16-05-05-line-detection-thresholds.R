
# DETECTION OF SLIGHTLY BRIGHT OR NON-RESPONSIVE LINES
# (creating plots for ./Notes/Line-detection/Line-detection.tex)

library("IO.Pixels"); library("CB.Misc")

load.pixel.means()
fpath <- "./Notes/Line-detection/fig/"

th.cols <- c("white", "paleturquoise1",  "cyan3", "skyblue", "green3",
             "gold", "orange", "red", "magenta3", "blue", "black")

# identify regions to look at in detail
im1 <- list(xrng = c(359:550), yrng = c(1100:1300), col = 427, row = 1199)
im2 <- list(xrng = c(766:896), yrng = c(50:200), col = 809)

line.px <- rbind(cbind(x = 427, y = c(1201:1996)),
                 cbind(x = 809, y = c(1:177)))

####################################################################################################
# specify kernels
k <- matrix(c(-1, 2, -1), nrow = 1)
k3 <- matrix(c(rep(-1, 3), rep(2, 3), rep(-1, 3)), ncol = 3)
k5 <- matrix(c(rep(-1, 10), rep(4, 5), rep(-1, 10)), ncol = 5)

# convolve all black & grey images
conv.lin.b <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k), r2m)
conv.lin.g <- lapply(lapply(apply(pw.m[,,"grey", ], 3, m2r), focal, k), r2m)

conv.sq3.b <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k3), r2m)
conv.sq3.g <- lapply(lapply(apply(pw.m[,,"grey", ], 3, m2r), focal, k3), r2m)

conv.sq5.b <- lapply(lapply(apply(pw.m[,,"black", ], 3, m2r), focal, k5), r2m)
conv.sq5.g <- lapply(lapply(apply(pw.m[,,"grey", ], 3, m2r), focal, k5), r2m)

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

# DISTRIBUTION OF EDGE PIXELS IN ALL IMAGES                                                     ####

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

####################################################################################################
# CLASSIFICATION COMPARISON                                                                     ####
# treat as classification problem: table of T/F pos/neg

# combinations to check:
# linear    at 4, 5, 6, 7 (black & grey)    expect 5.5?
# 3x3       at 5, 6, 7, 8                   expect 5.5/6.5?
# 5x5       at 5, 6, 7, 8, 9                expect 6.5/7.5
# sobel

