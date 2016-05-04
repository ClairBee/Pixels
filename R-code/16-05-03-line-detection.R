
# DETECTION OF SLIGHTLY BRIGHT OR NON-RESPONSIVE LINES
# (creating plots for ./Notes/Line-detection/Line-detection.tex)

library("IO.Pixels")
load.pixel.means()
fpath <- "./Notes/Line-detection/fig/"

im <- pw.m[,,"black", "160314"]

# identify regions to look at in detail
im1 <- list(xrng = c(359:550), yrng = c(1100:1300), col = 427, row = 1199)
im2 <- list(xrng = c(766:896), yrng = c(512:640), col = 809)

####################################################################################################
# NUMERICAL SUMMARIES                                                                           ####

median(im[im1$col, (im1$row + 1):max(im1$yrng)]) - 
    median(im[im1$col, min(im1$yrng):(im1$row - 1)]) # 336.15

median(im[im1$col, (im1$row + 1):max(im1$yrng)] - im[im1$col - 1, (im1$row + 1):max(im1$yrng)]) # 315.95

median(im[im1$col, (im1$row + 1):max(im1$yrng)] - im[im1$col + 1, (im1$row + 1):max(im1$yrng)]) # 315.7

median(im[im2$col, im2$yrng] - im[im2$col + 1, im2$yrng])         # 330.2
median(im[im2$col, im2$yrng] - im[im2$col - 1, im2$yrng])         # 347.15

####################################################################################################

# PLOT RAW DATA                                                                                 ####

# subset of black image - upper panel
{
    pdf(paste0(fpath, "im-raw-upper.pdf"))
    par(mar = c(2,2,0,0))
        pixel.image(im, xlim = range(im1$xrng), ylim = range(im1$yrng))
    dev.off()
}

# subset of black image - lower panel
{
    pdf(paste0(fpath, "im-raw-lower.pdf"))
    par(mar = c(2,2,0,0))
        pixel.image(im, xlim = range(im2$xrng), ylim = range(im2$yrng))
    dev.off()
}

# transect along black image - upper panel
{
    pdf(paste0(fpath, "trans-raw-upper.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(im[im1$col, im1$yrng], ylim = c(4500,5500), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(im[im1$col - 1, im1$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(im[im1$col + 1, im1$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    
    abline(h = median(im[im1$col, min(im1$yrng):(im1$row - 1)]), col = "green3")
    abline(h = median(im[im1$col, (im1$row + 1):max(im1$yrng)]), col = "red")
    dev.off()
}

# transect along black image - lower panel
{
    pdf(paste0(fpath, "trans-raw-lower.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(im[im2$col, im2$yrng], ylim = c(4500,5500), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(im[im2$col - 1, im2$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(im[im2$col + 1, im2$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    
    abline(h = median(im[im2$col - 1, im2$yrng]), col = "green3")
    abline(h = median(im[im2$col, im2$yrng]), col = "red")
    dev.off()
}

mad(im)     # 283.992
sd(im)      # 541.2748

####################################################################################################

# LINEAR (-1, 2, -1) KERNEL                                                                     ####
k <- matrix(c(-1, 2, -1), nrow = 1)

# perform convolution & convert image back to array format (instead of raster)
conv.lin <- r2m(focal(m2r(im), k))              # k horizontal

# numerical summaries
{
    mad(conv.lin, na.rm = T)        # 123.8712
    sd(conv.lin, na.rm = T)         # 952.1302
    median(conv.lin, na.rm = T)     # -4.85
    mean(conv.lin, na.rm = T)       # 0.06597833
}

# upper panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-horiz-upper-im.pdf"))
    par(mar = c(2,2,0,0))
        pixel.image(conv.lin, xlim = range(im1$xrng), ylim = range(im1$yrng))
    dev.off()
}

# lower panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-horiz-lower-im.pdf"))
    par(mar = c(2,2,0,0))
        pixel.image(conv.lin, xlim = range(im2$xrng), ylim = range(im2$yrng))
    dev.off()
}

# transect across upper panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-horiz-upper-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.lin[im1$col, im1$yrng], ylim = c(-1000,1000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.lin[im1$col - 1, im1$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.lin[im1$col + 1, im1$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.lin[im1$col + 2, im1$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.lin[im1$col - 2, im1$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.lin, na.rm = T) + mad(conv.lin, na.rm = T) * c(1,2), col = "blue", lty = 2)
    dev.off()
}

# transect across lower panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-horiz-lower-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.lin[im2$col, im2$yrng], ylim = c(-1000,1000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.lin[im2$col - 1, im2$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.lin[im2$col + 1, im2$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.lin[im2$col + 2, im2$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.lin[im2$col - 2, im2$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.lin, na.rm = T) + mad(conv.lin, na.rm = T) * c(1,2), col = "blue", lty = 2)
    dev.off()
}

####################################################################################################

# LINEAR(1, 1, 1) KERNEL                                                                        ####

k2 <- matrix(c(1, 1, 1), ncol = 1)

# perform convolution & convert image back to array format (instead of raster)
conv.lin.solid <- r2m(focal(m2r(im), k2))              # k2 vertical

# numerical summaries
{
    mad(conv.lin.solid, na.rm = T)        # 833.666
    sd(conv.lin.solid, na.rm = T)         # 1266.388
    median(conv.lin.solid, na.rm = T)     # 15110.85
    mean(conv.lin.solid, na.rm = T)       # 15299.78
}

# upper panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-solid-horiz-upper-im.pdf"))
    par(mar = c(2,2,0,0))
    pixel.image(conv.lin.solid, xlim = range(im1$xrng), ylim = range(im1$yrng))
    dev.off()
}

# lower panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-solid-horiz-lower-im.pdf"))
    par(mar = c(2,2,0,0))
    pixel.image(conv.lin.solid, xlim = range(im2$xrng), ylim = range(im2$yrng))
    dev.off()
}

# transect across upper panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-solid-horiz-upper-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.lin.solid[im1$col, im1$yrng], ylim = c(13000, 16000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.lin.solid[im1$col - 1, im1$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.lin.solid[im1$col + 1, im1$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.lin.solid[im1$col + 2, im1$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.lin.solid[im1$col - 2, im1$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.lin.solid, na.rm = T) + mad(conv.lin.solid, na.rm = T) * c(-2, -1, 0, 1,2), col = "blue", lty = 2)
    dev.off()
}

# transect across lower panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-lin-solid-horiz-lower-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.lin.solid[im2$col, im2$yrng], ylim = c(13000, 16000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.lin.solid[im2$col - 1, im2$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.lin.solid[im2$col + 1, im2$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.lin.solid[im2$col + 2, im2$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.lin.solid[im2$col - 2, im2$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.lin.solid, na.rm = T) + mad(conv.lin.solid, na.rm = T) * c(-2, -1, 0, 1,2), col = "blue", lty = 2)
    dev.off()
}

####################################################################################################

# SQUARE KERNEL                                                                                 ####

k3 <- matrix(c(rep(-1, 3), rep(2, 3), rep(-1, 3)), ncol = 3)

# perform convolution & convert image back to array format (instead of raster)
conv.sq <- r2m(focal(m2r(im), k3))              # k3 vertical

# numerical summaries
{
    mad(conv.sq, na.rm = T)        # 1689.211
    sd(conv.sq, na.rm = T)         # 1689.211
    median(conv.sq, na.rm = T)     # -11.05
    mean(conv.sq, na.rm = T)       # 0.1980
}

# upper panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-sq-horiz-upper-im.pdf"))
    par(mar = c(2,2,0,0))
    pixel.image(conv.sq, xlim = range(im1$xrng), ylim = range(im1$yrng))
    dev.off()
}

# lower panel, after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-sq-horiz-lower-im.pdf"))
    par(mar = c(2,2,0,0))
    pixel.image(conv.sq, xlim = range(im2$xrng), ylim = range(im2$yrng))
    dev.off()
}

# transect across upper panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-sq-horiz-upper-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.sq[im1$col, im1$yrng], ylim = c(-2000,3000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.sq[im1$col - 1, im1$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.sq[im1$col + 1, im1$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.sq[im1$col + 2, im1$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.sq[im1$col - 2, im1$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.sq, na.rm = T) + mad(conv.sq, na.rm = T) * c(-3:3), col = "blue", lty = 2)
    dev.off()
}

# transect across lower panel after convolution with horizontal kernel
{
    pdf(paste0(fpath, "conv-sq-horiz-lower-trans.pdf"))
    par(mar = c(2,2,0,0))
    
    o.plot(conv.sq[im2$col, im2$yrng], ylim = c(-2000,3000), xlab = "Transect sample: outer edge ->", ylab = "Pixel value")
    o.plot(conv.sq[im2$col - 1, im2$yrng], add = T, col = adjustcolor("skyblue", alpha = 0.4))
    o.plot(conv.sq[im2$col + 1, im2$yrng], add = T, col = adjustcolor("gold", alpha = 0.4))
    o.plot(conv.sq[im2$col + 2, im2$yrng], add = T, col = adjustcolor("orange", alpha = 0.4))
    o.plot(conv.sq[im2$col - 2, im2$yrng], add = T, col = adjustcolor("green3", alpha = 0.4))
    
    abline(h = median(conv.sq, na.rm = T) + mad(conv.sq, na.rm = T) * c(-2, -1, 0, 1,2), col = "blue", lty = 2)
    dev.off()
}

####################################################################################################

# CLASSIFICATION BY SQUARE KERNEL                                                               ####

# threshold to find only high points
xt <- matrix(findInterval(conv.sq, 
                          median(conv.sq, na.rm = T) + c(2:9) * mad(conv.sq, na.rm = T)), 
             ncol = 1996)

# plot thresholds with & without panel lines
{
    pdf(paste0(fpath, "conv-sq-mad-thresholds.pdf"))
    par(mar = c(2,2,0,0))
    image(c(1:1996), c(1:1996), xt, col = c("white", "green", "blue", "red"), 
          breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5), xlab = "", ylab = "")
    dev.off()
    
    pdf(paste0(fpath, "conv-sq-mad-thresholds-with-panels.pdf"))
    par(mar = c(2,2,0,0))
    image(c(1:1996), c(1:1996), xt, col = c("white", "green", "blue", "red"), 
          breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5), xlab = "", ylab = ""); draw.panels()
    dev.off()
}

# use a second focal window to identify line segments
mad2 <- matrix(findInterval(conv.sq, 
                                median(conv.sq, na.rm = T) + 2 * mad(conv.sq, na.rm = T)), 
                   ncol = 1996)
mad2.lines <- r2m(focal(m2r(mad2), w = matrix(rep(1, 3), ncol = 1)))

mad2.lines[mad2.lines < 3] <- 0
image(mad2.lines, col = c("white", "red"))
draw.panels()

image(c(1:1996), c(1:1996), r2m(mad2.lines), col = c("white", "green3", "blue", "red"), 
      breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5), xlab = "", ylab = "", 
      xlim = im1$col + c(-10,10), ylim = im1$row + c(-10,10), asp = T)
