
# FINDING BAD PIXELS BY CONVOLUTION

# try convolution with diagonal kernel (in both directions) to enhance the pattern

library("IO.Pixels"); library("CB.Misc")

fpath <- "./Notes/Point-convolution/fig/"
load.pixel.means()

im <- pw.m[,,, "160430"]                        # working with latest image only for now
r <- list(x = c(400:550), y = c(1150:1300))     # plotting range
focus <- get.focus(data.frame(x = c(427, 427), y = c(1199, 1200)), surround = 4)
interested <- matrix(c(rep(425, 2), rep(426, 2), rep(427, 3), rep(428, 2), rep(429, 2),
                       1199:1200, 1199:1200, 1198:1200, 1199:1200, 1199:1200), ncol = 2)

####################################################################################################

# INITIAL PLOTS                                                                                 ####

# plot raw images
pdf(paste0(fpath, "raw-data-image-upper-line.pdf")); {
    par(mar = c(2, 2, 1, 1))
    pixel.image(pw.m[,,"black", "160430"], xlim = range(r$x), ylim = range(r$y), panels = T)
    pixel.image(pw.m[,,"grey", "160430"], xlim = range(r$x), ylim = range(r$y), panels = T)
    pixel.image(pw.m[,,"white", "160430"], xlim = range(r$x), ylim = range(r$y), panels = T)
    dev.off()
}

# plot superclusters in raw images
pdf(paste0(fpath, "raw-data-sc.pdf")); {
    par(mar = c(2, 2, 1, 1), cex = 1.2)
    
    plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
    text(focus, labels = round(pw.m[,,"black", "160430"][focus]/1000, 1))
    
    plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
    text(focus, labels = round(pw.m[,,"grey", "160430"][focus]/1000, 1))
    
    plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
    text(focus, labels = round(pw.m[,,"white", "160430"][focus]/1000, 1))
    
    dev.off()
}

####################################################################################################

# CONVOLVE IMAGES WITH POINT KERNEL                                                             ####

# set kernels
{
    ks <- matrix(c(0, -1, 0, -1, 4, -1, 0, -1, 0), ncol = 3)/4    # 3x3 sharpening filter v1
    ks2 <- matrix(c(rep(-1, 4), 8, rep(-1, 4)), ncol = 3)/8               # 3x3 sharpening filter v2
    khp3 <- matrix(c(rep(-1, 4), 9, rep(-1, 4)), ncol = 3)              # 3x3 high-pass filter
    khp5 <- matrix(c(-1, -3, -4, -3, -1, -3, 0, 6, 0, -3, -4, 6, 21, 6, -4, -3, 0, 6, 0, -3, -1, -3, -4, -3, -1), ncol = 5)
    kd <- matrix(c(-1, 0, -1, 0, 4, 0, -1, 0, -1), ncol = 3)/4    # 3x3 diagonal
}

# convolution with smallest sharpening kernel (rook's dist)
{
    conv.s1 <- lapply(lapply(apply(im, 3, m2r), focal, ks), r2m)
    
    # produce plots
    {
        pdf(paste0(fpath, "s1-conv-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            pixel.image(conv.s1$black, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.s1$grey, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.s1$white, xlim = range(r$x), ylim = range(r$y), panels = T)
            dev.off()
        }
        
        pdf(paste0(fpath, "s1-th-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            image(1:1996, 1:1996, threshold(conv.s1$black, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s1$black, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$black, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$black, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.s1$grey, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s1$grey, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$grey, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$grey, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.s1$white, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s1$white, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$white, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s1$white, level = 600), col = c(NA, "blue"), add = T)
            dev.off()
        }
    }
    
    # check supercluster
    {
        pdf(paste0(fpath, "s1-conv-sc.pdf")); {
            par(mar = c(2, 2, 1, 1), cex = 1.2)
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s1$black[focus]/1000, 0))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s1$grey[focus]/1000, 1))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s1$white[focus]/1000, 1))
            
            dev.off()
        }
    }
}

# convolution with smallest sharpening kernel (queen's dist)
{
    conv.s2 <- lapply(lapply(apply(im, 3, m2r), focal, ks2), r2m)
    
    # produce plots
    {
        pdf(paste0(fpath, "s2-conv-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            pixel.image(conv.s2$black, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.s2$grey, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.s2$white, xlim = range(r$x), ylim = range(r$y), panels = T)
            dev.off()
        }
        
        pdf(paste0(fpath, "s2-th-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            image(1:1996, 1:1996, threshold(conv.s2$black, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s2$black, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$black, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$black, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.s2$grey, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s2$grey, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$grey, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$grey, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.s2$white, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.s2$white, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$white, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.s2$white, level = 600), col = c(NA, "blue"), add = T)
            dev.off()
        }
    }
    
    # check supercluster
    {
        pdf(paste0(fpath, "s2-conv-sc.pdf")); {
            par(mar = c(2, 2, 1, 1), cex = 1.2)
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s2$black[focus]/1000, 0))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s2$grey[focus]/1000, 1))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.s2$white[focus]/1000, 1))
            
            dev.off()
        }
    }
}

# convolution with diagonal kernel
{
    conv.diag <- lapply(lapply(apply(im, 3, m2r), focal, kd), r2m)
    
    # produce plots
    {
        pdf(paste0(fpath, "diag-conv-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            pixel.image(conv.diag$black, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.diag$grey, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.diag$white, xlim = range(r$x), ylim = range(r$y), panels = T)
            dev.off()
        }
        
        pdf(paste0(fpath, "diag-th-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            image(1:1996, 1:1996, threshold(conv.diag$black, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.diag$black, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$black, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$black, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.diag$grey, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.diag$grey, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$grey, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$grey, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.diag$white, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.diag$white, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$white, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.diag$white, level = 600), col = c(NA, "blue"), add = T)
            dev.off()
        }
    }
    
    # check supercluster
    {
        pdf(paste0(fpath, "diag-conv-sc.pdf")); {
            par(mar = c(2, 2, 1, 1), cex = 1.2)
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.diag$black[focus]/1000, 0))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.diag$grey[focus]/1000, 1))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.diag$white[focus]/1000, 1))
            
            dev.off()
        }
    }
}

# convolution with high-pass 5-kernel
{
    conv.hp5 <- lapply(lapply(apply(im, 3, m2r), focal, khp5), r2m)
    
    # produce plots
    {
        pdf(paste0(fpath, "hp5-conv-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            pixel.image(conv.hp5$black, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.hp5$grey, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(conv.hp5$white, xlim = range(r$x), ylim = range(r$y), panels = T)
            dev.off()
        }
        
        pdf(paste0(fpath, "hp5-th-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            image(1:1996, 1:1996, threshold(conv.hp5$black, level = 10000), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.hp5$black, level = 15000), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$black, level = 20000), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$black, level = 25000), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.hp5$grey, level = 25000), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.hp5$grey, level = 30000), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$grey, level = 35000), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$grey, level = 40000), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(conv.hp5$white, level = 70000), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(conv.hp5$white, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$white, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(conv.hp5$white, level = 600), col = c(NA, "blue"), add = T)
            dev.off()
        }
    }
    
    # check supercluster
    {
        pdf(paste0(fpath, "hp5-conv-sc.pdf")); {
            par(mar = c(2, 2, 1, 1), cex = 1.2)
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.hp5$black[focus]/1000, 0))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.hp5$grey[focus]/1000, 1))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(conv.hp5$white[focus]/1000, 1))
            
            dev.off()
        }
    }
}

# high-pass filter is weird. Do not like.

# try subtracting median filter from image & checking residuals
{
    conv.med <- lapply(lapply(apply(im, 3, m2r), focal, matrix(rep(1, 9), ncol = 3), fun = median), r2m)
    med.diff <- list(black = im[,,"black"] - conv.med$black,
                     grey = im[,,"grey"] - conv.med$grey,
                     white = im[,,"white"] - conv.med$white)
    
    # produce plots
    {
        pdf(paste0(fpath, "med-diff-conv-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            pixel.image(med.diff$black, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(med.diff$grey, xlim = range(r$x), ylim = range(r$y), panels = T)
            pixel.image(med.diff$white, xlim = range(r$x), ylim = range(r$y), panels = T)
            dev.off()
        }
        
        pdf(paste0(fpath, "med-diff-th-image-upper-line.pdf")); {
            par(mar = c(2, 2, 1, 1))
            image(1:1996, 1:1996, threshold(med.diff$black, level = 200), col = c(NA, "gold"), asp = T, xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(med.diff$black, level = 300), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$black, level = 400), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$black, level = 500), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(med.diff$grey, level = 300), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(med.diff$grey, level = 400), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$grey, level = 500), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$grey, level = 600), col = c(NA, "blue"), add = T)
            
            image(1:1996, 1:1996, threshold(med.diff$white, level = 400), col = c(NA, "gold"), xlim = range(r$x), ylim = range(r$y))
            image(1:1996, 1:1996, threshold(med.diff$white, level = 500), col = c(NA, "red"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$white, level = 600), col = c(NA, "magenta3"), add = T)
            image(1:1996, 1:1996, threshold(med.diff$white, level = 700), col = c(NA, "blue"), add = T)
            dev.off()
        }
    }
    
    # plot full image
    {
        pdf(paste0(fpath, "med-diff-black.pdf")); {
            pixel.image(med.diff$black, panels = T)
            dev.off()
        }
        pdf(paste0(fpath, "med-diff-grey.pdf")); {
            pixel.image(med.diff$grey, panels = T)
            dev.off()
        }
        pdf(paste0(fpath, "med-diff-white.pdf")); {
            pixel.image(med.diff$white, panels = T)
            dev.off()
        }
        
        # thresholded
        th <- 300
        pdf(paste0(fpath, "med-diff-th-", th, "-black.pdf")); {
            image(1:1996, 1:1996, threshold(med.diff$black, level = th), col = c(NA, "magenta3"))
            draw.panels()
            dev.off()
        }
        pdf(paste0(fpath, "med-diff-th-", th, "-grey")); {
            image(1:1996, 1:1996, threshold(med.diff$grey, level = th), col = c(NA, "magenta3"))
            draw.panels()
            dev.off()
        }
        pdf(paste0(fpath, "med-diff-th-", th, "-white.pdf")); {
            image(1:1996, 1:1996, threshold(med.diff$white, level = th), col = c(NA, "magenta3"))
            draw.panels()
            dev.off()
        }
    }

    # check supercluster
    {
        pdf(paste0(fpath, "med-diff-conv-sc.pdf")); {
            par(mar = c(2, 2, 1, 1), cex = 1.2)
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(med.diff$black[focus]/1000, 0))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(med.diff$grey[focus]/1000, 1))
            
            plot(interested, pch = 15, cex = 7, col = "gold", xlim = range(focus[,1]), ylim = range(focus[,2]), asp = T, xlab = "", ylab = "")
            text(focus, labels = round(med.diff$white[focus]/1000, 1))
            
            dev.off()
        }
    }
}

####################################################################################################

# THAT BLOODY DIAGONAL PATTERN IS BACK                                                          ####

th.w <- threshold(med.diff$white, level = 400)

# diagonal filter (sums to 0 for easier thresholding)
f.size <- 21; df <- matrix(rep(-1, f.size^2), ncol = f.size); diag(df) <- f.size - 1

pixel.image(r2m(focal(m2r(th.w), df)), title = paste0("Original filter \\, length ", f.size))
pixel.image(r2m(focal(m2r(th.w), df[,ncol(df):1])), title = paste0("Rotated filter /, length ", f.size))

conv.med <- lapply(lapply(apply(im, 3, m2r), focal, matrix(rep(1, 9), ncol = 3), fun = median), r2m)

md1 <- pw.m[,,"white", "141009"] - r2m(focal(m2r(pw.m[,,"white", "141009"]), matrix(rep(1, 9), ncol = 3), fun = median))
pdf("./Plots/Background-noise.pdf") {
    image(1:1996, 1:1996, threshold(md1, level = 500), col = c(NA, "red"))
    image(1:1996, 1:1996, threshold(pw.m[,,"white", "150529"] - r2m(focal(m2r(pw.m[,,"white", "141009"]), matrix(rep(1, 9), ncol = 3), fun = median)),
                                    level = 500), col = c(NA, adjustcolor("green3", alpha = 0.4)), add = F)
    
    image(1:1996, 1:1996, th.w, add = T, col = c(NA, adjustcolor("blue", alpha = 0.4)))
}



med.diff.white <- alply(pw.m[,,"white", ], 3, 
      function(x) x - r2m(focal(m2r(x), matrix(rep(1, 9), ncol = 3), fun = median)), .dims = T) 


