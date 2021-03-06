
library("IO.Pixels"); library("CB.Misc")

pw.m <- load.pixel.means(c("loan", "MCT225"))

fpath <- "./Image-plots/"

####################################################################################################

# DIFFERENCES BETWEEN SUCCESSIVE ACQUISITIONS                                                   ####

diff <- pw.m[,,, 2:dim(pw.m)[[4]]] - pw.m[,,, 1:(dim(pw.m)[[4]]-1)]

for (col in dimnames(diff)[[3]]) {
    for (dt in dimnames(diff)[[4]]) {
        bmp(paste0("./Diffs-between-acquisitions/image-diffs-", col, "-", dt, ".bmp"), width = 2048, height = 2048)
            pixel.image(diff[ , , col, dt], title = paste0(col, " - ", dt))
        dev.off()
    }
}



####################################################################################################

# RAW IMAGES                                                                                    ####

# pixelwise mean images
{
    lapply(dimnames(pw.m)[[4]],
           function(acq) {
               lapply(dimnames(pw.m)[[3]], 
                      function(cc) {
                          bmp(paste0(fpath, "Pixelwise-means/pw-mean-", cc, "-", acq, ".bmp"), width = 2048, height = 2048)
                          pixel.image(pw.m[ , , cc, acq], title = paste0(cc, " - ", acq))
                          dev.off()
                      })
           })
}

# histograms of raw images
{
    lapply(dimnames(pw.m)[[4]],
           function(acq) {
               lapply(dimnames(pw.m)[[3]], 
                      function(cc) {
                          cb.out("bmp", paste0(fpath, "Pixelwise-means/pw-mean-hist-", cc, "-", acq), width = 2048, height = 1024)
                          hist(pw.m[ , , cc, acq], breaks = "fd", title = paste0(cc, " - ", acq), main = "", col = "black", xlab = "", ylab = "", xlim = c(0,65535))
                          dev.off()
                          cb.out("bmp", paste0(fpath, "Pixelwise-means/pw-mean-hist-cropped-", cc, "-", acq), width = 2048, height = 1024)
                          hist(pw.m[ , , cc, acq], breaks = "fd", title = paste0(cc, " - ", acq), col = "black", ylim = c(0,30), main = "", xlab = "", ylab = "", xlim = c(0,65535))
                          dev.off()
                      })
           })
}




####################################################################################################

# SHADING-CORRECTED IMAGES                                                                      ####

for (dt in dimnames(pw.m)[[4]]) {
        bmp(paste0("./Shading-corrections/shading-correction-", dt, ".bmp"), width = 2048, height = 2048)
        pixel.image(shading.corrected(pw.m[,,, dt]), title = paste0(col, " - ", dt))
        dev.off()
}


####################################################################################################

# SHADING CORRECTION BOXPLOTS                                                                   ####

    for (dt in dimnames(sc)[[3]]) {
        bmp(paste0("./Shading-correction-boxplots/sc-boxplot-", col, "-", dt, ".bmp"), width = 1920, height = 1920)
        par(mfrow = c(1, 2))
        sp <- subpanels(sc[ , , dt])
        dimnames(sp) <- list(NULL, NULL, apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                                     rep(c(rep("0", 9), rep("", 7)), 2),
                                                     rep(c(1:16), 2)), 1, paste, collapse = ""))
        
        boxplot(value ~ X3, melt(sp[,,1:16]), pch = 20, ylim = c(0, 65535))
        abline(h = c(mean(sc[ , , dt]), 
                     median(sc[ , , dt]) + (65535 - median(sc[ , , dt])) / 4 * c(1, 2), 65535), 
               col = c("cyan3", "orange", "red", "magenta3"), lty = c(1, 2, 2, 2))
        boxplot(value ~ X3, melt(sp[,,17:32]), pch = 20, ylim = c(0, 65535))
        abline(h = c(mean(sc[ , , dt]), 
                     median(sc[ , , dt]) + (65535 - median(sc[ , , dt])) / 4 * c(1, 2), 65535), 
               col = c("cyan3", "orange", "red", "magenta3"), lty = c(1, 2, 2, 2))
        dev.off()
    }


####################################################################################################

# SUBPANEL BOXPLOTS                                                                             ####

for (col in dimnames(pw.m)[[3]]) {
    for (dt in dimnames(pw.m)[[4]]) {
        bmp(paste0("./Boxplot-mean/pw-mean-boxplot-", col, "-", dt, ".bmp"), width = 1920, height = 1920)
        par(mfrow = c(1, 2))
        sp <- subpanels(pw.m[ , , col, dt])
        dimnames(sp) <- list(NULL, NULL, apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                                     rep(c(rep("0", 9), rep("", 7)), 2),
                                                     rep(c(1:16), 2)), 1, paste, collapse = ""))
        
        boxplot(value ~ X3, melt(sp[,,1:16]), pch = 20, ylim = c(0, 65535))
        abline(h = c(mean(pw.m[ , , col, dt]), 
                          median(pw.m[ , , col, dt]) + (65535 - median(pw.m[ , , col, dt])) / 4 * c(1, 2), 65535), 
                     col = c("cyan3", "orange", "red", "magenta3"), lty = c(1, 2, 2, 2))
        boxplot(value ~ X3, melt(sp[,,17:32]), pch = 20, ylim = c(0, 65535))
        abline(h = c(mean(pw.m[ , , col, dt]), 
                     median(pw.m[ , , col, dt]) + (65535 - median(pw.m[ , , col, dt])) / 4 * c(1, 2), 65535), 
               col = c("cyan3", "orange", "red", "magenta3"), lty = c(1, 2, 2, 2))
        dev.off()
    }
}
