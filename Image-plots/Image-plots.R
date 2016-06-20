
library("IO.Pixels"); library("CB.Misc")

load.pixel.means()

####################################################################################################

# DIFFERENCES BETWEEN SUCCESSIVE ACQUISITIONS                                                   ####

diff <- pw.m[,,, 2:12] - pw.m[,,, 1:11]

for (col in dimnames(diff)[[3]]) {
    for (dt in dimnames(diff)[[4]]) {
        bmp(paste0("./Diffs-between-acquisitions/image-diffs-", col, "-", dt, ".bmp"), width = 1920, height = 1920)
            pixel.image(diff[ , , col, dt], title = paste0(col, " - ", dt))
        dev.off()
    }
}



####################################################################################################

# RAW IMAGES                                                                                    ####

for (col in dimnames(pw.m)[[3]]) {
    for (dt in dimnames(pw.m)[[4]]) {
        bmp(paste0("./Pixelwise-means/pw-mean-", col, "-", dt, ".bmp"), width = 1920, height = 1920)
        pixel.image(pw.m[ , , col, dt], title = paste0(col, " - ", dt))
        dev.off()
    }
}

####################################################################################################

# SHADING-CORRECTED IMAGES                                                                      ####

sc <- readRDS("../Other-data/Shading-corrections.rds")

sc[which(is.infinite(sc))] <- 0

for (dt in dimnames(sc)[[3]]) {
        bmp(paste0("./Shading-corrections/shading-correction-", col, "-", dt, ".bmp"), width = 1920, height = 1920)
        pixel.image(sc[ , , dt], title = paste0(col, " - ", dt))
        dev.off()
}

apply(sc, 3, mean, na.rm = T)

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
