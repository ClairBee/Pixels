
library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)
fpath <- "./Notes/MCT225/fig/"

.median <- hijack(median, na.rm = TRUE)

####################################################################################################

# XML PROFILE SUMMARY                                                                           ####

c(apply(pw.m, 3:4, mean, na.rm = T))

err <- merge(data.frame(which(pw.m[,,"grey", "MCT225"] < 10000, arr.ind = T), src = "g"),
             data.frame(which(pw.m[,,"white", "MCT225"] < 20000, arr.ind = T), src = "w"),
             by = c(1:2), all = T)
table(err$src.x, err$src.y, useNA = "ifany") # same pixels in both: 3254 of all

active.pwm <- pw.m[,,, c("160430", "MCT225")]; {
    active.pwm[,,"black", "160430"][which(pw.m[,,"white", "160430"] < 10000)] <- NA
    active.pwm[,,"grey", "160430"][which(pw.m[,,"white", "160430"] < 10000)] <- NA
    active.pwm[,,"white", "160430"][which(pw.m[,,"white", "160430"] < 10000)] <- NA
    
    active.pwm[,,"black", "MCT225"][which(pw.m[,,"white", "MCT225"] < 20000)] <- NA
    active.pwm[,,"grey", "MCT225"][which(pw.m[,,"white", "MCT225"] < 20000)] <- NA
    active.pwm[,,"white", "MCT225"][which(pw.m[,,"white", "MCT225"] < 20000)] <- NA
}

df <- summarise.profiles()
zz <- cbind(setNames(df[df$date %in% c("160430", "MCT225"), 
                        c("date", "batch", "kV", "uA", "Power", "ExposureTime",
                          "CameraProperties..attrs.gain", "ImageDetector")],
                     nm = c("Date", "Image", "kV", "uA", "Power", "Exposure", "Gain", "Model")), 
            "Mean value" = round(c(apply(active.pwm, 3:4, mean, na.rm = T)), 0),
            "Median value" = round(c(apply(active.pwm, 3:4, median, na.rm = T)), 0),
            "SD (active px)" = round(c(apply(active.pwm, 3:4, sd, na.rm = T)), 0))

write.csv(zz[,-c(1:2)], paste0(fpath, "Image-summary.csv"), quote = F, 
          row.names = apply(zz[,1:2], 1, paste, collapse = " - "))

####################################################################################################

# EXPLORATORY                                                                                   ####

# images
jpeg(paste0(fpath, "pw-mean-black.jpg")); {
    par(mar = c(2, 2, 1, 1))
    pixel.image(pw.m[,,"black", "MCT225"])
    dev.off()
}
jpeg(paste0(fpath, "pw-mean-grey.jpg")); {
    par(mar = c(2, 2, 1, 1))
    pixel.image(pw.m[,,"grey", "MCT225"])
    dev.off()
}
jpeg(paste0(fpath, "pw-mean-white.jpg")); {
    par(mar = c(2, 2, 1, 1))
    pixel.image(pw.m[,,"white", "MCT225"])
    dev.off()
}

####################################################################################################

# HISTOGRAMS                                                                                    ####

th.hist <- function(im, crop = NA) {
    par(mar = c(2,2,1,1))
    if (is.na(crop)) {
        hh <- hist(im, breaks = "fd", main = "", xlab = "", ylab = "", xlim = c(0,65535))
    } else {
        hh <- hist(im, breaks = "fd", main = "", xlab = "", ylab = "", ylim = c(0, crop), xlim = c(0,65535))
    }
    
    ym <- max(pretty(hh$counts)) * 1.5
    
    med <- median(im, na.rm = T)
    rect(med + (65535 - med)/2, 0, 65535, ym, col = adjustcolor("red", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, ym, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(0, 0, med /2, ym, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
    rect(med /2, 0, med * 0.75, ym, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
}

pdf(paste0(fpath, "MCT225-hist-black.pdf"), height = 4); {
    th.hist(pw.m[,,"black", "MCT225"])
    dev.off()
}
pdf(paste0(fpath, "MCT225-hist-grey.pdf"), height = 4); {
th.hist(pw.m[,,"grey", "MCT225"])
    dev.off()
}
pdf(paste0(fpath, "MCT225-hist-white.pdf"), height = 4); {
th.hist(pw.m[,,"white", "MCT225"])
    dev.off()
}

pdf(paste0(fpath, "MCT225-hist-black-cropped.pdf"), height = 4); {

th.hist(pw.m[,,"black", "MCT225"], crop = 30)
    dev.off()
}
pdf(paste0(fpath, "MCT225-hist-grey-cropped.pdf"), height = 4); {
th.hist(pw.m[,,"grey", "MCT225"], crop = 30)
    dev.off()
}
pdf(paste0(fpath, "MCT225-hist-white-cropped.pdf"), height = 4); {
th.hist(pw.m[,,"white", "MCT225"], crop = 30)
    dev.off()
}


####################################################################################################

# DARK PIXELS: NON-RESPONSIVE/BLOCKED                                                           ####

# probably necessary to do this first: can then exclude when defining 'normal' pixel behaviour
hist(pw.m[41:2028,41:2028,"grey", "160430"] - pw.m[41:2028,41:2028,"black", "160430"], 
     breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = ""); abline(v = 5000, col = "red", lty = 2)
hist(pw.m[41:2028,41:2028,"grey", "MCT225"] - pw.m[41:2028,41:2028,"black", "MCT225"], 
     breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = ""); abline(v = 5000, col = "red", lty = 2)

hist(pw.m[41:2028,41:2028,"white", "MCT225"] - pw.m[41:2028,41:2028,"black", "MCT225"], 
     breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = ""); abline(v = 5000, col = "red", lty = 2)
hist(pw.m[41:2028,41:2028,"white", "160430"] - pw.m[41:2028,41:2028,"black", "160430"], 
     breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = ""); abline(v = 5000, col = "red", lty = 2)

bpx <- list("MCT225" = data.frame(which(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] < 5000, arr.ind = T), 
                                type = "dark"),
            "160430" = data.frame(which(pw.m[,,"grey", "160430"] - pw.m[,,"black", "160430"] < 5000, arr.ind = T), 
                                 type = "dark"),
            "131122" = data.frame(which(pw.m[,,"grey", "131122"] - pw.m[,,"black", "131122"] < 5000, arr.ind = T), 
                                 type = "dark"))

plot(bpx$"MCT225"[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0,2048))
plot(bpx$"160430"[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0,2048))
plot(bpx$"131122"[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0,2048))

####################################################################################################

# OTHER COLUMN DEFECTS                                                                          ####

ll <- list("160430" = rbind(which(find.lines(pw.m[,,"grey", "160430"]) > 0, arr.ind = T),
                            which(find.lines(pw.m[,,"grey", "160430"], dim.lines = T) > 0, arr.ind = T)),
           "MCT225" = rbind(which(find.lines(pw.m[,,"grey", "MCT225"]) > 0, arr.ind = T),
                            which(find.lines(pw.m[,,"grey", "MCT225"], dim.lines = T) > 0, arr.ind = T)),
           "131122" = rbind(which(find.lines(pw.m[,,"grey", "131122"]) > 0, arr.ind = T),
                            which(find.lines(pw.m[,,"grey", "131122"], dim.lines = T) > 0, arr.ind = T)))
           
plot(ll$"MCT225", xlim = c(0,2048), ylim = c(0,2048), pch = ".")

pixel.image(pw.m[,,"grey", "MCT225"], xlim = c(240, 270), ylim = c(1870, 1900))
points(ll$"MCT225", xlim = c(240, 270), ylim = c(1850, 1900), pch = 0)

le.plot(256, 1880)
le.plot(448, 1137)
le.plot(882, 1119)
le.plot(1691, 1359)
le.plot(1926, 1627)

# create bad pixel map manually for now

####################################################################################################

# LOCALLY BRIGHT/DIM PIXELS                                                                     ####

# replace dark pixels with median values
med.replace <- function(im, px, w = 5) {
    
    get.px <- cbind(px[1] + c(-floor(w/2):floor(w/2)), px[2])
    median(im[get.px], na.rm = T)
    
}

dark.px <- which(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] < 5000, arr.ind = T)
mr <- apply(dark.px, 1, med.replace, im = pw.m[,,"grey", "MCT225"])

grey.adj <- pw.m[,,"grey", "MCT225"]
grey.adj[dark.px] <- mr
grey.adj.md <- med.diffs(grey.adj)

hist(grey.adj.md, breaks = "fd", ylim = c(0,30))

# median differences - immediate neighbourhood
md <- abind("MCT225" = readRDS("./02_Objects/med-diffs/md-MCT225.rds"), 
            "160430" = readRDS("./02_Objects/med-diffs/md-160430.rds"),
            along = 4)

hist(md[,,"black", "160430"], breaks = "fd", ylim = c(0,30))
hist(md[,,"black", "MCT225"], breaks = "fd", ylim = c(0,30))

hist(md[,,"grey", "160430"], breaks = "fd", ylim = c(0,30))
hist(md[,,"grey", "MCT225"], breaks = "fd", ylim = c(0,30))

# try thresholding at 1000
local.px <- which(abs(grey.adj.md) > 1000, arr.ind = T)
plot(which(abs(grey.adj.md) > 1000, arr.ind = T), pch = ".", xlim = c(0,2048), ylim = c(0,2048))
points(which(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] < 5000, arr.ind = T), col = "red", pch = ".")

# median-difference plots of ring artefacts
{
    focal.plot(md[,,"grey", "MCT225"], centre = c(892, 345), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"grey", "MCT225"], centre = c(356, 1829), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"grey", "MCT225"], centre = c(1208, 1955), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"grey", "MCT225"], centre = c(105, 640), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"grey", "MCT225"], centre = c(236, 271), bad.px = bpx.MCT225, bpx.cex = 5)
    
    focal.plot(md[,,"white", "MCT225"], centre = c(892, 345), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"white", "MCT225"], centre = c(356, 1829), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"white", "MCT225"], centre = c(1208, 1955), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"white", "MCT225"], centre = c(105, 640), bad.px = bpx.MCT225, bpx.cex = 5)
    focal.plot(md[,,"white", "MCT225"], centre = c(236, 271), bad.px = bpx.MCT225, bpx.cex = 5)
}

lim <- 1000
length(which(abs(md[,,"black", "160430"]) > lim))
length(which(abs(md[,,"grey", "160430"]) > lim))
length(which(abs(md[,,"white", "160430"]) > lim))
length(which(abs(grey.adj.md) > lim, arr.ind = T))
length(which(abs(md[,,"black", "MCT225"]) > lim))
length(which(abs(md[,,"grey", "MCT225"]) > lim))
length(which(abs(md[,,"white", "MCT225"]) > lim))

####################################################################################################

# INTERSECTION OF LOCAL PX WITH NONLINEAR?                                                      ####

wlm.MCT225 <- fit.w.lm(pw.m[,,,"MCT225"])

smoothScatter(wlm.MCT225$df$w, wlm.MCT225$df$fv,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
              xlab = "Observed", ylab = "Fitted value")
abline(line(wlm.MCT225$df$w, wlm.MCT225$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)

smoothScatter(wlm.MCT225$df$fv, wlm.MCT225$df$res,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)

#--------------------------------------
# REMOVE DARK PIXELS

active.MCT225 <- pw.m[,,,"MCT225"]
active.MCT225[,,"black"][dark.px] <- NA
active.MCT225[,,"grey"][dark.px] <- NA
active.MCT225[,,"white"][dark.px] <- NA

wlm.active <- fit.w.lm(active.MCT225)

smoothScatter(wlm.active$df$w, wlm.active$df$fv,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
              xlab = "Observed", ylab = "Fitted value")
abline(line(wlm.active$df$w, wlm.active$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)

smoothScatter(wlm.active$df$fv, wlm.active$df$res,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)

#--------------------------------------
# REMOVE EXTREME PIXELS

apply(active.MCT225, 3, 
      function(im) {
          hist(im, breaks = "fd", ylim = c(0,30), xlab = "Pixelwise mean (grey values)", ylab = "Frequency", main = "")
          med <- median(im, na.rm = T)
          rect(med + (65535 - med)/2, 0, 65535, 60000, col = adjustcolor("red", alpha = 0.3), border = NA)
          rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
          rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, 60000, col = adjustcolor("orange", alpha = 0.3), border = NA)
          rect(0, 0, med /2, 60000, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
          rect(med /2, 0, med * 0.75, 60000, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
      })


extreme.values <- function(im) {
    
    med.b <- median(im[,,"black"], na.rm = T)
    xt.b <- which(im[,,"black"] < med.b * 0.75 | im[,,"black"] > med.b + (65535 - med.b)/4, arr.ind = T)

    med.g <- median(im[,,"grey"], na.rm = T)
    xt.g <- which(im[,,"grey"] < med.g * 0.75 | im[,,"grey"] > med.g + (65535 - med.g)/4, arr.ind = T)
    
    med.w <- median(im[,,"white"], na.rm = T)
    xt.w <- which(im[,,"white"] < med.w * 0.75 | im[,,"white"] > med.w + (65535 - med.w)/4, arr.ind = T)
    
    unique(rbind(xt.b, xt.g, xt.w))
}

# 22 globally extreme pixels. Woo.
extreme <- extreme.values(active.MCT225)

global.MCT225 <- active.MCT225
global.MCT225[,,"black"][extreme] <- NA
global.MCT225[,,"grey"][extreme] <- NA
global.MCT225[,,"white"][extreme] <- NA

wlm.global <- fit.w.lm(global.MCT225)

smoothScatter(wlm.global$df$w, wlm.global$df$fv,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
              xlab = "Observed", ylab = "Fitted value")
abline(line(wlm.global$df$w, wlm.global$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)

smoothScatter(wlm.global$df$fv, wlm.global$df$res,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)

#--------------------------------------
# REMOVE LOCAL PIXELS

local.px <- unique(rbind(which(abs(md[,,"black", "MCT225"]) > 1000, arr.ind = T),
                         which(abs(md[,,"grey", "MCT225"]) > 1000, arr.ind = T),
                         which(abs(md[,,"white", "MCT225"]) > 1000, arr.ind = T)))
      
local.MCT225 <- global.MCT225
local.MCT225[,,"black"][local.px] <- NA
local.MCT225[,,"grey"][local.px] <- NA
local.MCT225[,,"white"][local.px] <- NA

wlm.local <- fit.w.lm(local.MCT225)

smoothScatter(wlm.local$df$w, wlm.local$df$fv,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
              xlab = "Observed", ylab = "Fitted value")
abline(line(wlm.local$df$w, wlm.local$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)

smoothScatter(wlm.local$df$fv, wlm.local$df$res,
              nrpoints = 0, xlim = c(0,65535),
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))
abline(0,0, col = adjustcolor("darkred", alpha = 0.4), lty = 2)

#--------------------------------------
# WHAT'S LEFT?

wres.px <- wlm.local$df[which(abs(wlm.local$df$res) > 500),]

plot(wres.px[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0, 2048))
points(dark.px, col = "red", pch = ".")

# remaining problematic pixels seem to be neighbours of dead lines
plot(wres.px[,1:2], pch = 15, col = "red", cex = 0.7, xlim = c(1900, 2000), ylim = c(1800, 2048))
points(dark.px, pch = 15, cex = 0.7)

plot(wres.px[,1:2], pch = 15, col = "red", cex = 0.7, xlim = c(400, 500), ylim = c(1000, 1300))
points(dark.px, pch = 15, cex = 0.7)

####################################################################################################

# WHITE RESIDUAL VS MEDIAN DIFFERENCE                                                           ####

df <- wlm.MCT225$df
df$md.b <- md[,,"black", "MCT225"][as.matrix(df[,1:2])]
df$md.g <- md[,,"grey", "MCT225"][as.matrix(df[,1:2])]
df$md.w <- md[,,"white", "MCT225"][as.matrix(df[,1:2])]

smoothScatter(df$md.g, df$res, nrpoints = 0,
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))

smoothScatter(df$md.g, df$res, nrpoints = 0,
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))

smoothScatter(df$md.g[abs(df$res) <= 500], df$md.b[abs(df$res) <= 500], , nrpoints = 0,
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
              nbin = 256)

####################################################################################################

# COMPARE ENDS OF DEAD LINES                                                                    ####

pixel.image(pw.m[,,"white", "MCT225"], xlim = c(200, 300), ylim = c(1800, 1900))

# column 448:449
{
    cc <- 447; xl <- c(1000, 1300)
    # column plots (white)
    {
        plot(pw.m[cc,,"white", "MCT225"], type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+2,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
        
        legend("right", col = c("black", adjustcolor(c("blue", "red", "green"), alpha = 0.4)), lty = 1,
               legend = paste("Column ", cc + c(0:3), sep = ""), bty = "n")
    }
    
    # power plots
    {
        plot(pw.m[cc,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc,,"black", "MCT225"])
        
        plot(pw.m[cc+1,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+1,,"black", "MCT225"])
        
        plot(pw.m[cc+2,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+2,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+2,,"black", "MCT225"])
        
        plot(pw.m[cc+3,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc+3,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+3,,"black", "MCT225"])
    }
    
    # neighbour plots
    {
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], ylim = c(-50000, 3000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
        
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], xlim = c(0,1000), ylim = c(-1000, 1000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
        
        plot(pw.m[cc+1,,"black", "MCT225"] - pw.m[cc,,"black", "MCT225"], xlim = c(0,1000), ylim = c(-500, 500), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"black", "MCT225"] - pw.m[cc+1,,"black", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"black", "MCT225"] - pw.m[cc+2,,"black", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"black", "MCT225"] - pw.m[cc+3,,"black", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
        
        mean(pw.m[cc+c(0:4),0:1000,"black", "MCT225"], na.rm = T)
    }
}

# column 255:257
{
    cc <- 255; xl <- c(1700, 2048)
    # column plots (white)
    {
        plot(pw.m[cc,,"white", "MCT225"], type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc+1,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+2,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
        
        legend("right", col = c("black", adjustcolor(c("blue", "red", "green"), alpha = 0.4)), lty = 1,
               legend = paste("Column ", cc + c(0:3), sep = ""), bty = "n")
    }
    
    # power plots
    {
        plot(pw.m[cc,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc,,"black", "MCT225"])
        
        plot(pw.m[cc+1,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc+1,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+1,,"black", "MCT225"])
        
        plot(pw.m[cc+2,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc+2,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+2,,"black", "MCT225"])
        
        plot(pw.m[cc+3,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = c(1700, 2048), xlab = "", ylab = "")
        lines(pw.m[cc+3,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+3,,"black", "MCT225"])
    }
    
    # neighbour plots
    {
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], ylim = c(-50000, 3000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
        
   }
}

# column 880:
{
    pixel.image(pw.m[,,"white", "MCT225"], xlim = c(800, 960), ylim = c(1050, 1250))
    cc <- 881; xl <- c(1050, 2048)
    
    # column plots (white)
    {
        plot(pw.m[cc,,"white", "MCT225"], type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+2,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
        
        legend("right", col = c("black", adjustcolor(c("blue", "red", "green"), alpha = 0.4)), lty = 1,
               legend = paste("Column ", cc + c(0:3), sep = ""), bty = "n")
    }
    
    # power plots
    {
        plot(pw.m[cc,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc,,"black", "MCT225"])
        
        plot(pw.m[cc+1,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+1,,"black", "MCT225"])
        
        plot(pw.m[cc+2,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+2,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+2,,"black", "MCT225"])
        
        plot(pw.m[cc+3,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+3,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+3,,"black", "MCT225"])
    }
    
    # neighbour plots
    {
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], ylim = c(-50000, 3000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
    }
}

# column 1690:
{
    pixel.image(pw.m[,,"white", "MCT225"], xlim = c(1600, 1800), ylim = c(1200, 1400))
    cc <- 1690; xl <- c(1050, 2048)
    
    # column plots (white)
    {
        plot(pw.m[cc,,"white", "MCT225"], type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+2,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
        
        legend("right", col = c("black", adjustcolor(c("blue", "red", "green"), alpha = 0.4)), lty = 1,
               legend = paste("Column ", cc + c(0:3), sep = ""), bty = "n")
    }
    
    # power plots
    {
        plot(pw.m[cc,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc,,"black", "MCT225"])
        
        plot(pw.m[cc+1,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+1,,"black", "MCT225"])
        
        plot(pw.m[cc+2,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+2,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+2,,"black", "MCT225"])
        
        plot(pw.m[cc+3,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+3,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+3,,"black", "MCT225"])
    }
    
    # neighbour plots
    {
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], ylim = c(-50000, 3000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
     }
}

pixel.image(pw.m[,,"white", "MCT225"], xlim = c(1800, 2000), ylim = c(1550, 1700))
# column 1926:
{
    pixel.image(pw.m[,,"white", "MCT225"], xlim = c(1600, 1800), ylim = c(1200, 1400))
    cc <- 1925; xl <- c(1400, 2048)
    
    # column plots (white)
    {
        plot(pw.m[cc,,"white", "MCT225"], type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+2,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
        
        legend("right", col = c("black", adjustcolor(c("blue", "red", "green"), alpha = 0.4)), lty = 1,
               legend = paste("Column ", cc + c(0:3), sep = ""), bty = "n")
    }
    
    # power plots
    {
        plot(pw.m[cc,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc,,"black", "MCT225"])
        
        plot(pw.m[cc+1,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+1,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+1,,"black", "MCT225"])
        
        plot(pw.m[cc+2,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+2,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+2,,"black", "MCT225"])
        
        plot(pw.m[cc+3,,"white", "MCT225"], col = "gold", type = "l", ylim = c(0,65535), xlim = xl, xlab = "", ylab = "")
        lines(pw.m[cc+3,,"grey", "MCT225"], col = "cyan3")
        lines(pw.m[cc+3,,"black", "MCT225"])
    }
    
    # neighbour plots
    {
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], ylim = c(-50000, 3000), type = "l", xlab = "", ylab = "")
        lines((pw.m[cc+2,,"white", "MCT225"] - pw.m[cc+1,,"white", "MCT225"]), 
              col = adjustcolor("blue", alpha = 0.4))
        lines(pw.m[cc+3,,"white", "MCT225"] - pw.m[cc+2,,"white", "MCT225"],
              col = adjustcolor("red", alpha = 0.4))
        lines(pw.m[cc+4,,"white", "MCT225"] - pw.m[cc+3,,"white", "MCT225"],
              col = adjustcolor("green3", alpha = 0.4))
        
        plot(pw.m[cc+1,,"white", "MCT225"] - pw.m[cc,,"white", "MCT225"], type = "l", xlim = c(0,1140), ylim = c(-1000,1000))
        mean((pw.m[cc+1,0:1140,"white", "MCT225"] - pw.m[cc+0,0:1140,"white", "MCT225"]), na.rm = T)
    }
}

# plot all dead lines together
{
    plot(pw.m[1691,,"white", "MCT225"], col = adjustcolor("blue", alpha = 0.4), type = "l", xlab = "", ylab = "",
         xlim = c(1000, 2048), ylim = c(4000,10000))
    lines(pw.m[883,,"white", "MCT225"], col = adjustcolor("red", alpha = 0.4))
    lines(pw.m[256,,"white", "MCT225"], col = adjustcolor("green3", alpha = 0.4))
    lines(pw.m[449,,"white", "MCT225"], col = adjustcolor("gold", alpha = 0.4))
    lines(pw.m[1927,,"white", "MCT225"], col = adjustcolor("purple", alpha = 0.4))
    
    
    legend("bottomleft", col = adjustcolor(c("blue", "red", "green3", "gold", "purple"), alpha = 0.7), lty = 1, bty = "n",
           legend = paste("Column ", c(1691, 883, 256, 449, 1927), sep = ""))
}

# pdf plots of line ends
{
    le.plot <- function(col, l.end) {
        pdf(paste0("./Notes/MCT225/fig/line-end-", col, ".pdf"), width = 2, height = 2)
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "MCT225"], xlim = col + c(-10, 10), ylim = l.end + c(-10, 10))
        dev.off()
    }
    
    le.plot(256, 1880)
    le.plot(448, 1137)
    le.plot(882, 1119)
    le.plot(1691, 1359)
    le.plot(1926, 1627)
}

####################################################################################################

# RING ARTEFACTS?!                                                                              ####

focal.plot(pw.m[,,"white", "MCT225"], centre = c(893, 345), title = "white")
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(893, 345), title = "grey")
focal.plot(pw.m[,,"black", "MCT225"], centre = c(893, 345))

focal.plot(pw.m[,,"white", "MCT225"], c(360, 1830), surround = 10, dp = 1)
pixel.image(pw.m[,,"white", "MCT225"], xlim = c(340, 370), ylim = c(1815, 1845))

# plots of ring-shaped clusters
{
    rc.plot <- function(col, row) {
        pdf(paste0("./Notes/MCT225/fig/ring-cluster-", col, ".pdf"), width = 2, height = 2)
        par(mar = c(2,2,1,1))
        pixel.image(pw.m[,,"white", "MCT225"], xlim = col + c(-10, 10), ylim = row + c(-10, 10))
        dev.off()
    }
    
    pixel.image(pw.m[,,"white", "MCT225"], xlim = 360 + c(-10, 10), ylim = 1830 + c(-10, 10))
    
    rc.plot(892, 345)
    rc.plot(356, 1829)
    rc.plot(1208, 1955)
    rc.plot(105, 640)
    rc.plot(236, 271)
}

# shading correction of ring-shaped clusters
{
    sc.MCT225 <- shading.corrected(pw.m[,,,"MCT225"])
    sc.plot <- function(col, row) {
        par(mar = c(2,2,1,1))
        pixel.image(sc.MCT225, xlim = col + c(-10, 10), ylim = row + c(-10, 10))
    }

    sc.plot(892, 345)
    sc.plot(356, 1829)
    sc.plot(1208, 1955)
    sc.plot(105, 640)
    sc.plot(236, 271)
}

# plot in all images for quick comparison
{
    sc <- abind("MCT225" = shading.corrected(pw.m[,,,"MCT225"]),
                "160430" = shading.corrected(pw.m[,,,"160430"]), along = 3)
    compare.plots <- function(centre, acq = "MCT225", ...) {
        par(mfrow = c(2, 2))
        focal.plot(pw.m[,,"black", acq], centre, ...)
        focal.plot(pw.m[,,"grey", acq], centre,  ...)
        focal.plot(pw.m[,,"white", acq], centre,  ...)
        focal.plot(sc[,, acq], centre,  ...)
        par(mfrow = c(1,1))
    }
    
    compare.plots(c(892, 345))
    compare.plots(c(356, 1829))
    compare.plots(c(1208, 1955))
    compare.plots(c(105, 640))
    compare.plots(c(236, 271))
}

# find odd by visual scan of sections of image
{
    for (y in 0:5) {
        for (x in 0:5) {
            pixel.image(pw.m[,,"white", "MCT225"], xlim = c(x, x+1) * 350, ylim = c(y, y+1) * 350)
        }
    }
    
    pixel.image(pw.m[,,"white", "MCT225"], xlim = c(650, 1050), ylim = c(200, 600))
    cc <- as.matrix(read.csv(paste0(fpath, "MCT225-anomalies.csv")))
    
    pdf(paste0(fpath, "MCT225-anomalies.pdf")); {
        par(mar = c(2, 2, 1, 1), mfrow = c(4, 3))
        for (i in 1:nrow(cc)) {
            focal.plot(pw.m[,,"white", "MCT225"], centre = cc[i,])
        }
        dev.off()
    }
    
    write.csv(cc, paste0(fpath, "MCT225-anomalies.csv"), quote = F, row.names = F)
    
}

####################################################################################################

# BAD PIXEL MAPS                                                                                ####

Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", NA, "grey", "green3", "blue", "skyblue")
b.px <- basic.bpx(pw.m[,,,"MCT225"])
{
    plot(b.px[, 1:2], pch = ".", col = Cat.cols[b.px$type])
    draw.panels(col = "cyan3")
    abline(h = 512 * c(1:3), col = "green3")
    abline(h = 512 * c(1:4) - 256, col = "gold")
}

pixel.image(pw.m[,,"white", "MCT225"], xlim = c(225, 275), ylim = c(1800, 1900))
# basic thresholding is picking up horizonal subpanel edges

nl.px <- nl.bpx(pw.m[,,,"MCT225"])
nl.cols <- c("blue", "red", "orange", "cyan3", "skyblue", "green3")

plot(nl.px[, 1:2], pch = ".", col = nl.cols[nl.px$type])
plot(nl.px[, 1:2], col = nl.cols[nl.px$type], 
     xlim = c(200, 300), ylim = c(1860, 1960), pch = 15)
pixel.image(pw.m[,,"white", "MCT225"], xlim = c(200, 300), ylim = c(1860, 1960))
focal.plot(pw.m[,,"black", "MCT225"], c(258, 1916))

apply(pw.m[,,,"MCT225"], 3, function(im) {
    med <- median(im, na.rm = T)
    c(v.dim = med * 0.5, dim = med * 0.75, bright = med + 
          (65535 - med)/4, v.bright = med + (65535 - med)/2)
})

####################################################################################################

# PARAMETRIC DESCRIPTION OF SPOT RESPONSE - CUBIC                                               ####

# transect plots across centre of image

# check if poly(2) gives same model as manually inputting trend

# if it does, fit poly(3) model

# alternatively: fit quadratic trend, then fit second trend to residuals, to model loss of sensitivity
####################################################################################################

# TRANSECT ALONG BRIGHT COLUMNS                                                                 ####

# actual values (429)
pdf(paste0(fpath, "col-429-vals.pdf"), height = 4); {
    par(mar = c(2, 2, 1, 1))
    plot(pw.m[429,,"black", "151015"], type = "l", col = "darkblue", 
         xlim = c(1025, 2048), ylim = c(4500,7000), xlab = "", ylab = "")
    lines(pw.m[429,,"black", "160314"], col = "purple")
    lines(pw.m[429,,"black", "160430"], col = "red")
    lines(pw.m[429,,"black", "160705"], col = "orange")
    legend("bottomright", col = c("darkblue", "purple", "red", "orange"), lty = 1,
        legend = c("151015", "160314", "160430", "160705"), bty = "n")
    dev.off()
}

# difference from neighbouring column in successive acquisitions (429)
pdf(paste0(fpath, "col-429-diffs.pdf"), height = 4); {
    par(mar = c(2, 2, 1, 1))
    plot(pw.m[429,,"black", "151015"] - pw.m[428,,"black", "151015"], type = "l", col = "darkblue", 
         xlim = c(1025, 2048), ylim = c(-500, 2000), xlab = "", ylab = "")
    lines(pw.m[429,,"black", "160314"] - pw.m[428,,"black", "160314"], col = "purple")
    lines(pw.m[429,,"black", "160430"] - pw.m[428,,"black", "160430"], col = "red")
    lines(pw.m[429,,"black", "160705"] - pw.m[428,,"black", "160705"], col = "orange")
    dev.off()
}

# difference at each power setting (429)
pdf(paste0(fpath, "col-429-powers.pdf"), height = 4); {
    par(mar = c(2, 2, 1, 1))
    plot(pw.m[429,,"white", "160705"] - pw.m[428,,"white", "160705"], type = "l", col = "gold", 
         xlim = c(1025, 2048), ylim = c(-500, 2000), xlab = "", ylab = "")
    lines(pw.m[429,,"grey", "160705"] - pw.m[428,,"grey", "160705"], col = "green3")
    lines(pw.m[429,,"black", "160705"] - pw.m[428,,"black", "160705"], col = "black")
    legend("bottomright", col = c("black", "green3", "gold"), lty = 2,
           legend = dimnames(pw.m)[[3]], bty = "n")
    dev.off()
}

