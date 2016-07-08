
library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)
fpath <- "./Notes/MCT225/fig/"

####################################################################################################

# XML PROFILE SUMMARY                                                                           ####

c(apply(pw.m, 3:4, mean, na.rm = T))

df <- summarise.profiles()
zz <- cbind(setNames(df[df$date %in% c("160430", "MCT225"), 
                        c("date", "batch", "kV", "uA", "Power", "ExposureTime",
                          "CameraProperties..attrs.gain", "ImageDetector")],
                     nm = c("Date", "Image", "kV", "uA", "Power", "Exposure", "Gain", "Model")), 
            "Mean value" = round(c(apply(pw.m[,,,c("160430", "MCT225")], 3:4, mean, na.rm = T)), 0),
            "Median value" = round(c(apply(pw.m[,,,c("160430", "MCT225")], 3:4, median, na.rm = T)), 0),
            "St. Dev." = round(c(apply(pw.m[,,,c("160430", "MCT225")], 3:4, sd, na.rm = T)), 0))

write.csv(zz[,-c(1:2)], paste0(fpath, "Image-summary.csv"), quote = F, 
          row.names = apply(zz[,1:2], 1, paste, collapse = " - "))

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

####################################################################################################

# RING ARTEFACTS?!                                                                              ####

focal.plot(pw.m[,,"white", "MCT225"], centre = c(893, 345), title = "white")
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(893, 345), title = "grey")
focal.plot(pw.m[,,"black", "MCT225"], centre = c(893, 345))

focal.plot(pw.m[,,"white", "MCT225"], c(360, 1830), surround = 10, dp = 1)
pixel.image(pw.m[,,"white", "MCT225"], xlim = c(340, 370), ylim = c(1815, 1845))

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

