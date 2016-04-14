
# PARAMETRIC MODEL & SD USED TO IDENTIFY BAD PIXELS WITH JOHNSON QUANTILES

library("IO.Pixels")
library(data.tree)      # produce tree diagrams plotting progress 

load.pixel.means(); load.pixel.sds()



# simple parametric model: o2 circular spot, linear panels (x + y)
    fpath <- "./Models/Simple-parametric/"
    res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))

JF.res <- read.csv(paste0(fpath, "JF-residuals.csv"), as.is = T)
JF.sd <- read.csv(paste0(fpath, "JF-sd.csv"), as.is = T)
bp <- readRDS(paste0(fpath, "Bad-pixels.rds"))

bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)
bpx <- readRDS(paste0(fpath, "Bad-pixel-df.rds"))


# alternative model for comparison later
{
    fpath <- "./Models/Complex-parametric/"
    res <- readRDS(paste0(fpath, "Residuals-panel-int-parametric.rds"))
}

#=====================================================================================
# bad pixel function fails on 150108 data. Black images are duplicated.
####################################################################################################
# Also need to find way to identify lines of bright pixels, which may not be picked up here

# change ordering on 'dead' & 'hot' s.t. dim/bright takes priority (more stringent condition)?
# use ks.test (Kolmogorov-Smirnov) to test equality of sample & distribution.

#--------------------------------------------------------------------------------------
# model fit on black-141009 looks particularly poor (minipanels)
# SD thresholding seems too strict - picking up a lot of points that aren't really problematic 
#   -> switched to Q0.0001 and Q0.9999 (otherwise v close to central peak)
# hot/dead pixels have 0 SD.

####################################################################################################

# FUNCTIONS                                                                                     ####
get.bad.px <- function(dt) {
    
    # identify bad pixels using Johnson quantiles of residuals (Q.001 and Q.999)
    bp.val <- rbind(data.frame(which(pw.m[,,"black", dt] == 0, arr.ind = T), src = "black", type = "dead"),
                    data.frame(which(pw.m[,,"grey", dt] == 0, arr.ind = T), src = "grey", type = "dead"),
                    data.frame(which(pw.m[,,"white", dt] == 0, arr.ind = T), src = "white", type = "dead"),
                    data.frame(which(pw.m[,,"black", dt] == 65535, arr.ind = T), src = "black", type = "hot"),
                    data.frame(which(pw.m[,,"grey", dt] == 65535, arr.ind = T), src = "grey", type = "hot"),
                    data.frame(which(pw.m[,,"white", dt] == 65535, arr.ind = T), src = "white", type = "hot"),
                    data.frame(which(res[,,"black", dt] > JF.res$q.999[JF.res$batch == "black" & JF.res$img.date == dt], arr.ind = T),
                               src = "black", type = "bright"),
                    data.frame(which(res[,,"grey", dt] > JF.res$q.999[JF.res$batch == "grey" & JF.res$img.date == dt], arr.ind = T),
                               src = "grey", type = "bright"),
                    data.frame(which(res[,,"white", dt] > JF.res$q.999[JF.res$batch == "white" & JF.res$img.date == dt], arr.ind = T),
                               src = "white", type = "bright"),
                    data.frame(which(res[,,"black", dt] < JF.res$q.001[JF.res$batch == "black" & JF.res$img.date == dt], arr.ind = T),
                               src = "black", type = "dim"),
                    data.frame(which(res[,,"grey", dt] < JF.res$q.001[JF.res$batch == "grey" & JF.res$img.date == dt], arr.ind = T),
                               src = "grey", type = "dim"),
                    data.frame(which(res[,,"white", dt] < JF.res$q.001[JF.res$batch == "white" & JF.res$img.date == dt], arr.ind = T),
                               src = "white", type = "dim"))
        
    # identify bad pixels using Johnson quantiles of pixelwise SD (Q.001 and Q.999)
    bp.sd <- rbind(data.frame(which(pw.sd[,,"black", dt] == 0, arr.ind = T), src = "black", sd.type = "static"),
                   data.frame(which(pw.sd[,,"grey", dt] == 0, arr.ind = T), src = "grey", sd.type = "static"),
                   data.frame(which(pw.sd[,,"white", dt] == 0, arr.ind = T), src = "white", sd.type = "static"),
                   data.frame(which(pw.sd[,,"black", dt] < JF.sd$q.0001[JF.sd$batch == "black" & JF.sd$img.date == dt], arr.ind = T),
                              src = "black", sd.type = "quiet"),
                   data.frame(which(pw.sd[,,"grey", dt] < JF.sd$q.0001[JF.sd$batch == "grey" & JF.sd$img.date == dt], arr.ind = T),
                              src = "grey", sd.type = "quiet"),
                   data.frame(which(pw.sd[,,"white", dt] < JF.sd$q.0001[JF.sd$batch == "white" & JF.sd$img.date == dt], arr.ind = T),
                              src = "white", sd.type = "quiet"),
                   data.frame(which(pw.sd[,,"black", dt] > JF.sd$q.9999[JF.sd$batch == "black" & JF.sd$img.date == dt], arr.ind = T),
                              src = "black", sd.type = "noisy"),
                   data.frame(which(pw.sd[,,"grey", dt] > JF.sd$q.9999[JF.sd$batch == "grey" & JF.sd$img.date == dt], arr.ind = T),
                              src = "grey", sd.type = "noisy"),
                   data.frame(which(pw.sd[,,"white", dt] > JF.sd$q.9999[JF.sd$batch == "white" & JF.sd$img.date == dt], arr.ind = T),
                              src = "white", sd.type = "noisy"))
    
    # remove duplicates
    bp.val <- bp.val[!duplicated(bp.val[order(bp.val$type),1:3]),]
    bp.sd <- bp.sd[!duplicated(bp.sd[order(bp.sd$sd.type),1:3]),]
    
    bp <- merge(bp.val, bp.sd, all = T)
    
    bp$type <- ordered(bp$type, levels = c( "dead","dim", "bright", "hot", "-"))
    bp$type[is.na(bp$type)] <- "-"
    
    bp$sd.type <- ordered(bp$sd.type, levels = c("static", "quiet", "noisy",  "-"))
    bp$sd.type[is.na(bp$sd.type)] <- "-"
    
    return(bp[order(bp$type, bp$sd.type),])
}

####################################################################################################
# FIT MODELS, GET ALL RESIDUALS (~20mins)                                                       ####
    # all models fitted are robust
    res <- array(dim = c(1996, 1996, 3, 11), dimnames = dimnames(pw.m))

    # took around 20 mins to run (1042 elapsed)
    pb <- txtProgressBar(min = 0, max = 11, style = 3)
    pb.col <- txtProgressBar(min = 0, max = 3, style = 3)

    for (batch in dimnames(pw.m)[[3]]) {
        for (dt in dimnames(pw.m)[[4]]) {
            
            # fit circular spot (quadratic in z)
            spot <- spot.lm(pw.m[ , , batch, dt], o = 2, robust = T)
            spot.res <- matrix(spot$residuals, ncol = 1996)
            
            panel <- panel.lm(spot.res, "poly(x,2) * poly(y,2)", robust = T)
            res[ , , batch, dt] <- spot.res - panel$fitted.values
            
            setTxtProgressBar(pb, which(dimnames(pw.m)[[4]] == dt))
        }
        setTxtProgressBar(pb.col, which(dimnames(pw.m)[[3]] == batch))
    }

    close(pb); close(pb.col)
    remove(batch, spot, spot.res, panel, pb, pb.col, dt)

    saveRDS(res, paste0(fpath, "Residuals-panel-int-parametric.rds"))

####################################################################################################
# FIT JOHNSON DISTRIBUTIONS                                                                     ####
    
    # fit models to residuals & SD for each acquisition batch
    JF.res <- cbind("batch" = rep(dimnames(res)[[3]], dim(res)[[4]]), 
                    "img.date" = sort(rep(dimnames(res)[[4]], dim(res)[[3]])),
                    do.call("rbind", lapply(apply(res, c(3:4), JohnsonFit, moment = "quant"), data.frame)))
    
    JF.sd <- cbind("batch" = rep(dimnames(pw.sd)[[3]], dim(pw.sd)[[4]]), 
                   "img.date" = sort(rep(dimnames(pw.sd)[[4]], dim(pw.sd)[[3]])),
                   do.call("rbind", lapply(apply(pw.sd, c(3:4), JohnsonFit, moment = "quant"), data.frame)))
    
    # get quantiles & add to table
    JF.res <- cbind(JF.res,
                    do.call("rbind", lapply(split(JF.res[,3:7], seq(nrow(JF.res))),
                                            qJohnson, p = c(0.001, 0.01, 0.99, 0.999))))
    colnames(JF.res) <- c(colnames(JF.sd), c("q.001", "q.01", "q.99", "q.999"))

    JF.sd <- cbind(JF.sd,
                   do.call("rbind", lapply(split(JF.sd[,3:7], seq(nrow(JF.sd))),
                                           qJohnson, p = c(0.0001, 0.001, 0.999, 0.9999))))
    colnames(JF.sd) <- c(colnames(JF.sd[,1:7]), c("q.0001", "q.001", "q.999", "q.9999"))
    
    write.csv(JF.res, paste0(fpath, "JF-residuals.csv"), row.names = F)
    write.csv(JF.sd, paste0(fpath, "JF-sd.csv"), row.names = F)
   
#################################################################################################### 
# GET BAD PIXEL MAPS                                                                            ####
    bp <- list()
    for (d in dimnames(pw.m)[[4]][-4]) {
        bp[[d]] <- get.bad.px(d)
    }
    saveRDS(bp, paste0(fpath, "Bad-pixels.rds"))
    
bp <- readRDS(paste0(fpath, "Bad-pixels.rds"))

# combine all bad pixel maps into one enormo map
bpx <- do.call("rbind", lapply(bp, "[", c(1:3)))
bpx <- bpx[!duplicated(bpx),]

for (df in names(bp)) {
    bpx <- merge(bpx,
                 cbind(bp[[df]][,c("row", "col", "src")],
                       "cat" = apply(bp[[df]][,c("type", "sd.type")], 1, paste, collapse = "-")),
                 all = T, by = c("row", "col", "src"))
    colnames(bpx)[ncol(bpx)] <- paste0("cat.", df)
}

bpx[,4:ncol(bpx)] <- lapply(bpx[,4:ncol(bpx)], factor, 
                            levels = c(unique(unlist(lapply(bpx[,4:ncol(bpx)], levels))),"---"))
bpx[is.na(bpx)] <- "---"

saveRDS(bpx, paste0(fpath, "Bad-pixel-df.rds"))
    
# 14.10.09                                                                                      ####

    # plot histograms
    {
        jpeg(paste0(fpath, "Histograms-141009.jpg"), width = 900, height = 150, pointsize = 12)
        par(mfrow = c(1,6), mar = c(2,2,3,1))
        hist(res[,,"black", "141009"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - black images", ylab = "", xlab = "")
        abline(v = JF.res[JF.res$batch == "black" & JF.res$img.date == 141009, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "black" & JF.res$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)

        hist(res[,,"grey", "141009"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - grey images")
        abline(v = JF.res[JF.res$batch == "grey" & JF.res$img.date == 141009, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "grey" & JF.res$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(res[,,"white", "141009"], breaks = "fd", prob = T, xlim = c(-2000,2000), main = "Residuals - white images")
        abline(v = JF.res[JF.res$batch == "white" & JF.res$img.date == 141009, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "white" & JF.res$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)
        
        
        hist(pw.sd[,,"black", "141009"], breaks = "fd", prob = T, xlim = c(0,50), main = "SD - black images")
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141009, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141009, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"grey", "141009"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - grey images")
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141009, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141009, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"white", "141009"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - white images")
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141009, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141009, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141009, 3:7])), col = "chartreuse3", lwd = 2)
        dev.off()
        
    }

    # get unique coordinates & classification
    bp.141009 <- bp[["141009"]]
    px.141009 <- bp[["141009"]][!duplicated(bp[["141009"]][,1:2]),c(1:2, 4:5)]

    apply(res[,,,"141009"], 3, sd); apply(res[,,,"141009"], 3, mad)
    apply(pw.sd[,,,"141009"], 3, sd); apply(pw.sd[,,,"141009"], 3, mad)
    #           sd(res)     mad(res)        sd(pw.sd)       mad(pw.sd)
    # black       408          90               8               5                     
    # grey        376         147              36              35
    # white       536         305              54              53
    
    # get bad pixel map and values at each power setting
    {
        table(px.141009$type, px.141009$sd.type)

        # results
        {
            # simple parametric model: quadratic spot, linear panels (x + y)
            #        static   quiet   noisy      -
            # dead        5       0       0      0
            # dim         0     166       6  18344
            # bright      0       4    1803  13677
            # hot       131       0       0      0
            # -           0     560    2165      0
            
            #--------------------------------------------------------------------------------------
            
            # complex parametric model: quadratic spot, o2 panels with int (poly(x,2) * poly(y,2))
            #        static   quiet   noisy       -
            # dead        5       0       0       0
            # dim         0     164       2    9588
            # bright      0       4    1782   10637
            # hot       131       0       0       0
            # -           0     564    2199       0
            
            # not much change: small transfer from 'bright'/'dim' to '-'
            # much stricter classification of 'dim' (less points identified)
            
        }

    }

    # are these cutpoints appropriate for SD discrimination? (better than q.001 and q.999)
        hist(pw.sd[,,"black", "141009"], breaks = "fd", prob = T, xlim = c(0, 500))
        lines(c(0:500), dJohnson(c(0:500), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141009,3:7])), col = "red", lwd = 2)
        abline(v = qJohnson(c(0.0001, 0.9999), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141009,3:7])), col = "red", lty = 2)
        
        jpeg(paste0(fpath, "Residuals.jpg"), width = 900, height = 300, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
            pixel.image(res[,,"black", "141009"], title = "Black image")
            pixel.image(res[,,"grey", "141009"], title = "Grey image")
            pixel.image(res[,,"white", "141009"], title = "White image")
        dev.off()
        
        
    # spatial arrangement of bad pixels detected at each power level
    {
        pdf(paste0(fpath, "Bad-px-detected-141009.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
            plot(bp.141009[bp.141009$src == "black",1:2], pch = 20, asp = T,
                 col = c("blue", "green3", "gold", "red", NA)[bp.141009$type[bp.141009$src == "black"]], 
                 main = "Distribution of bad pixels - black - 14-10-09")
            
            plot(bp.141009[bp.141009$src == "grey",1:2], pch = 20, asp = T,
                 col = c("blue", "green3", "gold", "red", NA)[bp.141009$type[bp.141009$src == "grey"]], 
                 main = "Distribution of bad pixels - grey - 14-10-09")
        
            plot(bp.141009[bp.141009$src == "white",1:2], pch = 20, asp = T,
                 col = c("blue", "green3", "gold", "red", NA)[bp.141009$type[bp.141009$src == "white"]], 
                 main = "Distribution of bad pixels - white - 14-10-09")
        dev.off()
        #----------------------------------------------------------------------------
        pdf(paste0(fpath, "Bad-sd-detected-141009.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
            plot(bp.141009[bp.141009$src == "black",1:2], pch = 20, asp = T,
                 col = c("magenta3", "slateblue1", "orange", NA)[bp.141009$sd.type[bp.141009$src == "black"]], 
                 main = "Distribution of bad SDs - black")
        
            plot(bp.141009[bp.141009$src == "grey",1:2], pch = 20, asp = T,
                 col = c("magenta3", "slateblue1", "orange", NA)[bp.141009$sd.type[bp.141009$src == "grey"]], 
                 main = "Distribution of bad SDs - grey")
        
            plot(bp.141009[bp.141009$src == "white",1:2], pch = 20, asp = T,
                 col = c("magenta3", "slateblue1", "orange", NA)[bp.141009$sd.type[bp.141009$src == "white"]], 
                 main = "Distribution of bad SDs - white")
            dev.off()
    }
    
    # look at value & SD of pixels identified
    px <- px.141009[order(px.141009$row, px.141009$col),]
    gp <- sample.healthy(px)
    
    # black images
    {
        plot(pw.m[,,"black","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"black","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"black","141009"][as.matrix(px[,1:2])],
             pw.sd[,,"black","141009"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - black images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "black" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "black" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
    }

    # grey images
    {
        plot(pw.m[,,"grey","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141009"][as.matrix(px[,1:2])],
             pw.sd[,,"grey","141009"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - grey images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "grey" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "grey" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
    }
    
    # white images
    {
        plot(pw.m[,,"white","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"white","141009"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141009"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141009"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"white","141009"][as.matrix(px[,1:2])],
             pw.sd[,,"white","141009"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "white"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - white images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "white" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "white" & JF.sd$img.date == 141009], col = "skyblue", lty = 3)
    }
    
    # should also compare to official map
    
        
####################################################################################################
# 14.11.18                                                                                      ####
    
    # plot histograms
    {
        jpeg(paste0(fpath, "Histograms-141118.jpg"), width = 900, height = 150, pointsize = 12)
        par(mfrow = c(1,6), mar = c(2,2,3,1))
        hist(res[,,"black", "141118"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - black images", ylab = "", xlab = "")
        abline(v = JF.res[JF.res$batch == "black" & JF.res$img.date == 141118, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "black" & JF.res$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(res[,,"grey", "141118"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - grey images")
        abline(v = JF.res[JF.res$batch == "grey" & JF.res$img.date == 141118, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "grey" & JF.res$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(res[,,"white", "141118"], breaks = "fd", prob = T, xlim = c(-2000,2000), main = "Residuals - white images")
        abline(v = JF.res[JF.res$batch == "white" & JF.res$img.date == 141118, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "white" & JF.res$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        
        
        hist(pw.sd[,,"black", "141118"], breaks = "fd", prob = T, xlim = c(0,50), main = "SD - black images")
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141118, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141118, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"grey", "141118"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - grey images")
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141118, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141118, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"white", "141118"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - white images")
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141118, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141118, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141118, 3:7])), col = "chartreuse3", lwd = 2)
        dev.off()
        
    }
    
    # get unique coordinates & classification
    bp.141118 <- bp[["141118"]]
    px.141118 <- bp[["141118"]][!duplicated(bp[["141118"]][,1:2]),c(1:2, 4:5)]
    
    apply(res[,,,"141118"], 3, sd); apply(res[,,,"141118"], 3, mad)
    apply(pw.sd[,,,"141118"], 3, sd); apply(pw.sd[,,,"141118"], 3, mad)
    #           sd(res)     mad(res)        sd(pw.sd)       mad(pw.sd)
    # black       408          90               8               5                     
    # grey        376         147              36              35
    # white       536         305              54              53
    
    # get bad pixel map and values at each power setting
    {
        table(px.141118$type, px.141118$sd.type)
        
        # results
        {
            # simple parametric model: quadratic spot, linear panels (x + y)
            #        static   quiet   noisy      -
            # dead        5       0       0      0
            # dim         0     166       6  18344
            # bright      0       4    1803  13677
            # hot       131       0       0      0
            # -           0     560    2165      0
            
            #--------------------------------------------------------------------------------------
            
            # complex parametric model: quadratic spot, o2 panels with int (poly(x,2) * poly(y,2))
            #        static   quiet   noisy       -
            # dead        5       0       0       0
            # dim         0     164       2    9588
            # bright      0       4    1782   10637
            # hot       131       0       0       0
            # -           0     564    2199       0
            
            # not much change: small transfer from 'bright'/'dim' to '-'
            # much stricter classification of 'dim' (less points identified)
            
        }
        
    }
    
    # are these cutpoints appropriate for SD discrimination? (better than q.001 and q.999)
    hist(pw.sd[,,"black", "141118"], breaks = "fd", prob = T, xlim = c(0, 500))
    lines(c(0:500), dJohnson(c(0:500), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141118,3:7])), col = "red", lwd = 2)
    abline(v = qJohnson(c(0.0001, 0.9999), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141118,3:7])), col = "red", lty = 2)
    
    jpeg(paste0(fpath, "Residuals.jpg"), width = 900, height = 300, pointsize = 12)
    par(mfrow = c(1,3), mar = c(2,2,3,1))
    pixel.image(res[,,"black", "141118"], title = "Black image")
    pixel.image(res[,,"grey", "141118"], title = "Grey image")
    pixel.image(res[,,"white", "141118"], title = "White image")
    dev.off()
    
    
    # spatial arrangement of bad pixels detected at each power level
    {
        pdf(paste0(fpath, "Bad-px-detected-141118.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
        plot(bp.141118[bp.141118$src == "black",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141118$type[bp.141118$src == "black"]], 
             main = "Distribution of bad pixels - black - 14-11-18")
        
        plot(bp.141118[bp.141118$src == "grey",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141118$type[bp.141118$src == "grey"]], 
             main = "Distribution of bad pixels - grey - 14-11-18")
        
        plot(bp.141118[bp.141118$src == "white",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141118$type[bp.141118$src == "white"]], 
             main = "Distribution of bad pixels - white - 14-11-18")
        dev.off()
        #----------------------------------------------------------------------------
        pdf(paste0(fpath, "Bad-sd-detected-141118.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
        plot(bp.141118[bp.141118$src == "black",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141118$sd.type[bp.141118$src == "black"]], 
             main = "Distribution of bad SDs - black")
        
        plot(bp.141118[bp.141118$src == "grey",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141118$sd.type[bp.141118$src == "grey"]], 
             main = "Distribution of bad SDs - grey")
        
        plot(bp.141118[bp.141118$src == "white",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141118$sd.type[bp.141118$src == "white"]], 
             main = "Distribution of bad SDs - white")
        dev.off()
    }
    
    # look at value & SD of pixels identified
    px <- px.141118[order(px.141118$row, px.141118$col),]
    gp <- sample.healthy(px)
    
    # black images
    {
        plot(pw.m[,,"black","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"black","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"black","141118"][as.matrix(px[,1:2])],
             pw.sd[,,"black","141118"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - black images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "black" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "black" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
    }
    
    # grey images
    {
        plot(pw.m[,,"grey","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141118"][as.matrix(px[,1:2])],
             pw.sd[,,"grey","141118"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - grey images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "grey" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "grey" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
    }
    
    # white images
    {
        plot(pw.m[,,"white","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"white","141118"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141118"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141118"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"white","141118"][as.matrix(px[,1:2])],
             pw.sd[,,"white","141118"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "white"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - white images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "white" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "white" & JF.sd$img.date == 141118], col = "skyblue", lty = 3)
    }
    
    # should also compare to official map
    
    
####################################################################################################
# 14.11.18 vs 14.10.09                                                                          ####

    # get unique coordinates & classification
    px.141118 <- bp[["141118"]][!duplicated(bp[["141118"]][,1:2]),c(1:2, 4:5)]
    
    apply(res[,,,"141118"], 3, sd); apply(res[,,,"141118"], 3, mad)
    apply(pw.sd[,,,"141118"], 3, sd); apply(pw.sd[,,,"141118"], 3, mad)
    #           sd(res)     mad(res)        sd(pw.sd)       mad(pw.sd)
    # black       409          89              10               7                     
    # grey        384         171              68              67
    # white       525         305              51              50
    
    table(px.141118$type, px.141118$sd.type)
            
        #        static quiet noisy      -      (o2 circular spot, panels x + y)
        # dead        5     0     0      0
        # dim         0   109     1  16983
        # bright      0     4  1653  15820
        # hot       138     0     0      0
        # -           0   268  1130      0

# state changes at each power setting
    {
        td <- merge(cbind(bp[["141009"]][,c("row", "col", "src")],
                          "cat" = apply(bp[["141009"]][,c("type", "sd.type")], 1, paste, collapse = "-")),
                    cbind(bp[["141118"]][,c("row", "col", "src")],
                          "cat" = apply(bp[["141118"]][,c("type", "sd.type")], 1, paste, collapse = "-")),
                    by = c("row", "col", "src"), all = T, suffix = c(".141009", ".141118"))
        td[,c("cat.141009", "cat.141118")] <- lapply(td[,c("cat.141009", "cat.141118")], factor, 
                                            levels = c(unique(unlist(lapply(td[,c("cat.141009", "cat.141118")], levels))),"---"))
        td[is.na(td)] <- "---"
        
        table(td[td$src == "black", 4:5])
        {
            #               cat.141118
            # cat.141009     bright,- bright,noisy bright,quiet dead,static dim,- dim,noisy dim,quiet hot,static -,noisy -,quiet   -,-
            # bright,-         9709          576            0           0     0         0         0          0      20       0    2759
            # bright,noisy      740          873            0           0     0         0         0         10      24       0     156
            # bright,quiet        0            0            0           0     0         0         0          0       0       0       0
            # dead,static         0            0            0           5     0         0         0          0       0       0       0
            # dim,-               0            1            0           0  5335         0         0          0       0       0     506
            # dim,noisy           0            0            0           0     4         0         0          0       1       0       0
            # dim,quiet           0            0            0           0     0         0         0          0       0       0       0
            # hot,static          0            3            0           0     0         0         0        128       0       0       0
            # -,noisy            70           37            0           0     1         0         0          0     126       0    1207
            # -,quiet             0            0            0           0     0         0         0          0       0       0      53
            # -,-              4154          169            0           0  1959         0         0          0     558       6       0
        }

        table(td[td$src == "grey", 4:5])
        {
            #              cat.141118
            # cat.141009     bright,- bright,noisy bright,quiet dead,static dim,- dim,noisy dim,quiet hot,static -,noisy -,quiet   -,-
            # bright,-         2287            5            0           0     0         0         0          2       0       0    1213
            # bright,noisy       93           46            0           0     0         0         0          4       0       0       4
            # bright,quiet        1            0            0           0     0         0         0          0       0       0       0
            # dead,static         0            0            0           4     0         0         0          0       0       0       0
            # dim,-               0            0            0           0  8529         1         0          0       0       0    6799
            # dim,noisy           0            0            0           0     0         0         0          0       0       0       0
            # dim,quiet           0            0            0           0   114         0         6          0       0       0       0
            # hot,static          1            1            0           0     0         0         0        143       0       0       0
            # -,noisy             3            0            0           0     0         0         0          0       1       0     367
            # -,quiet             0            0            0           0     0         0         0          0       0       0     260
            # -,-              1718            2            0           0  3467         0         0          0      22       2       0
        }
        
        td.w <- table(td[td$src == "white", 4:5])
        {
            #               cat.141118
            # cat.141009     bright,- bright,noisy bright,quiet dead,static dim,- dim,noisy dim,quiet hot,static -,noisy -,quiet  -,-
            # bright,-          363           18            2           0     0         0         0          4       0       0     55
            # bright,noisy       12           17            0           0     0         0         0          0       1       0      1
            # bright,quiet        2            0            1           0     0         0         0          1       0       0      0
            # dead,static         0            0            0           3     0         0         0          0       0       0      0
            # dim,-               0            0            0           0  4242         1        51          0       0       0   1713
            # dim,noisy           0            0            0           0     0         0         0          0       0       0      1
            # dim,quiet           0            0            0           0    52         0        58          0       0       0      0
            # hot,static          5            0            1           0     0         0         0        192       0       0      0
            # -,noisy             1            0            0           0     0         0         0          0       2       0    380
            # -,quiet             0            0            0           0     0         0         0          0       0       0    249
            # -,-                49            0            0           0   891         0         0          1     395     260      0
            
            round(sweep(table(td[td$src == "white", 4:5]),2,colSums(table(td[td$src == "white", 4:5])),`/`),2)
            #cat.141118
            #cat.141009     bright,- bright,noisy bright,quiet dead,static dim,- dim,noisy dim,quiet hot,static -,noisy -,quiet  -,-
            #    bright,-         0.84         0.51         0.50        0.00  0.00      0.00      0.00       0.02    0.00    0.00 0.02
            #bright,noisy     0.03         0.49         0.00        0.00  0.00      0.00      0.00       0.00    0.00    0.00 0.00
            #bright,quiet     0.00         0.00         0.25        0.00  0.00      0.00      0.00       0.01    0.00    0.00 0.00
            #dead,static      0.00         0.00         0.00        1.00  0.00      0.00      0.00       0.00    0.00    0.00 0.00
            #dim,-            0.00         0.00         0.00        0.00  0.82      1.00      0.47       0.00    0.00    0.00 0.71
            #dim,noisy        0.00         0.00         0.00        0.00  0.00      0.00      0.00       0.00    0.00    0.00 0.00
            #dim,quiet        0.00         0.00         0.00        0.00  0.01      0.00      0.53       0.00    0.00    0.00 0.00
            #hot,static       0.01         0.00         0.25        0.00  0.00      0.00      0.00       0.97    0.00    0.00 0.00
            #-,noisy          0.00         0.00         0.00        0.00  0.00      0.00      0.00       0.00    0.01    0.00 0.16
            #-,quiet          0.00         0.00         0.00        0.00  0.00      0.00      0.00       0.00    0.00    0.00 0.10
            #-,-              0.11         0.00         0.00        0.00  0.17      0.00      0.00       0.01    0.99    1.00 0.00
        }

    }
    
# compare all bad px to previous map
    {
        td <- merge(cbind(px.141009[,1:2], "cat" = apply(px.141009[,3:4], 1, paste, collapse = ",")),
                    cbind(px.141118[,1:2], "cat" = apply(px.141118[,3:4], 1, paste, collapse = ",")),
                    by = c("row", "col"), all = T, suffix = c(".141009", ".141118"))
        td[sapply(td, is.factor)] <- lapply(td[sapply(td, is.factor)], factor, 
                                            levels = c(unique(unlist(lapply(td[sapply(td, is.factor)], levels))),"-,-"))
        td[is.na(td)] <- "-,-"
        
        td$status <- factor("changed", levels = c("new", "gone", "same", "changed"))
        td$status[td$cat.141009 == td$cat.141118] <- "same"
        td$status[td$cat.141009 == "-,-"] <- "new"
        td$status[td$cat.141118 == "-,-"] <- "gone"
        
        table(td$status)    #     New     Fixed    Same   Changed 
                            #    10934    11684   23548     1629 
        
        # state changes by value/SD
        table(data.frame("cat.141009" = read.table(text = as.character(td[td$status == "changed","cat.141009"]), sep = ",", colClasses = "character")[,1],
                         "cat.141118" = read.table(text = as.character(td[td$status == "changed","cat.141118"]), sep = ",", colClasses = "character")[,1]))
        table(data.frame("cat.141009" = read.table(text = as.character(td[td$status == "changed","cat.141009"]), sep = ",", colClasses = "character")[,2],
                         "cat.141118" = read.table(text = as.character(td[td$status == "changed","cat.141118"]), sep = ",", colClasses = "character")[,2]))
        
        table(td[td$status == "gone",3:4])
        sum(colSums(table(td[td$status == "changed",3:4]))); sum(rowSums(table(td[td$status == "changed",3:4])))
        
        }
    

    
####################################################################################################
# 14.12.17                                                                                      ####
    
    # plot histograms
    {
        jpeg(paste0(fpath, "Histograms-141217.jpg"), width = 900, height = 150, pointsize = 12)
        par(mfrow = c(1,6), mar = c(2,2,3,1))
        hist(res[,,"black", "141217"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - black images", ylab = "", xlab = "")
        abline(v = JF.res[JF.res$batch == "black" & JF.res$img.date == 141217, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "black" & JF.res$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(res[,,"grey", "141217"], breaks = "fd", prob = T, xlim = c(-1000,1000), main = "Residuals - grey images")
        abline(v = JF.res[JF.res$batch == "grey" & JF.res$img.date == 141217, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "grey" & JF.res$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(res[,,"white", "141217"], breaks = "fd", prob = T, xlim = c(-2000,2000), main = "Residuals - white images")
        abline(v = JF.res[JF.res$batch == "white" & JF.res$img.date == 141217, c(8,11)], col = "red", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.res[JF.res$batch == "white" & JF.res$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        
        
        hist(pw.sd[,,"black", "141217"], breaks = "fd", prob = T, xlim = c(0,50), main = "SD - black images")
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141217, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141217, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"grey", "141217"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - grey images")
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141217, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141217, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "grey" & JF.sd$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        
        hist(pw.sd[,,"white", "141217"], breaks = "fd", prob = T, xlim = c(0,500), main = "SD - white images")
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141217, c(8,11)], col = "red", lty = 2)
        abline(v = JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141217, c(9,10)], col = "blue", lty = 2)
        lines(c(-2000:2000), dJohnson(c(-2000:2000), list(JF.sd[JF.sd$batch == "white" & JF.sd$img.date == 141217, 3:7])), col = "chartreuse3", lwd = 2)
        dev.off()
        
    }
    
    # get unique coordinates & classification
    bp.141217 <- bp[["141217"]]
    px.141217 <- bp[["141217"]][!duplicated(bp[["141217"]][,1:2]),c(1:2, 4:5)]
    
    apply(res[,,,"141217"], 3, sd); apply(res[,,,"141217"], 3, mad)
    apply(pw.sd[,,,"141217"], 3, sd); apply(pw.sd[,,,"141217"], 3, mad)
    #           sd(res)     mad(res)        sd(pw.sd)       mad(pw.sd)
    # black       408          90               8               5                     
    # grey        376         147              36              35
    # white       536         305              54              53
    
    # get bad pixel map and values at each power setting
    {
        table(px.141217$type, px.141217$sd.type)
        
        # results
        {
            # simple parametric model: quadratic spot, linear panels (x + y)
            #        static   quiet   noisy      -
            # dead        5       0       0      0
            # dim         0     166       6  18344
            # bright      0       4    1803  13677
            # hot       131       0       0      0
            # -           0     560    2165      0
            
            #--------------------------------------------------------------------------------------
            
            # complex parametric model: quadratic spot, o2 panels with int (poly(x,2) * poly(y,2))
            #        static   quiet   noisy       -
            # dead        5       0       0       0
            # dim         0     164       2    9588
            # bright      0       4    1782   10637
            # hot       131       0       0       0
            # -           0     564    2199       0
            
            # not much change: small transfer from 'bright'/'dim' to '-'
            # much stricter classification of 'dim' (less points identified)
            
        }
        
    }
    
    # are these cutpoints appropriate for SD discrimination? (better than q.001 and q.999)
    hist(pw.sd[,,"black", "141217"], breaks = "fd", prob = T, xlim = c(0, 500))
    lines(c(0:500), dJohnson(c(0:500), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141217,3:7])), col = "red", lwd = 2)
    abline(v = qJohnson(c(0.0001, 0.9999), list(JF.sd[JF.sd$batch == "black" & JF.sd$img.date == 141217,3:7])), col = "red", lty = 2)
    
    jpeg(paste0(fpath, "Residuals.jpg"), width = 900, height = 300, pointsize = 12)
    par(mfrow = c(1,3), mar = c(2,2,3,1))
    pixel.image(res[,,"black", "141217"], title = "Black image")
    pixel.image(res[,,"grey", "141217"], title = "Grey image")
    pixel.image(res[,,"white", "141217"], title = "White image")
    dev.off()
    
    
    # spatial arrangement of bad pixels detected at each power level
    {
        pdf(paste0(fpath, "Bad-px-detected-141217.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
        plot(bp.141217[bp.141217$src == "black",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141217$type[bp.141217$src == "black"]], 
             main = "Distribution of bad pixels - black - 14-11-18")
        
        plot(bp.141217[bp.141217$src == "grey",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141217$type[bp.141217$src == "grey"]], 
             main = "Distribution of bad pixels - grey - 14-11-18")
        
        plot(bp.141217[bp.141217$src == "white",1:2], pch = 20, asp = T,
             col = c("blue", "green3", "gold", "red", NA)[bp.141217$type[bp.141217$src == "white"]], 
             main = "Distribution of bad pixels - white - 14-11-18")
        dev.off()
        #----------------------------------------------------------------------------
        pdf(paste0(fpath, "Bad-sd-detected-141217.pdf"), width = 12, height = 4, pointsize = 12)
        par(mfrow = c(1,3), mar = c(2,2,3,1))
        plot(bp.141217[bp.141217$src == "black",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141217$sd.type[bp.141217$src == "black"]], 
             main = "Distribution of bad SDs - black")
        
        plot(bp.141217[bp.141217$src == "grey",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141217$sd.type[bp.141217$src == "grey"]], 
             main = "Distribution of bad SDs - grey")
        
        plot(bp.141217[bp.141217$src == "white",1:2], pch = 20, asp = T,
             col = c("magenta3", "slateblue1", "orange", NA)[bp.141217$sd.type[bp.141217$src == "white"]], 
             main = "Distribution of bad SDs - white")
        dev.off()
    }
    
    # look at value & SD of pixels identified
    px <- px.141217[order(px.141217$row, px.141217$col),]
    gp <- sample.healthy(px)
    
    # black images
    {
        plot(pw.m[,,"black","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"black","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - black, 14-10-09"), pch = 20)
        points(pw.m[,,"black","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"black","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - black, 14-10-09"), pch = 20)
        points(pw.sd[,,"black","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"black","141217"][as.matrix(px[,1:2])],
             pw.sd[,,"black","141217"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - black images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "black" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "black" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
    }
    
    # grey images
    {
        plot(pw.m[,,"grey","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - grey, 14-10-09"), pch = 20)
        points(pw.m[,,"grey","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"grey","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - grey, 14-10-09"), pch = 20)
        points(pw.sd[,,"grey","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"grey","141217"][as.matrix(px[,1:2])],
             pw.sd[,,"grey","141217"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "grey"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - grey images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "grey" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "grey" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
    }
    
    # white images
    {
        plot(pw.m[,,"white","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad pixel SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("blue", "green", "gold", "red", NA)[px$type], alpha = 0.2))
        
        plot(pw.m[,,"white","141217"][gp], ylim = c(0,65535), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD values - white, 14-10-09"), pch = 20)
        points(pw.m[,,"white","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.sd[,,"white","141217"][gp], ylim = c(0,4000), col = adjustcolor("lightgrey", alpha = 0.2),
             ylab = "Pixel value", main = paste0("Bad SD SDs - white, 14-10-09"), pch = 20)
        points(pw.sd[,,"white","141217"][as.matrix(px[,1:2])], pch = 20,
               col = adjustcolor(c("magenta3", "slateblue1", "orange", NA)[px$sd.type], alpha = 0.2))
        
        plot(pw.m[,,"white","141217"][as.matrix(px[,1:2])],
             pw.sd[,,"white","141217"][as.matrix(px[,1:2])],
             col = adjustcolor(c("blue", "green", "orange", "red", "white"), alpha = 0.2)[px$type],
             pch = c(20, 17 ,18, 4)[px$sd.type],
             ylab = "Pixelwise SD", xlab = "Pixelwise mean value",
             main = "Pixelwise mean vs SD - white images, 14-10-09")
        abline(h = JF.sd$q.001[JF.sd$batch == "white" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
        abline(h = JF.sd$q.999[JF.sd$batch == "white" & JF.sd$img.date == 141217], col = "skyblue", lty = 3)
    }
    
    # should also compare to official map
    
    
####################################################################################################
    
# TREE DIAGRAM OF PROGRESS OF ALL POINTS (BY COLOUR)                                            ####
    
    tree <- count(bpx[bpx$src == "black",], names(bpx)[4:ncol(bpx)])
    ps <- apply(cbind("bad.px", tree[,2:(ncol(tree)-1)]), 1, paste, collapse = "/")
    tree$pathString <- ps
    
    g <- as.Node(tree)
    print(g, "cat.141009", "cat.141009", "freq")
    

    # need to figure out how to print all values.
    # Also unable to plot the whole thing.
    
    
    
# ASIDE: BLACK IMAGE FROM 150108 HAS STRANGE SD DISTRIBUTION (TRUNCATED AT 0)                   ####
    hist(pw.sd[,,"black", "150108"], breaks = "fd", xlim = c(0,200))
    hist(pw.sd[,,"black", "150828"], breaks = "fd", add = T, col = adjustcolor("magenta3", alpha = 0.2),
        border = adjustcolor("magenta3", alpha = 0.2))
    hist(pw.sd[,,"black", "141009"], breaks = "fd", add = T, col = adjustcolor("blue", alpha = 0.2),
         border = adjustcolor("blue", alpha = 0.2))
    hist(pw.sd[,,"black", "160314"], breaks = "fd", add = T, col = adjustcolor("orange", alpha = 0.2),
         border = adjustcolor("orange", alpha = 0.2))
    hist(pw.sd[,,"black", "150113"], breaks = "fd", add = T, col = adjustcolor("gold", alpha = 0.2),
         border = adjustcolor("gold", alpha = 0.2))
    
    # image set from 15-01-08 has 22 frames. Does importing the other 2 make a difference?
    pixel.image(pw.sd[,,"black", "150108"], break.levels = sd.levels(pw.sd[,,"black", "150828"]))
    pixel.image(pw.sd[,,"black", "150828"])
    draw.panels()
    
####################################################################################################
    
# ASIDE: OUTLIER DETECTION USING Z-SCORES                                                       ####
    # black residuals - looks good
    {
        z.res.black <- (res[,,"black", "141009"] - median(res[,,"black", "141009"])) / sd(res[,,"black", "141009"])
        m.res.black <- (res[,,"black", "141009"] - median(res[,,"black", "141009"])) / mad(res[,,"black", "141009"])
        
        hist(z.res.black, breaks = "fd", xlim = c(-5,5))
        hist(m.res.black, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }

    # grey residuals - picking up slightly too many
    {
        z.res.grey <- (res[,,"grey", "141009"] - median(res[,,"grey", "141009"])) / sd(res[,,"grey", "141009"])
        m.res.grey <- (res[,,"grey", "141009"] - median(res[,,"grey", "141009"])) / mad(res[,,"grey", "141009"])
        
        hist(z.res.grey, breaks = "fd", xlim = c(-5,5))
        hist(m.res.grey, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }
    
    # white residuals - picking up more slightly too many
    {
        z.res.white <- (res[,,"white", "141009"] - median(res[,,"white", "141009"])) / sd(res[,,"white", "141009"])
        m.res.white <- (res[,,"white", "141009"] - median(res[,,"white", "141009"])) / mad(res[,,"white", "141009"])
        
        hist(z.res.white, breaks = "fd", xlim = c(-5,5))
        hist(m.res.white, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }

    # black SDs - looks ok
    {
        z.sd.black <- (pw.sd[,,"black", "141009"] - median(pw.sd[,,"black", "141009"])) / sd(pw.sd[,,"black", "141009"])
        m.sd.black <- (pw.sd[,,"black", "141009"] - median(pw.sd[,,"black", "141009"])) / mad(pw.sd[,,"black", "141009"])
        
        hist(z.sd.black, breaks = "fd", xlim = c(-5,5))
        hist(m.sd.black, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }
    
    # grey SDs - looks ok
    {
        z.sd.grey <- (pw.sd[,,"grey", "141009"] - median(pw.sd[,,"grey", "141009"])) / sd(pw.sd[,,"grey", "141009"])
        m.sd.grey <- (pw.sd[,,"grey", "141009"] - median(pw.sd[,,"grey", "141009"])) / mad(pw.sd[,,"grey", "141009"])
        
        hist(z.sd.grey, breaks = "fd", xlim = c(-5,5))
        hist(m.sd.grey, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }
    
    # white SDs - looks ok
    {
        z.sd.white <- (pw.sd[,,"white", "141009"] - median(pw.sd[,,"white", "141009"])) / sd(pw.sd[,,"white", "141009"])
        m.sd.white <- (pw.sd[,,"white", "141009"] - median(pw.sd[,,"white", "141009"])) / mad(pw.sd[,,"white", "141009"])
        
        hist(z.sd.white, breaks = "fd", xlim = c(-5,5))
        hist(m.sd.white, breaks = "fd", add = T, col= adjustcolor("blue", alpha = 0.2), border = adjustcolor("blue", alpha = 0.2))
        abline(v = c(-3.5, 3.5), col = "red", lty = 2)
    }
    
# try normalising the data with Johnson transformation, then using z-scores
    dat <- pw.m[,,"white", "141009"]
    dat.sd <- pw.sd[,,"white", "141009"]

    jf.dat <- JohnsonFit(dat, moment = "quant") 
    jf.sd <- JohnsonFit(dat.sd, moment = "quant") 
    
    dat.u <- (dat - jf.dat$xi) / jf.dat$lambda
    dat.f <- dat.u / (1 - dat.u)
    dat.z <- jf.dat$gamma + (jf.dat$delta * dat.f)
        
    pixel.image(dat.z) 
    hist(dat.z, breaks = "fd")
    
    z.trans.white <- (dat.z - median(dat.z)) / sd(dat.z)
    m.trans.white <- (dat.z - median(dat.z)) / mad(dat.z)
    
    hist(z.trans.white, breaks = "fd")
    