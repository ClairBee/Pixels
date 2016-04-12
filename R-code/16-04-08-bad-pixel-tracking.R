
# PARAMETRIC MODEL & SD USED TO IDENTIFY BAD PIXELS WITH JOHNSON QUANTILES

library("IO.Pixels")
load.pixel.means(); load.pixel.sds()

fpath <- "./Models/Simple-parametric/"
res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))

fpath <- "./Models/Complex-parametric/"
res <- readRDS(paste0(fpath, "Residuals-panel-int-parametric.rds"))

JF.res <- read.csv(paste0(fpath, "JF-residuals.csv"), as.is = T)
JF.sd <- read.csv(paste0(fpath, "JF-sd.csv"), as.is = T)

bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)

#=====================================================================================
# bad pixel function fails on 150108 data. Hist & plot cutoffs to check
####################################################################################################
# Also need to find way to identify lines of bright pixels, which may not be picked up here

# change ordering on 'dead' & 'hot' s.t. dim/bright takes priority (more stringent condition)?
# use ks.test (Kolmogorov-Smirnov) to test equality of sample & distribution.

#--------------------------------------------------------------------------------------
# model fit on black-141009 looks particularly poor.
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
# FIT MODELS, GET ALL RESIDUALS                                                                 ####
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
# BAD PIXELS - 14.10.09                                                                         ####
    
    # get bad pixel map and values at each power setting
    {
        bp.141009 <- get.bad.px("141009")
        
        # get unique coordinates & classification
        px.141009 <- bp.141009[!duplicated(bp.141009[,1:2]),c(1:2, 4:5)]
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

        saveRDS(px.141009, paste0(fpath, "Bad-px-141009.rds"))
    }
    px.141009 <- readRDS(paste0(fpath, "Bad-px-141009.rds"))
    
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
# BAD PIXELS - 14.11.18                                                                         ####

    # get bad pixel map and values at each power setting
        {
            bp.141118 <- get.bad.px("141118")
            
            # get unique coordinates & classification
            px.141118 <- bp.141118[!duplicated(bp.141118[,1:2]),c(1:2, 4:5)]
            table(px.141118$type, px.141118$sd.type)
            
    # results
            {
                #        static quiet noisy     -
                # dead        5     0     0     0
                # dim         0   219     5 16869
                # bright      0     7  2133 15337
                # hot       138     0     0     0
                # -           0  4935  7878     0
                
                #        static quiet noisy     -
                # dead        5     0     0     0
                # dim         0   108     5  8609
                # bright      0     4  1642 12455
                # hot       138     0     0     0
                # -           0   268  1143     0
            }

            saveRDS(px.141118, paste0(fpath, "Bad-px-141118.rds"))
        }
    px.141118 <- readRDS(paste0(fpath, "Bad-px-141118.rds"))

# compare to previous map
    {
        td <- merge(cbind(px.141009[,1:2], "cat" = apply(px.141009[,3:4], 1, paste, collapse = ",")),
                    cbind(px.141118[,1:2], "cat" = apply(px.141118[,3:4], 1, paste, collapse = ",")),
                    by = c("row", "col"), all = T, suffix = c(".141009", ".141118"))
        td[sapply(td, is.factor)] <- lapply(td[sapply(td, is.factor)], factor, 
                                            levels = c(unique(unlist(lapply(td[sapply(td, is.factor)], levels))),"-,-"))
        td[is.na(td)] <- "-,-"
        
        td$status <- factor("changed", levels = c("new", "fixed", "same", "changed"))
        td$status[td$cat.141009 == td$cat.141118] <- "same"
        td$status[td$cat.141009 == "-,-"] <- "new"
        td$status[td$cat.141118 == "-,-"] <- "fixed"
        
        table(td$status)    #     New   Fixed    Same Changed 
                            #    7951    8655   14844    1577 
        
        chng <- table(td[td$status == "changed",3:4])
        chng <- rbind(chng, colSums(chng))
        chng <- cbind(chng, rowSums(chng))
        
        table(data.frame("cat.141009" = read.table(text = as.character(td[td$status == "changed","cat.141009"]), sep = ",", colClasses = "character")[,1],
                         "cat.141118" = read.table(text = as.character(td[td$status == "changed","cat.141118"]), sep = ",", colClasses = "character")[,1]))
        table(data.frame("cat.141009" = read.table(text = as.character(td[td$status == "changed","cat.141009"]), sep = ",", colClasses = "character")[,2],
                         "cat.141118" = read.table(text = as.character(td[td$status == "changed","cat.141118"]), sep = ",", colClasses = "character")[,2]))
        
        table(td[td$status == "fixed",3:4])
        sum(colSums(table(td[td$status == "changed",3:4]))); sum(rowSums(table(td[td$status == "changed",3:4])))
        
        }
    
# development of bad pixels - by value only

    
####################################################################################################
    
    
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
    