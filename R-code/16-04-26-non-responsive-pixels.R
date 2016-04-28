
####################################################################################################
# ADD NON-RESPONSIVE PIXELS TO SUPERCLUSTERS

# look at pixels identified as bright within screen spots. What should be priority?

library("IO.Pixels")

fpath <- "./Notes/Revised-thresholds/fig/"
load.pixel.means()
res <- readRDS("./Models/Simple-parametric/Residuals-simple-parametric.rds")

####################################################################################################
# NON-RESPONSIVE PIXELS                                                                         ####

# back-to-back histogram of black & grey pixels, showing limits of 'normal' behaviour
{
    hb <- hist(pw.m[,,"black", "160314"], breaks = "fd", plot = F)
    
    hg <- hist(pw.m[,,"grey", "160314"], breaks =  "fd", plot = F)
    hg.c <- hist(pw.m[11:1985, 11:1985, "grey", "160314"], breaks = "fd", plot = F)
    hb$counts <- -hb$counts
    
    hw <- hist(pw.m[,,"white", "160314"], breaks =  "fd", plot = F)
    hw.c <- hist(pw.m[11:1985, 11:1985, "white", "160314"], breaks = "fd", plot = F)

    JF <- JohnsonFit(pw.m[,, "black", "160314"])
    
    pdf(paste0(fpath, "grey-hist-no-response.pdf"), height = 4, width = 7)
    {
        par(mar = c(2,2,1,1))
        plot(hg, ylim = c(-30, 30), col = "lightseagreen", border = "lightseagreen", xlab = "", ylab = "", main = "")
        plot(hb, add = T)
        rect(qJohnson(0.005, JF), -50, qJohnson(0.995, JF), 50, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
        plot(hg.c, add = T, col = "blue", border = "blue")
        legend("bottomright", pch = 15, bty = "n", cex = 0.8, inset = 0.05,
               col = c("lightseagreen", "blue", "black", adjustcolor("skyblue", alpha = 0.3)),
               legend = c("All grey pixels", "Grey pixels, cropped", "Black pixels", "Normal black response"))
    }
    dev.off()
    
    
    pdf(paste0(fpath, "white-hist-no-response.pdf"), height = 4, width = 7)
    {
        par(mar = c(2,2,1,1))
        plot(hw, ylim = c(-30, 30), col = "lightseagreen", border = "lightseagreen", xlab = "", ylab = "", main = "")
        plot(hb, add = T)
        rect(qJohnson(0.005, JF), -50, qJohnson(0.995, JF), 50, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
        plot(hw.c, add = T, col = "blue", border = "blue")
        legend("bottomright", pch = 15, bty = "n", cex = 0.8, inset = 0.05,
               col = c("lightseagreen", "blue", "black", adjustcolor("skyblue", alpha = 0.3)),
               legend = c("All white pixels", "White pixels, cropped", "Black pixels", "Normal black response"))
    }
    dev.off()
    
    
}

# FIND CONVEX HULL OF SPOTS & 'FILL IN'                                                         ####
xy <- xyFromCell(zz, which(getValues(zz) == 23))

plot(xyFromCell(zz, which(getValues(zz) == 38)), asp = T, pch = 15)
qq <- dilate(xyFromCell(zz, which(getValues(zz) == 38)), shapeKernel(c(5,5), "disc"))

tmp <- dilate(zz, shapeKernel(c(5,5), "disc"))
tmp[xyFromCell(zz, which(getValues(zz) == 38))] <- 1

points(xyFromCell(zz, which(getValues(zz) == 37))[chull(xyFromCell(zz, which(getValues(zz) == 37))),], col = "red", pch = 15)

points(qq, col = "red")
plot(qq)
xyFromCell(zz, which(!is.na(getValues(zz)))
           
           table(getValues(zz))

####################################################################################################
# BAD PIXEL MAP - NO PARAMETRIC SMOOTHING                                                       ####
zz <- screen.spots(160314)
JF.w <- JohnsonFit(pw.m[11:1985, 11:1985, "white", "160314"], moment = "quant")
JF.g <- JohnsonFit(pw.m[11:1985, 11:1985, "grey", "160314"], moment = "quant")

bp <- setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                   rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                          ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "screen spot", "edge", "v.dim", "dim"))),
               c("row", "col", "type"))

bp <- rbind(bp,
            data.frame(no.response(160314), type = "no response"),
            data.frame(which(pw.m[, , "black", "160314"] == 65535, arr.ind = T), type = "hot"),
            data.frame(which(pw.m[, , "white", "160314"] == 0, arr.ind = T), type = "dead"),
            data.frame(which(pw.m[, , "grey", "160314"] > qJohnson(0.999, JF.g), arr.ind = T), type = "bright"),
            data.frame(which(pw.m[, , "grey", "160314"] < qJohnson(0.001, JF.g), arr.ind = T), type = "dim"),
            data.frame(which(pw.m[, , "white", "160314"] > qJohnson(0.999, JF.w), arr.ind = T), type = "bright"),
            data.frame(which(pw.m[, , "white", "160314"] < qJohnson(0.001, JF.w), arr.ind = T), type = "dim"),
            data.frame(which(pw.m[, , "white", "160314"] > max(pretty(qJohnson(0.999, JF.w))), arr.ind = T), type = "v.bright"),
            data.frame(which(pw.m[, , "grey", "160314"] > max(pretty(qJohnson(0.999, JF.g))), arr.ind = T), type = "v.bright"),
            data.frame(which(pw.m[, , "white", "160314"] < min(pretty(qJohnson(0.001, JF.w))), arr.ind = T), type = "v.dim"),
            data.frame(which(pw.m[, , "grey", "160314"] < min(pretty(qJohnson(0.001, JF.g))), arr.ind = T), type = "v.dim"),
            setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                     c("row", "col", "type")))

bp <- bp[order(bp$type),]
bp <- bp[!duplicated(bp[,1:2]),]
table(bp$type)

# exploratory spatial plot
{
    t.cols <- c("blue", "black", "red", "gold", NA, "lightgrey", NA, "cyan3", "green3")
    plot(bp[, 1:2], pch = 20, asp = T, cex = 0.8, xlab = "", ylab = "",
         col = adjustcolor(t.cols, alpha = 0.4)[bp$type])
    plot.new(); legend("top", col = t.cols, legend = levels(bp$type), pch = 20)
}

# BAD PIXEL MAP WITH PARAMETRIC SMOOTHING                                                       ####

zz <- screen.spots(160314, auto.threshold = T)
JF.w <- JohnsonFit(res[11:1985, 11:1985, "white", "160314"], moment = "quant")
JF.g <- JohnsonFit(res[11:1985, 11:1985, "grey", "160314"], moment = "quant")

vb.w <- (max(res[11:1985, 11:1985, "white", "160314"]) / 2)
vb.g <- (max(res[11:1985, 11:1985, "grey", "160314"]) / 2)
vd.w <- (min(res[11:1985, 11:1985, "white", "160314"]) / 2)
vd.g <- (min(res[11:1985, 11:1985, "grey", "160314"]) / 2)

#-----------------------------------------------------------------------------------------
# Johnson 0.001- & 0.999- quantiles
{
    bp <- setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                       rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                              ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "screen spot", "edge", "v.dim", "dim"))),
                   c("row", "col", "type"))
    
    bp <- rbind(bp,
                data.frame(no.response(160314), type = "no response"),
                data.frame(which(pw.m[, , "black", "160314"] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", "160314"] == 0, arr.ind = T), type = "dead"),
                data.frame(which(res[, , "grey", "160314"] > qJohnson(0.999, JF.g), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "grey", "160314"] < qJohnson(0.001, JF.g), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > qJohnson(0.999, JF.w), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "white", "160314"] < qJohnson(0.001, JF.w), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > vb.w, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "grey", "160314"] > vb.g, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "white", "160314"] < vd.w, arr.ind = T), type = "v.dim"),
                data.frame(which(res[, , "grey", "160314"] < vd.g, arr.ind = T), type = "v.dim"),
                setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                         c("row", "col", "type")))
    
    bp <- bp[order(bp$type),]
    bp <- bp[!duplicated(bp[,1:2]),]
    table(bp$type)
    
    # exploratory spatial plot
    {
        t.cols <- c("blue", "black", "red", "gold", NA, "lightgrey", NA, "cyan3", NA)
        plot(bp[, 1:2], pch = 20, asp = T, cex = 0.8, xlab = "", ylab = "",
             col = adjustcolor(t.cols, alpha = 0.4)[bp$type])
        # plot.new(); legend("top", col = t.cols, legend = levels(bp$type), pch = 20)
    }
    
    # histogram of white pixels, showing unclassified values
    {
        excl <- c("no response", "edge", "hot", "dead", "screen spot", "v.bright")
        qq <- as.matrix(bp[bp$type %in% excl,1:2])
        
        cc <- res[,,"white","160314"]
        cc[qq] <- NA
        cc <- cc[!is.na(cc)]
        
        par(mar = c(2,2,1,1))
        hist(res[,,"white", "160314"], breaks = "fd", xlab = "", ylab = "", main = "")
        hist(cc, breaks = "fd", add = T, col = "cyan3", border = "cyan3")
        
        hist(res[,,"white", "160314"], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0,30))
        hist(cc, breaks = "fd", add = T, col = "cyan3", border = "cyan3")
        
        zz <- merge(data.frame(x = c(floor(min(res[,,"white", "160314"])):ceiling(max(res[,,"white", "160314"])))),
                    count(round(cc, 0)), all = T)
        zz[is.na(zz)] <- 0
        abline(v = min(zz$x[which(zz$freq == 0 & zz$x > mean(cc))]), col = "red")
        abline(v = max(zz$x[which(zz$freq == 0 & zz$x < mean(cc))]), col = "red")
        abline(v = qJohnson(c(0.0001, 0.9999), JF.w), col = "gold")
    } # better limit would be the Johnson 0.0001- & 0.9999- quantiles
}

#-----------------------------------------------------------------------------------------
# use Johnson 0.0001- & 0.9999- quantiles 
{
    bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                             rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                    ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "screen spot", "edge", "v.dim", "dim"))),
                         c("row", "col", "type")),
                data.frame(no.response(160314), type = "no response"),
                data.frame(which(pw.m[, , "black", "160314"] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", "160314"] == 0, arr.ind = T), type = "dead"),
                data.frame(which(res[, , "grey", "160314"] > qJohnson(0.9999, JF.g), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "grey", "160314"] < qJohnson(0.0001, JF.g), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > qJohnson(0.9999, JF.w), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "white", "160314"] < qJohnson(0.0001, JF.w), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > vb.w, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "grey", "160314"] > vb.g, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "white", "160314"] < vd.w, arr.ind = T), type = "v.dim"),
                data.frame(which(res[, , "grey", "160314"] < vd.g, arr.ind = T), type = "v.dim"),
                setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                         c("row", "col", "type")))
    
    bp <- bp[order(bp$type),]
    bp <- bp[!duplicated(bp[,1:2]),]
    table(bp$type)
    
    # spatial plot
    {
        t.cols <- c("blue", "black", "red", "gold", NA, "lightgrey", NA, "cyan3", "green3")
        plot(bp[, 1:2], pch = 20, asp = T, cex = 0.8, xlab = "", ylab = "",
             col = adjustcolor(t.cols, alpha = 0.4)[bp$type])
        # plot.new(); legend("top", col = t.cols, legend = levels(bp$type), pch = 20)
    }
}

#-----------------------------------------------------------------------------------------
# cut dim & bright categories at midpoint, compare numbers

bm.w <- qJohnson(0.9999, JF.w) + (vb.w - qJohnson(0.9999, JF.w))/2
bm.g <- qJohnson(0.9999, JF.g) + (vb.g - qJohnson(0.9999, JF.g))/2

dm.w <- qJohnson(0.0001, JF.w) + (vd.w - qJohnson(0.0001, JF.w))/2
dm.g <- qJohnson(0.0001, JF.g) + (vd.g - qJohnson(0.0001, JF.g))/2
{
    bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                             rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                    ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright2", "bright", "screen spot", "edge", "v.dim", "dim2", "dim"))),
                         c("row", "col", "type")),
                data.frame(no.response(160314), type = "no response"),
                data.frame(which(pw.m[, , "black", "160314"] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", "160314"] == 0, arr.ind = T), type = "dead"),
                data.frame(which(res[, , "grey", "160314"] > qJohnson(0.9999, JF.g), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "grey", "160314"] < qJohnson(0.0001, JF.g), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > qJohnson(0.9999, JF.w), arr.ind = T), type = "bright"),
                data.frame(which(res[, , "white", "160314"] < qJohnson(0.0001, JF.w), arr.ind = T), type = "dim"),
                data.frame(which(res[, , "white", "160314"] > vb.w, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "grey", "160314"] > vb.g, arr.ind = T), type = "v.bright"),
                data.frame(which(res[, , "white", "160314"] < vd.w, arr.ind = T), type = "v.dim"),
                data.frame(which(res[, , "grey", "160314"] < vd.g, arr.ind = T), type = "v.dim"),
                data.frame(which(res[, , "white", "160314"] > bm.w, arr.ind = T), type = "bright2"),
                data.frame(which(res[, , "grey", "160314"] > bm.g, arr.ind = T), type = "bright2"),
                data.frame(which(res[, , "white", "160314"] < dm.w, arr.ind = T), type = "dim2"),
                data.frame(which(res[, , "grey", "160314"] < dm.g, arr.ind = T), type = "dim2"),
                setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                         c("row", "col", "type")))
    
    bp <- bp[order(bp$type),]
    bp <- bp[!duplicated(bp[,1:2]),]
    table(bp$type)
    
    # spatial plot
    {
        t.cols <- c("purple", "black", "red", "gold", "orange", "red", "lightgrey", NA, "blue", "cyan3", "green3")
        plot(bp[, 1:2], pch = 20, asp = T, cex = 0.8, xlab = "", ylab = "",
             col = adjustcolor(t.cols, alpha = 0.4)[bp$type])
        # plot.new(); legend("top", col = t.cols, legend = levels(bp$type), pch = 20)
    }
    
    excl <- c("no response", "edge", "hot", "dead", "screen spot")
    
    lim <- bad.px.limits(res[,,"white", "160314"])
    t.cols <- c("purple", "black", "red", "gold", "orange", "red", "lightgrey", NA, "blue", "cyan3", "green3")
    
    # histogram of residuals (white)
    {
        lim <- bad.px.limits(res[,,"white", "160314"])
        # pdf(paste0(fpath, "white-residual-thresholds.pdf"), height = 4, width = 7)
        par(mar = c(2,2,1,1))
        
        hist(res[11:1985,11:1985,"white", "160314"], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0,30))
        
        # shade cutoffs
        rect(floor(min(res[,,"white", "160314"])), 0, lim$dv, 35, col = adjustcolor("blue", alpha = 0.3), border = NA)
        rect(lim$dv, 0, lim$dm, 35, col = adjustcolor("cyan3", alpha = 0.3), border = NA) 
        rect(lim$dm, 0, lim$ds, 35, col = adjustcolor("green3", alpha = 0.3), border = NA)
        
        rect(lim$bs, 0, lim$bm, 35, col = adjustcolor("red", alpha = 0.3), border = NA)
        rect(lim$bm, 0, lim$bv, 35, col = adjustcolor("orange", alpha = 0.3), border = NA)
        rect(lim$bv, 0, ceiling(max(res[,,"white", "160314"])), 35, col = adjustcolor("gold", alpha = 0.3), border = NA)
        
        #qq <- res[,,"white", "160314"][c(1:1996^2)[-matrix(1:1996^2, nrow = 1996)[as.matrix(bp[bp$type %in% excl, 1:2])]]]
        
        #hist(qq, breaks = "fd", add = T, col = "blue", border = "blue")
        abline(v = unlist(lim), col = "red", lty = 2)

        legend("topleft", col = adjustcolor(t.cols[c(9:11, 6:4)], alpha = 0.4), pch = 15, cex = 0.8, inset = 0.05, bg = "white",
               legend = c("Very dim", "Dim", "Slightly dim", "Slightly bright", "Bright", "Very bright"))
            # dev.off()
    }
    
    # histogram of residuals (grey)
    {
        lim <- bad.px.limits(res[,,"grey", "160314"])
        
        # pdf(paste0(fpath, "grey-residual-thresholds.pdf"), height = 4, width = 7)
        par(mar = c(2,2,1,1))
        hist(res[,,"grey", "160314"], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0,30))
        
        # shade cutoffs
        rect(floor(min(res[,,"white", "160314"])), 0, lim$dv, 35, col = adjustcolor("blue", alpha = 0.3), border = NA)
        rect(lim$dv, 0, lim$dm, 35, col = adjustcolor("cyan3", alpha = 0.3), border = NA) 
        rect(lim$dm, 0, lim$ds, 35, col = adjustcolor("green3", alpha = 0.3), border = NA)
        
        rect(lim$bs, 0, lim$bm, 35, col = adjustcolor("red", alpha = 0.3), border = NA)
        rect(lim$bm, 0, lim$bv, 35, col = adjustcolor("orange", alpha = 0.3), border = NA)
        rect(lim$bv, 0, ceiling(max(res[,,"white", "160314"])), 35, col = adjustcolor("gold", alpha = 0.3), border = NA)
        
        #hist(res[,,"grey", "160314"][c(1:1996^2)[-matrix(1:1996^2, nrow = 1996)[as.matrix(bp[bp$type %in% excl, 1:2])]]], 
        #     breaks = "fd", add = T, col = "darkred", border = "darkred")      
        abline(v = c(dm.g, bm.g, qJohnson(c(0.0001, 0.9999), JF.g), vd.g, vb.g), col = "red", lty = 2)
        
        legend("topleft", col = adjustcolor(t.cols[c(9:11, 6:4)], alpha = 0.4), pch = 15, cex = 0.8, inset = 0.05, bg = "white",
               legend = c("Very dim", "Dim", "Slightly dim", "Slightly bright", "Bright", "Very bright"))
        # dev.off()
        
    }
    
    # recalculate boundaries based on raw values instead of residuals
    {
        JF.w <- JohnsonFit(pw.m[11:1985, 11:1985, "white", "160314"], moment = "quant")
        JF.g <- JohnsonFit(pw.m[11:1985, 11:1985, "grey", "160314"], moment = "quant")
        
        vb.w <- mean(pw.m[11:1985, 11:1985, "white", "160314"]) +
            (max(pw.m[11:1985, 11:1985, "white", "160314"]) - 
                 mean(pw.m[11:1985, 11:1985, "white", "160314"])) / 2
        
        vb.g <- mean(pw.m[11:1985, 11:1985, "grey", "160314"]) + 
            (max(pw.m[11:1985, 11:1985, "grey", "160314"]) -
                 mean(pw.m[11:1985, 11:1985, "grey", "160314"])) / 2
        
        vd.w <- mean(pw.m[11:1985, 11:1985, "white", "160314"]) - 
            (mean(pw.m[11:1985, 11:1985, "white", "160314"]) - 
                 min(pw.m[11:1985, 11:1985, "white", "160314"])) / 2
        
        vd.g <- mean(pw.m[11:1985, 11:1985, "grey", "160314"]) - 
            (mean(pw.m[11:1985, 11:1985, "white", "160314"]) - 
                 min(pw.m[11:1985, 11:1985, "grey", "160314"])) / 2
        
        bm.w <- qJohnson(0.9999, JF.w) + (vb.w - qJohnson(0.9999, JF.w))/2
        bm.g <- qJohnson(0.9999, JF.g) + (vb.g - qJohnson(0.9999, JF.g))/2
        
        dm.w <- qJohnson(0.0001, JF.w) + (vd.w - qJohnson(0.0001, JF.w))/2
        dm.g <- qJohnson(0.0001, JF.g) + (vd.g - qJohnson(0.0001, JF.g))/2
    }
    
    # histogram of raw values (white)
    {
        # pdf(paste0(fpath, "white-pwm-thresholds.pdf"), height = 4, width = 7)
        par(mar = c(2,2,1,1))
        
        hist(pw.m[,,"white", "160314"], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0,30))
        
        # shade cutoffs
        rect(floor(min(pw.m[,,"white", "160314"])), 0, vd.w, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "v.dim")], alpha = 0.3), border = NA)
        rect(vd.w, 0, dm.w, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "dim2")], alpha = 0.3), border = NA) 
        rect(dm.w, 0, qJohnson(0.0001, JF.w), 35, col = adjustcolor(t.cols[which(levels(bp$type) == "dim")], alpha = 0.3), border = NA)
        
        rect(qJohnson(0.9999, JF.w), 0, bm.w, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "bright")], alpha = 0.3), border = NA)
        rect(bm.w, 0, vb.w, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "bright2")], alpha = 0.3), border = NA)
        rect(vb.w, 0, ceiling(max(pw.m[,,"white", "160314"])), 35, col = adjustcolor(t.cols[which(levels(bp$type) == "v.bright")], alpha = 0.3), border = NA)
        
        #qq <- pw.m[,,"white", "160314"][c(1:1996^2)[-matrix(1:1996^2, nrow = 1996)[as.matrix(bp[bp$type %in% excl, 1:2])]]]
        
        #hist(qq, breaks = "fd", add = T, col = "blue", border = "blue")
        abline(v = c(dm.w, bm.w, qJohnson(c(0.0001, 0.9999), JF.w), vd.w, vb.w), col = "red", lty = 2)
        
        legend("topleft", col = adjustcolor(t.cols[c(9:11, 6:4)], alpha = 0.4), pch = 15, cex = 0.8, inset = 0.05, bg = "white",
               legend = c("Very dim", "Dim", "Slightly dim", "Slightly bright", "Bright", "Very bright"))
        # dev.off()
    }
    
    # histogram of raw values (grey)
    {
        # pdf(paste0(fpath, "grey-pwm-thresholds.pdf"), height = 4, width = 7)
        par(mar = c(2,2,1,1))
        
        hist(pw.m[,,"grey", "160314"], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0,30))
        
        # shade cutoffs
        rect(floor(min(pw.m[,,"grey", "160314"])), 0, vd.g, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "v.dim")], alpha = 0.3), border = NA)
        rect(vd.g, 0, dm.g, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "dim2")], alpha = 0.3), border = NA) 
        rect(dm.g, 0, qJohnson(0.0001, JF.g), 35, col = adjustcolor(t.cols[which(levels(bp$type) == "dim")], alpha = 0.3), border = NA)
        
        rect(qJohnson(0.9999, JF.g), 0, bm.g, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "bright")], alpha = 0.3), border = NA)
        rect(bm.g, 0, vb.g, 35, col = adjustcolor(t.cols[which(levels(bp$type) == "bright2")], alpha = 0.3), border = NA)
        rect(vb.g, 0, ceiling(max(pw.m[,,"white", "160314"])), 35, col = adjustcolor(t.cols[which(levels(bp$type) == "v.bright")], alpha = 0.3), border = NA)
        
        # qq <- pw.m[,,"white", "160314"][c(1:1996^2)[-matrix(1:1996^2, nrow = 1996)[as.matrix(bp[bp$type %in% excl, 1:2])]]]
        
        #hist(qq, breaks = "fd", add = T, col = "blue", border = "blue")
        abline(v = c(dm.g, bm.g, qJohnson(c(0.0001, 0.9999), JF.g), vd.g, vb.g), col = "red", lty = 2)
        
        legend("topleft", col = adjustcolor(t.cols[c(9:11, 6:4)], alpha = 0.4), pch = 15, cex = 0.8, inset = 0.05, bg = "white",
               legend = c("Very dim", "Dim", "Slightly dim", "Slightly bright", "Bright", "Very bright"))
        # dev.off()
    }
}
####################################################################################################
# FUNCTIONS                                                                                     ####

bad.px.limits <- function(im, cropped = T) {
    
    if (cropped) {im <- im[11:1985, 11:1985]}
    
    JF <- JohnsonFit(im, moment = "quant")

    bright.v <- mean(im) + (max(im) - mean(im)) / 2
    bright <- mean(im) + abs((bright.v - mean(im)) / 2)
    bright.s <- qJohnson(0.9999, JF)
    
    dim.v <- mean(im) - (mean(im) - min(im)) / 2
    dim <- mean(im) - abs((mean(im) - dim.v) / 2)
    dim.s <- qJohnson(0.0001, JF)
    
    
    list(dv = dim.v, dm = dim, ds = dim.s, bs = bright.s, bm = bright, bv = bright.v)
}

plot.limits <- function(im, cropped = T, hist.height = 30, legend = T, ...) {
    
    lim <- bad.px.limits(im)
    
    t.cols <- adjustcolor(c("blue", "cyan3", "green3", "red", "orange", "gold"), alpha = 0.4)
    
    hist(im[11:1985, 11:1985], breaks = "fd", xlab = "", ylab = "", main = "", ylim = c(0, hist.height))
    
    rect(floor(min(im)), 0, lim$dv, hist.height * 2, col = t.cols[1], border = NA)
    rect(lim$dv, 0, lim$dm, hist.height * 2, col = t.cols[2], border = NA)
    rect(lim$dm, 0, lim$ds, hist.height * 2, col = t.cols[3], border = NA)
    
    rect(lim$bs, 0, lim$bm, hist.height * 2, col = t.cols[4], border = NA)
    rect(lim$bm, 0, lim$bv, hist.height * 2, col = t.cols[5], border = NA)
    rect(lim$bv, 0, ceiling(max(im)), hist.height * 2, col = t.cols[6], border = NA)
    
    abline(v = unlist(lim), col = "red", lty = 2)
    
    if (legend) {
        legend("topleft", col = t.cols, pch = 15, cex = 0.8, inset = 0.05, bg = "white",
               legend = c("Very dim", "Dim", "Slightly dim", "Slightly bright", "Bright", "Very bright"))
    }
}

get.dim.bright.px <- function(im, cropped = T) {
    lim <- bad.px.limits(im, cropped)
    
    if (cropped) {im <- im[11:1985, 11:1985]}
    
    bp <- rbind(data.frame(which(im > lim$bv, arr.ind = T), type = "v.bright"),
                data.frame(which(im > lim$bm, arr.ind = T), type = "bright"),
                data.frame(which(im > lim$bs, arr.ind = T), type = "s.bright"),
                data.frame(which(im < lim$dv, arr.ind = T), type = "v.dim"),
                data.frame(which(im < lim$dm, arr.ind = T), type = "dim"),
                data.frame(which(im < lim$ds, arr.ind = T), type = "s.dim"))
    bp <- bp[!duplicated(bp[,1:2]),]
    return(bp)
}

bad.px <- function(dt, threshold.b = res) {
    dt <- toString(dt)
    bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                             rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                    ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "s.bright", "screen spot", "edge", "v.dim", "dim", "s.dim"))),
                         c("row", "col", "type")),
                data.frame(no.response(dt), type = "no response"),
                data.frame(which(pw.m[, , "black", dt] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", dt] == 0, arr.ind = T), type = "dead"),
                get.dim.bright.px(res[, , "grey", dt]),
                get.dim.bright.px(res[, , "white", dt]),
                screen.spots.xt(dt))
    
    bp <- bp[order(bp$type),]
    bp[!duplicated(bp[,1:2]),]
}
####################################################################################################

# BAD PX BY RESIDUALS VS BY RAW VALUES                                                          ####
# quick plots
{
    pdf(paste0(fpath, "grey-pwm-thresholds.pdf"), height = 4, width = 7)
    par(mar = c(2,2,1,1))
    plot.limits(pw.m[,,"grey", "160314"])
    dev.off()
    
    pdf(paste0(fpath, "white-pwm-thresholds.pdf"), height = 4, width = 7)
    par(mar = c(2,2,1,1))
    plot.limits(pw.m[,,"white", "160314"])
    dev.off()
    
    pdf(paste0(fpath, "grey-residual-thresholds.pdf"), height = 4, width = 7)
    par(mar = c(2,2,1,1))
    plot.limits(res[,,"grey", "160314"])
    dev.off()
    
    pdf(paste0(fpath, "white-residual-thresholds.pdf"), height = 4, width = 7)
    par(mar = c(2,2,1,1))
    plot.limits(res[,,"white", "160314"])
    dev.off()
}

bp.g.pwm <- get.dim.bright.px(pw.m[,,"grey", "160314"])
bp.w.pwm <- get.dim.bright.px(pw.m[,,"white", "160314"])

bp.g.res <- get.dim.bright.px(res[,,"grey", "160314"])
bp.w.res <- get.dim.bright.px(res[,,"white", "160314"])

# combine & compare (grey & white combined)
{
    bp.pwm <- rbind(bp.g.pwm, bp.w.pwm)
    bp.pwm <- bp.pwm[order(bp.pwm$type),]
    bp.pwm <- bp.pwm[!duplicated(bp.pwm[,1:2]),]
    
    bp.res <- rbind(bp.g.res, bp.w.res)
    bp.res <- bp.res[order(bp.res$type),]
    bp.res <- bp.res[!duplicated(bp.res[,1:2]),]
    
    bp.compare <- merge(bp.pwm, bp.res,  by = c("row", "col"), all = T, suffix = c(".val", ".res"))
    write.csv(as.matrix(table(bp.compare$type.val, bp.compare$type.res, useNA = "ifany")),
              paste0(fpath, "table-bp.csv"), quote = F)
}

# combine & compare (grey & white separate)
{
    bp.pwm <- merge(bp.g.pwm, bp.w.pwm, by = c("row", "col"), all = T, suffix = c(".grey", ".white"))
    table("grey" = bp.pwm$type.grey, "white" = bp.pwm$type.white, useNA = "ifany")
    
    bp.white <- merge(bp.w.pwm, bp.w.res,  by = c("row", "col"), all = T, suffix = c(".val", ".res"))
    table("value" = bp.white$type.val, "residual" = bp.white$type.res, useNA = "ifany")
    
    write.csv(table(bp.white$type.val, bp.white$type.res, useNA = "ifany"), 
              paste0(fpath, "table-bp-white.csv"), quote = F)
    
    table("val" = bp.white$type.val,"res" =  bp.white$type.res, useNA = "ifany")
    
    bp.grey <- merge(bp.g.pwm, bp.g.res,  by = c("row", "col"), all = T, suffix = c(".val", ".res"))
    table("value" = bp.grey$type.val, "residual" = bp.grey$type.res, useNA = "ifany")
}


####################################################################################################

# BAD PX CLASSIFICAION                                                                          ####

bp.summ <- rbind.fill(lapply(lapply(lapply(lapply(lapply(bp, "[", 3), table), as.matrix),t), as.data.frame))
rownames(bp.summ) <- dimnames(pw.m)[[4]]

write.csv(bp.summ, paste0(fpath, "bp-table-all.csv"), quote = F)

bp.cols <- c("purple", "black", "red", "orange", "gold", NA, "grey", NA, "green3", "green", "lightskyblue")

pdf(paste0(fpath, "plot-bp-1.pdf"))
par(mar = c(2,2,1,1))
plot(bp[[1]][,1:2], pch = 15, col = bp.cols[bp[[1]]$type], asp = T, cex = 0.6, xlab = "", ylab = "")
dev.off()

pdf(paste0(fpath, "plot-bp-2.pdf"))
par(mar = c(2,2,1,1))
plot(bp[[2]][,1:2], pch = 15, col = bp.cols[bp[[2]]$type], asp = T, cex = 0.6, xlab = "", ylab = "")
dev.off()

pdf(paste0(fpath, "plot-bp-10.pdf"))
par(mar = c(2,2,1,1))
plot(bp[[10]][,1:2], pch = 15, col = bp.cols[bp[[10]]$type], asp = T, cex = 0.6, xlab = "", ylab = "")
dev.off()

pdf(paste0(fpath, "plot-bp-11.pdf"))
par(mar = c(2,2,1,1))
plot(bp[[11]][,1:2], pch = 15, col = bp.cols[bp[[11]]$type], asp = T, cex = 0.6, xlab = "", ylab = "")
dev.off()

pdf(paste0(fpath, "legend.pdf"))
plot.new(); legend("center", pch = 15, col = bp.cols, legend = levels(bp[[2]]$type))
dev.off()


####################################################################################################

# DIM PIXEL BEHAVIOUR: ARE THEY DIM OR NON-RESPONSIVE?                                          ####
dim <- as.matrix(bp[[11]][bp[[11]]$type == "v.dim",1:2])
nr <- as.matrix(bp[[11]][bp[[11]]$type == "no response",1:2])

dim.vals <- apply(pw.m, 3:4, "[", dim)
nr.vals <- apply(pw.m, 3:4, "[", nr)

dim.res <- apply(res, 3:4, "[", dim)
nr.res <- apply(res, 3:4, "[", nr)

{
    pdf(paste0(fpath, "v-dim-values.pdf"))
    matplot(t(dim.vals[,1,]), type = "l", ylim = range(pretty(range(dim.vals))), col = "green3",
            xlab = "", ylab = "")
    matplot(t(dim.vals[,2,]), type = "l", add = T, col = "blue")
    matplot(t(dim.vals[,3,]), type = "l", add = T, col = "red")
    dev.off()
}

{
    pdf(paste0(fpath, "nonresponsive-values.pdf"))
    matplot(t(nr.vals[,1,]), type = "l", ylim = range(pretty(range(nr.vals))), col = "green3",
            xlab = "", ylab = "")
    matplot(t(nr.vals[,2,]), type = "l", add = T, col = "blue")
    matplot(t(nr.vals[,3,]), type = "l", add = T, col = "red")
    dev.off()
}

{
    pdf(paste0(fpath, "v-dim-res.pdf"))
    matplot(t(dim.res[,1,]), type = "l", ylim = range(pretty(range(dim.res))), col = "green3",
            xlab = "", ylab = "")
    matplot(t(dim.res[,2,]), type = "l", add = T, col = "blue")
    matplot(t(dim.res[,3,]), type = "l", add = T, col = "red")
    dev.off()
}

{
    pdf(paste0(fpath, "v-nr-res.pdf"))
    matplot(t(nr.res[,1,]), type = "l", ylim = range(pretty(range(nr.res))), col = "green3",
            xlab = "", ylab = "")
    matplot(t(nr.res[,2,]), type = "l", add = T, col = "blue")
    matplot(t(nr.res[,3,]), type = "l", add = T, col = "red")
    dev.off()
}