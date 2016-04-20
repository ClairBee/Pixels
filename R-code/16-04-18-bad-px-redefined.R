
# REDEFINING BAD PIXELS

library("IO.Pixels")

load.pixel.means()
ff <- readRDS("./Other-data/Flat-field-corrected.rds")

# use simplest parametric model (o2-spot, x + y panels)
# look at 16-03-14 data to start with (model fit is ok there)
fpath <- "./Models/Simple-parametric/"
res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))
bp <- readRDS(paste0(fpath, "Bad-pixels.rds"))

####################################################################################################
# use black image to identify hot pixels, white to identify bright?

# redefine supercluster: Cluster containing at least one hot pixel + one hot or bright pixel
#   (so may contain 2 hot or 1 hot 1 bright...)
# no particular pattern in development of superclusters
####################################################################################################

# BRIGHT/HOT PIXELS IN WHITE VS GREY IMAGES                                                     ####

bright <- merge(data.frame(bp[bp$src == "white" & bp$type %in% c("bright", "hot"), c("row", "col", "type")],
                           "white" = T),
                data.frame(bp[bp$src == "grey" & bp$type %in% c("bright", "hot"), c("row", "col", "type")],
                           "grey" = T),
                all = T)
bright[is.na(bright)] <- F
table(bright[,3:5])

bright$class <- factor("Bright, G only", levels = c("Hot, W only", "Hot, both", "Bright, W only", "Bright, both", "Bright, G only"))
bright$class[bright$type == "hot" & bright$grey == F] <- "Hot, W only"
bright$class[bright$type == "hot" & bright$grey == T] <- "Hot, both"
bright$class[bright$type == "bright" & bright$grey == F] <- "Bright, W only"
bright$class[bright$type == "bright" & bright$grey == T & bright$white == T] <- "Bright, both"

table(bright$class)
plot(bright[,1:2], pch = 20, asp = T, ylim = c(0,2200),
     col = c("red", "gold", "skyblue", "green", NA)[bright$class])
legend("top", horiz = T, bty = "n" ,pch = 20, col = c("red", "gold", "skyblue", "green", NA),
       legend = levels(bright$class), cex = 0.8)

# still quite a lot of bright pixels.
# what does this do to the superclusters?
# CLUSTERS & SUPERCLUSTERS IN EACH IMAGE SET                                                    ####
xy.w <- cluster.px(bp[bp$src == "white" & bp$type != "-", 1:4])
table(xy.w[,c("type", "shape")])
length(unique(xy.w$sc.id[xy.w$shape == "Supercluster"]))

# plot all hot pixels with surrounding area
{
    pdf(paste0(fpath, "Hot-pixels-160314-white.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (id in unique(xy.w$id[xy.w$type == "hot"])) {
        plot(xy.w[, c("x", "y")],
             xlim = range(xy.w$x[xy.w$id == id]) + c(-5,5),
             ylim = range(xy.w$y[xy.w$id == id]) + c(-5,5),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.w$type])
        points(xy.w[xy.w$shape == "Supercluster", c("x", "y")], pch = 0)
    }
    par(mfrow = c(1,1))
    dev.off()
}

#----------------------------------------------------------------------
xy.g <- cluster.px(bp[bp$src == "grey" & bp$type != "-", 1:4])
table(xy.g[,c("type", "shape")])
length(unique(xy.g$sc.id[xy.g$shape == "Supercluster"]))

# plot all hot pixels with surrounding area
{
    pdf(paste0(fpath, "Hot-pixels-160314-grey.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (id in unique(xy.g$id[xy.g$type == "hot"])) {
        plot(xy.g[, c("x", "y")],
             xlim = range(xy.g$x[xy.g$id == id]) + c(-5,5),
             ylim = range(xy.g$y[xy.g$id == id]) + c(-5,5),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.g$type])
        points(xy.g[xy.g$shape == "Supercluster", c("x", "y")], pch = 0)
    }
    par(mfrow = c(1,1))
    dev.off()
}

#----------------------------------------------------------------------
xy.b <- cluster.px(bp[bp$src == "black" & bp$type != "-", 1:4])
table(xy.b[,c("type", "shape")])
length(unique(xy.b$sc.id[xy.b$shape == "Supercluster"]))

# plot all hot pixels with surrounding area
{
    pdf(paste0(fpath, "Hot-pixels-160314-black.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (id in unique(xy.b$id[xy.b$type == "hot"])) {
        plot(xy.b[, c("x", "y")],
             xlim = range(xy.b$x[xy.b$id == id]) + c(-5,5),
             ylim = range(xy.b$y[xy.b$id == id]) + c(-5,5),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.b$type])
        points(xy.b[xy.b$shape == "Supercluster", c("x", "y")], pch = 0)
    }
    par(mfrow = c(1,1))
    dev.off()
}

#----------------------------------------------------------------------
# plot all hot pixels, with supercluster, in each image

hpx <- bp[bp$type == "hot", 1:2]
hpx <- hpx[!duplicated(hpx),]

{
    pdf(paste0(fpath, "Hot-pixels-160314.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (id in unique(xy.w$id[xy.w$type == "hot"])) {
        plot(xy.w[, c("x", "y")],
             xlim = range(xy.w$x[xy.w$id == id]) + c(-5,5),
             ylim = range(xy.w$y[xy.w$id == id]) + c(-5,5),
             main = "", xlab = "", ylab = "", pch = 20, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.w$type])

        points(xy.g[, c("x", "y")], pch = 1,
               col = c("blue", "green3", "gold", "red")[xy.g$type], lwd = 2, cex = 1.1)

        points(xy.b[, c("x", "y")], pch = 0,
               col = c("blue", "green3", "gold", "red")[xy.b$type], lwd = 2, cex = 1.1)
    }
    par(mfrow = c(1,1))
    dev.off()
}

hot.px <- merge(xy.w[xy.w$type == "hot", 1:4], xy.g[, 1:4],
                by = c("x", "y"), suffix = c(".w", ".g"), all.x = T)
hot.px <- merge(hot.px, xy.b[, 1:4],
                by = c("x", "y"), suffix = c(".", ".b"), all.x = T)

hot.px$b <- hot.px$g <- hot.px$w <- factor("SB", levels = c("SC", "CH", "CB", "SH", "SB"))

    hot.px$w[hot.px$shape.w == "Supercluster"] <- "SC"
    hot.px$b[hot.px$shape.b == "Supercluster"] <- "SC"
    hot.px$g[hot.px$shape.g == "Supercluster"] <- "SC"

    hot.px$w[hot.px$shape.w == "Cluster" & hot.px$type.w == "hot"] <- "SC"
    hot.px$g[hot.px$shape.g == "Cluster" & hot.px$type.g == "hot"] <- "SC"
    hot.px$b[hot.px$shape.b == "Cluster" & hot.px$type.b == "hot"] <- "SC"
    
    hot.px$w[hot.px$shape.w == "Cluster" & hot.px$type.w == "bright"] <- "CB"
    hot.px$g[hot.px$shape.g == "Cluster" & hot.px$type.g == "bright"] <- "CB"
    hot.px$b[hot.px$shape.b == "Cluster" & hot.px$type.b == "bright"] <- "CB"
    
    hot.px$w[hot.px$shape.w == "Singleton" & hot.px$type.w == "hot"] <- "SH"
    hot.px$g[hot.px$shape.g == "Singleton" & hot.px$type.g == "hot"] <- "SH"
    hot.px$b[hot.px$shape.b == "Singleton" & hot.px$type.b == "hot"] <- "SH"
    
count(hot.px, colnames(hot.px)[11:9])

par(mfrow = c(6,6), mar = c(2,2,1,1))

for (i in 1:nrow(hot.px)) {
    plot(xy.b[, c("x", "y")],
         xlim = range(xy.b$x[xy.b$id == id]) + c(-5,5),
         ylim = range(xy.b$y[xy.b$id == id]) + c(-5,5),
         main = "", xlab = "", ylab = "", pch = 15, asp = F,
         col = c("blue", "green3", "gold", "red")[xy.b$type])
    points(xy.b[xy.b$shape == "Supercluster", c("x", "y")], pch = 0)
}


#----------------------------------------------------------------------

corr <- flat.field.corrected(160314)

# plot white superclusters in flat-field corrected image
sc <- xy.w[xy.w$shape == "Supercluster", c(1:2, 7)]
sc.id <- sc[!duplicated(sc$sc.id), ]
{
    pdf(paste0(fpath, "Superclusters-160314-corrected-white.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (i in 1:nrow(sc.id)) {
        pixel.image(corr, xlim = sc.id[i,1] + c(-10, 10), ylim = sc.id[i,2] + c(-10, 10))
        points(sc[sc$sc.id == sc.id[i,3],1:2], pch = 0, cex = 1.2)
    }
    dev.off()
}

# plot grey superclusters in flat-field corrected image
sc <- xy.g[xy.g$shape == "Supercluster", c(1:2, 7)]
sc <- sc[!duplicated(sc$sc.id), 1:2]
{
    pdf(paste0(fpath, "Superclusters-160314-corrected-grey-1.pdf"), width = 7500, height = 7500)
    par(mfrow = c(n,n), mar = c(2,2,1,1))
    for (i in 1:(n^2)) {
        pixel.image(corr, xlim = sc[i,1] + c(-10, 10), ylim = sc[i,2] + c(-10, 10))
    }
    dev.off()
    
    pdf(paste0(fpath, "Superclusters-160314-corrected-grey-2.pdf"), width = 7500, height = 7500)
    par(mfrow = c(n,n), mar = c(2,2,1,1))
    for (i in (n^2 + c(1:(n^2)))) {
        pixel.image(corr, xlim = sc[i,1] + c(-10, 10), ylim = sc[i,2] + c(-10, 10))
    }
    dev.off()
}

# plot black superclusters in flat-field corrected image
sc <- xy.b[xy.b$shape == "Supercluster", c(1:2, 7)]
sc <- sc[!duplicated(sc$sc.id), 1:2]
{
    pdf(paste0(fpath, "Superclusters-160314-corrected-black.pdf"))
    par(mfrow = c(4,4), mar = c(2,2,1,1))
    for (i in 1:nrow(sc)) {
        pixel.image(corr, xlim = sc[i,1] + c(-10, 10), ylim = sc[i,2] + c(-10, 10))
    }
    dev.off()
}
# CLASSIFICATION ACROSS IMAGE SETS                                                              ####
xy.b <- cluster.px(bp[bp$src == "black" & bp$type != "-", 1:4])
xy.g <- cluster.px(bp[bp$src == "grey" & bp$type != "-", 1:4])
xy.w <- cluster.px(bp[bp$src == "white" & bp$type != "-", 1:4])

# plot by value to allow easy comparison/assessment
{
    pdf(paste0(fpath, "Hot-pixel-values-160314.pdf"))
    par(mfrow = c(8,8), mar = c(2,2,1,1))
    for (id in unique(xy.w$id[xy.w$type == "hot"])) {
        xy <- xy.w[xy.w$id == id, c("x", "y")]
        focus <- matrix(c(xy$x[1] + rep(c(-3:3), 7), xy$y[1] + sort(rep(c(-3:3), 7))), ncol = 2)
        focus[focus <= 0] <- 0
        focus[focus >= 1996] <- 1996

        plot(xy.w[, c("x", "y")],
             xlim = xy$x[1] + c(-3,3),
             ylim = xy$y[1] + c(-3,3),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.w$type])
        text(focus, cex = 0.4,
             labels = round(pw.m[,,"white", "160314"][focus]/1000,0))
        points(xy.w[xy.w$shape == "Supercluster", c("x", "y")], pch = 0)
        
        plot(xy.g[, c("x", "y")],
             xlim = xy$x[1] + c(-3,3),
             ylim = xy$y[1] + c(-3,3),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.g$type])
        text(focus, cex = 0.4,
             labels = round(pw.m[,,"grey", "160314"][focus]/1000,0))
        points(xy.g[xy.g$shape == "Supercluster", c("x", "y")], pch = 0)
        
        plot(xy.b[, c("x", "y")],
             xlim = xy$x[1] + c(-3,3),
             ylim = xy$y[1] + c(-3,3),
             main = "", xlab = "", ylab = "", pch = 15, asp = F,
             col = c("blue", "green3", "gold", "red")[xy.b$type])
        text(focus, cex = 0.4,
             labels = round(pw.m[,,"black", "160314"][focus]/1000,0))
        points(xy.b[xy.b$shape == "Supercluster", c("x", "y")], pch = 0)
        
        plot.new()
    }
    dev.off()
}

# numerical assessment of behaviour: are superclusters always superclusters?
sc <- setNames(xy.w[xy.w$type == "hot", c("x", "y", "count", "sc.count")],
               c("x", "y", "w.count", "w.sc.count"))
sc <- merge(sc,
            setNames(xy.g[, c("x", "y", "count", "sc.count")],
                     c("x", "y", "g.count", "g.sc.count")), all.x = T)
sc <- merge(sc,
            setNames(xy.b[, c("x", "y", "count", "sc.count")],
                     c("x", "y", "b.count", "b.sc.count")), all.x = T)

# is black supercluster ever smaller than white/grey?
table("b larger than g?" = sc$b.sc.count >= sc$g.sc.count,
      "g larger than w?" = sc$g.sc.count >= sc$w.sc.count,
      "b larger than w?" = sc$b.sc.count >= sc$w.sc.count, useNA = "ifany")
# all superclusters are same size or larger in black/grey images thatn in white
# some are larger in grey than black

sc[sc$b.sc.count < sc$g.sc.count,]

# find plots in pixel images
which(unique(xy.w$id[xy.w$type == "hot"]) == xy.w$id[xy.w$x == 48 & xy.w$y == 1233]) %/% 16
which(unique(xy.w$id[xy.w$type == "hot"]) == xy.w$id[xy.w$x == 48 & xy.w$y == 1233]) %% 16

# 4 are larger in grey images than in black - adjacent pixels lie v close to threshold
# DEVELOPMENT OF SUPERCLUSTERS                                                                  ####
bp.b <- lapply(bp, subset, src == "black" & type != "-")
xy.b <- lapply(bp.b, cluster.px)

xy <- do.call("rbind", lapply(xy.b, subset, shape == "Supercluster", c("x", "y")))
xy <- xy[!duplicated(xy),]
# 638 pixels have ever appeared in a supercluster

# need to re-cluster these to avoid plotting duplicates: only 139 superclusters ever identified
{
    px.vals <- array(dim = c(1996, 1996))
    px.vals[as.matrix(xy)] <- 1
    
    cc <- clump(raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5), dir = 4)
    sc <- ddply(data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                           id = getValues(cc)[!is.na(getValues(cc))]),
                .(id), summarise, xm = round(mean(x),0), ym = round(mean(y),0))
}

# plot progression of each cluster over successive acquisitions
{
    pdf(paste0(fpath, "Supercluster-progression.pdf"))
    par(mfrow = c(10,10), mar = c(2,2,1,1))
    
    for (i in 1:nrow(sc)) {
        focus <- matrix(c(sc$xm[i] + rep(c(-3:3), 7),
                          sc$ym[i] + sort(rep(c(-3:3), 7))), ncol = 2)
        focus[focus < 0] <- 0; focus[focus > 1996] <- 1996
        for (j in 1:10) {
            plot(bp.b[[j]][, c("row", "col")],
                 xlim = sc$xm[i] + c(-3,3),
                 ylim = sc$ym[i] + c(-3,3),
                 col = c("blue", "green3", "gold", "red")[bp.b[[j]]$type],
                 main = "", xlab = "", ylab = "", pch = 15, asp = F)
            text(focus, cex = 0.3,
                 labels = round(pw.m[,,"black", j][focus]/1000,0))
            points(xy.b[[j]][xy.b[[j]]$shape == "Supercluster", c("x", "y")], pch = 0)
        }
    }
    dev.off()
}
# SUPERCLUSTERS IN FF-CORRECTED IMAGE                                                           ####
# thresholding over FF-corrected image
s.hist(ff[,,"141009"], prob = T)
JF <- JohnsonFit(ff[,,"141009"])
lines(c(19000:20500), dJohnson(c(19000:20500), JF), lwd = 2, col = "orange")
lines(c(19000:20500), dnorm(c(19000:20500), mean(ff[,,"141009"]), mad(ff[,,"141009"])), lwd = 2, col = "red")

bp.ff <- rbind(data.frame(which(ff[,,"141009"] == 0, arr.ind = T), src = "ff", type = "dead"),
               data.frame(which(ff[,,"141009"] == 65535, arr.ind = T), src = "ff", type = "hot"),
               data.frame(which(ff[,,"141009"] > qnorm(0.9999, mean(ff[,,"141009"]), mad(ff[,,"141009"])), arr.ind = T), src = "ff", type = "v.bright"),
               data.frame(which(ff[,,"141009"] > qnorm(0.999, mean(ff[,,"141009"]), mad(ff[,,"141009"])), arr.ind = T), src = "ff", type = "bright"),
               data.frame(which(ff[,,"141009"] < qnorm(0.0001, mean(ff[,,"141009"]), mad(ff[,,"141009"])), arr.ind = T), src = "ff", type = "v.dim"),
               data.frame(which(ff[,,"141009"] < qnorm(0.001, mean(ff[,,"141009"]), mad(ff[,,"141009"])), arr.ind = T), src = "ff", type = "dim"))

bp.ff <- bp.ff[!duplicated(bp.ff[, c("row", "col")]),]
bp.ff$type <- ordered(bp.ff$type, levels = c("dead", "v.dim", "dim", "bright", "v.bright", "hot"))
table(bp.ff$type)

plot(bp.ff[,1:2], pch = 15, asp = T, main = "", xlab = "", ylab = "", cex = 0.4,
     col = c("black", "blue", "cornflowerblue", "gold", "red", "magenta3")[bp.ff$type])
# definite panelwise effect in bright pixels

plot(bp.ff[,1:2], pch = 15, asp = T, main = "", xlab = "", ylab = "", cex = 0.4,
     col = c("black", "blue", "NA", "NA", "red", "magenta3")[bp.ff$type])

spot <- spot.lm(ff[,,"141009"])
spot.res <- matrix(spot$residuals, ncol = 1996)
panels <- panel.lm(spot.res)
max(panels$fitted.values) - min(panels$fitted.values)
pixel.image(panels$fitted.values)
ff.res <- spot.res - panels$fitted.values
mad(ff.res)

s.hist(ff.res, prob = T)
lines(c(-700:700), dJohnson(c(-700:700), JohnsonFit(ff.res)), lwd = 2, col = "orange")
lines(c(-700:700), dnorm(c(-700:700), mean(ff.res), mad(ff.res)), lwd = 2, col = "red")

bp.ff2 <- rbind(data.frame(which(ff[,,"141009"] == 0, arr.ind = T), src = "ff", type = "dead"),
               data.frame(which(ff[,,"141009"] == 65535, arr.ind = T), src = "ff", type = "hot"),
               data.frame(which(ff.res > qnorm(0.9999, mean(ff.res), mad(ff.res)), arr.ind = T), src = "ff", type = "v.bright"),
               data.frame(which(ff.res > qnorm(0.999, mean(ff.res), mad(ff.res)), arr.ind = T), src = "ff", type = "bright"),
               data.frame(which(ff.res < qnorm(0.0001, mean(ff.res), mad(ff.res)), arr.ind = T), src = "ff", type = "v.dim"),
               data.frame(which(ff.res < qnorm(0.001, mean(ff.res), mad(ff.res)), arr.ind = T), src = "ff", type = "dim"))

bp.ff2$type <- ordered(bp.ff2$type, levels = c("dead", "v.dim", "dim", "bright", "v.bright", "hot"))
plot(bp.ff2[,1:2], pch = 15, asp = T, main = "", xlab = "", ylab = "", cex = 0.4,
     col = c("black", "blue", "cornflowerblue", "gold", "red", "magenta3")[bp.ff2$type])
table(bp.ff2$type)

qnorm(0.9999, mean(ff.res), mad(ff.res)); qJohnson(0.9999, JohnsonFit(ff.res))
qnorm(0.0001, mean(ff.res), mad(ff.res)); qJohnson(0.0001, JohnsonFit(ff.res))

# test complete spatial randomness
ff.ppp <-ppp(bp.ff2$row[bp.ff2$type == "v.dim"], bp.ff2$col[bp.ff2$type == "v.dim"], c(1,1996), c(1,1996)) 
ff.env <- envelope(ff.ppp, Kest, nsim = 100)
mad.test(ff.ppp); plot(ff.env, main = "")

{
    pdf(paste0(fpath, "Supercluster-progression-ffc.pdf"))
    par(mfrow = c(10,10), mar = c(2,2,1,1))
    
    for (i in 1:nrow(sc)) {
        focus <- matrix(c(sc$xm[i] + rep(c(-3:3), 7),
                          sc$ym[i] + sort(rep(c(-3:3), 7))), ncol = 2)
        focus[focus < 0] <- 0; focus[focus > 1996] <- 1996
        for (j in 1:10) {
            plot(bp.b[[j]][, c("row", "col")],
                 xlim = sc$xm[i] + c(-3,3),
                 ylim = sc$ym[i] + c(-3,3),
                 col = c("blue", "green3", "gold", "red")[bp.b[[j]]$type],
                 main = "", xlab = "", ylab = "", pch = 15, asp = F)
            text(focus, cex = 0.3,
                 labels = round(ff[,,j][focus]/1000,0))
            points(xy.b[[j]][xy.b[[j]]$shape == "Supercluster", c("x", "y")], pch = 0)
        }
    }
    dev.off()
}
# FIND UNRESPONSIVE PIXELS                                                                      ####

# get threshold for normal behaviour in black images
bn <- qJohnson(0.9999, JohnsonFit(pw.m[,,"black", "141009"]))

# identify pixels within that range in grey/white images
un.g <- which(matrix(findInterval(pw.m[,,"grey", "141009"], bn), ncol = 1996) == 0, arr.ind = T)
un.w <- which(matrix(findInterval(pw.m[,,"white", "141009"], bn), ncol = 1996) == 0, arr.ind = T)

plot(un.g, pch = 20, xlim = c(1,1996), ylim = c(1,1996), asp = T)
points(un.w, pch = 1, lwd = 2, col = "red")

cbind(apply(pw.m[,,,"141009"], 3, "[", un.g),
      ff[,,"141009"][un.g])

th <- c(0, bn, 65535)
table(findInterval(pw.m[,,"grey", "141009"], th))

resp.w <- pw.m[,,"white", "141009"] / pw.m[,,"black", "141009"]
resp.w[pw.m[,,"black", "141009"] == 0] <- 0
pixel.image(resp.w)
# FIND LINES OF WARM PIXELS                                                                     ####
# NEW DEFINITIONS                                                                               ####

# plot grey image with cutoffs for each pixel type

plot(0, type = "n", xlim = c(0,65535), ylim = c(0,10), bty = "n")
{
    # upper limit for UNRESPONSIVE pixels
    rect(0, 0, qJohnson(0.999, JohnsonFit(pw.m[,,"black", "141009"])), 11,
         col = adjustcolor("skyblue", alpha = 0.4), border = NA)
    
    # lower limit for BRIGHT pixels
    rect(qJohnson(0.9999, JohnsonFit(pw.m[,,"grey", "141009"])), 0, 65535, 11,
         col = adjustcolor("gold", alpha = 0.4), border = NA)
    rect(qJohnson(0.99999, JohnsonFit(pw.m[,,"grey", "141009"])), 0, 65535, 11,
         col = adjustcolor("orange", alpha = 0.4), border = NA)
    
    hist(pw.m[,,"grey", "141009"], breaks = "fd", add = T)
    # dead pixels
    lines(c(0,0), c(0,length(which(pw.m[,,"grey", "141009"] == 0))), col = "blue", lwd = 2)
    # hot pixels
    lines(c(65535, 65535), c(0,length(which(pw.m[,,"grey", "141009"] == 65535))), col = "red", lwd = 2)
}


par(mfrow = c(2,1))
{
    # obtain cutoff directly from density (p(x) = 0)
    m <- mean(pw.m[,,"grey", "141009"])
    JF <- JohnsonFit(pw.m[,,"grey", "141009"])
    dj <- 65535 * dJohnson(c(0:65535), JF)
    cols <- c("red", "red", "orange", "orange", "gold", "gold", 
              "greenyellow", "greenyellow", "green3", "green3")
    
    plot(c(0:65535), dj, pch = 20, xlim = c(10000, 25000))
    for (i in 1:10) {
        abline(v = max(which(round(dj, i) == 0 & c(0:65535) < m)), col = cols[i])
        abline(v = min(which(round(dj, i) == 0 & c(0:65535) > m)), col = cols[i])
    }

    # obtain cutoff from gradient of density
    zz <- data.frame(x = c(0:65535), 
                     y1 = dJohnson(c(0:65535), JF),
                     y2 = dJohnson(c(1:65536), JF))
    zz$grad <- zz$y2 - zz$y1
    plot(zz$x, zz$grad, pch = 20, xlim = c(10000, 25000)); abline(h = 0)

    for (i in 1:10) {
        abline(v = max(zz$x[which((round(zz$grad,i) == 0) & (zz$x < m - 1000))]), col = cols[i])
        abline(v = min(zz$x[which((round(zz$grad,i) == 0) & (zz$x > m + 1000))]), col = cols[i])
    }
}

# add histogram with various potential thresholding methods





####################################################################################################
# STATE SPACE DIAGRAM                                                                           ####
# ADJUST THRESHOLDS TO STABILISE STATE SPACE                                                    ####