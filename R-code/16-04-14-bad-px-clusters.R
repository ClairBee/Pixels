
# LOOKING AT ARRANGEMENT OF HOT PIXELS: DO THEY SPREAD VERTICALLY OR HORIZONTALLY?

library("IO.Pixels")

load.pixel.means()

fpath <- "./Models/Simple-parametric/"
res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))
bp <- readRDS(paste0(fpath, "Bad-pixels.rds"))

####################################################################################################
# pair of dead pixels at c. [250, 150]: worth further investigation?
# need to look at effect of queen's vs rook's dist (prefer to use rook's)

# almost all hot pixels are found in superclusters with pixels of another type
# allows classification of 'super-heated' px: hot pixel bleeding into nearby cells
# (find superclusters among hot/bright singles/clusters)

# possible change in shape classification: use morphological closing
#   size <= 9 by rook's distance: cluster
#   run morph. closing over blobs/other, merge resulting clumps?
#   (though less interested in blobs except in terms of excluding them from analysis)

####################################################################################################

# CONVERT BAD PIXEL MAP TO RASTER                                                               ####

# get list of hot pixel coords
hp <- bp$"160314"[bp$"160314"$type !="-", 1:4]
hp <- hp[!duplicated(hp[,1:2]), c(1:2,4)]

# create matrix to populate raster
px.vals <- array(dim = c(1996, 1996))
px.vals[as.matrix(hp[,1:2])] <- hp$type

    # check that image is correctly oriented!
    {
        tmp <- pw.m[,,"black", "150828"]
        tmp[500:550, c(1500:1550, 1800:1850)] <- 1
        pixel.image(tmp)
        
        # tag top-left corner
        px.vals[500:550, c(1500:1550, 1800:1850)] <- 1
        image(px.vals)
        
        # also checked by clumping over test data:
        # x and y coords from xyFromCell give same coords as in original matrix
    }

# convert to raster
r <- raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5)

# plots
{
    image(r, col = c("blue", "green3", "gold", "red"), asp = T)
    plot(r, col = c("blue", "green3", "gold", "red"), asp = T)
}

# CLUMP HOT PIXELS                                                                              ####

# ideally, want to get pixel type & clump size for each bad pixel identified. 
# clump, get clump sizes and types

xy <- data.frame()
for (val in unique(getValues(r)[!is.na(getValues(r))])) {
    
    # mask all but one type of bad pixel
    tmp <- r; 
    tmp[tmp != val] <- NA
    
    # clump remaining bad pixels together
    cc <- clump(tmp, dir = 4)
    
    xy <- rbind(xy,
                merge(data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                           type = val,
                           id = getValues(cc)[!is.na(getValues(cc))]),
                setNames(data.frame(val, freq(cc)), c("type", "id", "count"))))
}
xy$type <- factor(levels(hp$type)[xy$type], levels = levels(hp$type)[levels(hp$type) != "-"])

# add clump height & width to classify lines/blobs/singles
rng <- cbind(aggregate(x ~ type + id, data = xy, max),
             x.min = aggregate(x ~ type + id, data = xy, min)$x,
             y.max = aggregate(y ~ type + id, data = xy, max)$y,
             y.min = aggregate(y ~ type + id, data = xy, min)$y)

rng$width <- rng$x - rng$x.min + 1
rng$height <- rng$y.max - rng$y.min + 1

xy <- merge(xy, rng[,c("type", "id", "width", "height")])
zz <- apply(xy[, c("width", "height")], 1, max)
xy$density <- xy$count / (apply(xy[, c("width", "height")], 1, max)^2)
xy$ratio <- apply(xy[, c("width", "height")], 1, min) / apply(xy[, c("width", "height")], 1, max)

# CLASSIFY SHAPES                                                                               ####    
# classify shapes based on width & height
xy$shape <- factor("Other", levels = c("Singleton", "Horiz", "Vert", "Cluster", "Blob", "Other"))
xy$shape[xy$ratio >= 0.8] <- "Blob"
xy$shape[xy$width == 1] <- "Vert"
xy$shape[xy$height == 1] <- "Horiz"
xy$shape[xy$count <= 9] <- "Cluster"
xy$shape[xy$count == 1] <- "Singleton"

table(xy$type, xy$shape) 

# RELATIONSHIP BETWEEN TYPE & SHAPE                                                             ####

# clusters of each type (as opposed to # pixels in clusters of each type)
table(xy[!duplicated(xy[,c("type", "id")]), c("shape", "type")])

    #            type
    # shape      dead   dim bright   hot
    # Singleton     3  2493  17723   110        # individual bad pixels
    # Horiz         0     1      0     0
    # Vert          0     5      0     0
    # Cluster       1   860   1287     8        # spreading bad pixels
    # Blob          0    30     31     0        # screen spots/damaged regions
    # Other         0    52     32     0        # screen spots/damaged regions

plot(xy[xy$type == "hot", c("x", "y")], pch = 20, cex = 0.3, main = "hot pixels by shape", ylim = c(0,2200),
     col = c("orange", "blue", "red", "green3", "pink", "yellow")[xy$shape[xy$type == "hot"]], asp = T)
legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("orange", "blue", "red", "green3", "pink", "gold"), legend = levels(xy$shape))
table(xy$shape[xy$type == "hot"])

plot(xy[xy$type == "bright", c("x", "y")], pch = 20, cex = 0.3, main = "bright pixels by shape", ylim = c(0,2200),
     col = c("orange", "blue", "red", "green3", "pink", "yellow")[xy$shape[xy$type == "bright"]], asp = T)
legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("orange", "blue", "red", "green3", "pink", "gold"), legend = levels(xy$shape))
table(xy$shape[xy$type == "bright"])

plot(xy[xy$type == "dim", c("x", "y")], pch = 20, cex = 0.3, main = "dim pixels by shape", ylim = c(0,2200),
     col = c("orange", "blue", "red", "green3", "pink", "yellow")[xy$shape[xy$type == "dim"]], asp = T)
legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("orange", "blue", "red", "green3", "pink", "gold"), legend = levels(xy$shape))
table(xy$shape[xy$type == "dim"])

plot(xy[xy$type == "hot", c("x", "y")], pch = 20, cex = 0.3, main = "Hot pixels by shape", ylim = c(0,2200),
     col = c("orange", "blue", "red", "green3", "pink", "yellow")[xy$shape[xy$type == "hot"]], asp = T)
legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("orange", "blue", "red", "green3", "pink", "gold"), legend = levels(xy$shape))
table(xy$shape[xy$type == "hot"])

plot(xy[xy$shape == "Blob", c("x", "y")], pch = 20, cex = 0.3, 
     col = c("blue", "green3", "gold", "red")[xy$type[xy$shape == "Blob"]])

# SMALL CLUSTERS                                                                                ####
# less interested in blobs except in terms of excluding them from analysis
plot(xy[xy$shape == "Cluster", c("x", "y")], asp = T, pch = 20, cex = 0.3, main = "Clusters by type", ylim = c(0,2200), 
     col = c("blue", "green3", "gold", "red")[xy$type[xy$shape == "Cluster"]])
legend("top", pch = 20, col = c("blue", "green3", "gold", "red"), legend  = levels(xy$type), bty = "n", horiz = T, cex = 0.8)

plot(hp[,1:2], pch = 20, asp = T, cex = 0.2, col = hp[,3], xlim = c(570, 600), ylim = c(1700, 1730))
legend("top", pch = 20, col = c("blue", "green3", "gold", "red"), legend  = levels(xy$type), bty = "n", horiz = T, cex = 0.8)

# start by looking at dead & hot clusters
cl <- xy[xy$type %in% c("dead", "hot") & xy$shape == "Cluster",c("x", "y", "type", "id")]
cl <- cl[!duplicated(cl[,c("type", "id")]), c("x", "y")]

for (i in 1:nrow(cl)) {
    plot(hp[,1:2], pch = 20, asp = T, col = c("blue", "green3", "gold", "red")[hp$type], xlim = cl[i,1] + c(-10, 10), ylim = cl[i,2] + c(-10, 10))
    legend("top", pch = 20, col = c("blue", "green3", "gold", "red"), legend  = levels(xy$type), bty = "n", horiz = T, cex = 0.8)
}

# IDENTIFY SUPERCLUSTERS?                                                                       ####
xy2 <- xy
# create raster containing only hot/bright clusters and singletons
sc <- r
sc[cellFromXY(sc, xy = xy[(xy$count > 9), c("x", "y")])] <- NA      # to include dead/dim in supercluster
sc[cellFromXY(sc, xy = xy[(xy$count > 9) | (xy$type %in% c("dead", "dim")), c("x", "y")])] <- NA

# clump these smaller pixels
sc.cc <- clump(sc, dir = 4)
sc.xy <- merge(data.frame(xyFromCell(sc.cc, which(!is.na(getValues(sc.cc)))),
                          sc.id = getValues(sc.cc)[!is.na(getValues(sc.cc))]),
                  setNames(data.frame(freq(sc.cc)), c("sc.id", "sc.count")))

xy <- merge(xy, sc.xy)
xy$shape2 <- factor(xy$shape, levels = c(levels(xy$shape), "Supercluster", "Edge"))
xy$shape2[xy$count < xy$sc.count] <- "Supercluster"
xy$shape2[xy$x < 20 | xy$x > 1976 | xy$y < 20 | xy$y > 1976] <- "Edge"

table(xy$type, xy$shape2) 

# plots
{
    plot(xy[xy$type == "hot", c("x", "y")], pch = 20, cex = 0.3, main = "Hot pixels by shape2", ylim = c(0,2200),
         col = c("green", "blue", "purple", "green3", "pink", "yellow", "red")[xy$shape2[xy$type == "hot"]], asp = T)
    legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("green", "blue", "purple", "green3", "pink", "yellow", "red"), legend = levels(xy$shape2))
    table(xy$shape2[xy$type == "hot"])
    
    plot(xy[xy$type == "bright", c("x", "y")], pch = 20, cex = 0.3, main = "bright pixels by shape2", ylim = c(0,2200),
         col = c("green", "green3", "purple", "blue", "pink", "yellow", "red")[xy$shape2[xy$type == "bright"]], asp = T)
    legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("green", "green3", "purple", "blue", "pink", "yellow", "red"), legend = levels(xy$shape2))
    table(xy$shape2[xy$type == "bright"])
    
    plot(xy[xy$type == "dim", c("x", "y")], pch = 20, cex = 0.3, main = "dim pixels by shape2", ylim = c(0,2200),
         col = c("green", "blue", "purple", "green3", "pink", "yellow", "red")[xy$shape2[xy$type == "dim"]], asp = T)
    legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("green", "blue", "purple", "green3", "pink", "yellow", "red"), legend = levels(xy$shape2))
    table(xy$shape2[xy$type == "dim"])
    
    plot(xy[xy$type == "dead", c("x", "y")], pch = 20, cex = 0.3, main = "dead pixels by shape2", ylim = c(0,2200),
         col = c("green", "blue", "purple", "green3", "pink", "yellow", "red")[xy$shape2[xy$type == "dead"]], asp = T)
    legend("top", horiz = T, bty = "n", cex = 0.7, pch = 20, col = c("green", "blue", "purple", "green3", "pink", "yellow", "red"), legend = levels(xy$shape2))
    table(xy$shape2[xy$type == "dead"])
}

length(unique(xy$sc.id[xy$shape2 == "Supercluster"])) # 186 unique superclusters after edges removed
# 116 superclusters containing only hot & bright pixels (10 of which close to edge)
sc.id <- unique(xy$sc.id[xy$shape2 == "Supercluster"])

pdf(paste0(fpath,"Superclusters.pdf") )
par(mfrow = c(6,6), mar = c(2,2,1,1))
for (id in sc.id) {
    plot(xy[xy$sc.id == id, c("x", "y")],
         xlim = range(xy$x[xy$sc.id == id]) + c(-5,5),
         ylim = range(xy$y[xy$sc.id == id]) + c(-5,5), main = "", xlab = "", ylab = "", 
         pch = 15, col = c("blue", "green3", "gold", "red")[xy$type[xy$sc.id == id]], asp = F)
}
par(mfrow = c(1,1))
dev.off()

# BAD PIXELS IN OFFSET-CORRECTED IMAGES                                                         ####
corr <- (pw.m[,,"grey", "160314"] - pw.m[,,"black", "160314"]) / (pw.m[,,"white", "160314"] - pw.m[,,"black", "160314"]) * 60000
corr[is.na(corr)] <- 0
panel <- panel.lm(corr, "x+y", robust = T)
res <- (corr - panel$fitted.values)

for (i in 1:nrow(cl)) {
    pixel.image(res, xlim = cl[i,1] + c(-10, 10), ylim = cl[i,2] + c(-10, 10))
}

# clusters of hot & dead pixels are still visible in offset-corrected image.

# check if these bad pixels are flagged in official map - T
bpm <- read.csv("./Other-data/BadPixelMap-160314.csv", as.is = T)
merge(bpm, cl, by.x = c("X", "Y"), by.y = c("x", "y"))
