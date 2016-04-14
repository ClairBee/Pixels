
# LOOKING AT ARRANGEMENT OF HOT PIXELS: DO THEY SPREAD VERTICALLY OR HORIZONTALLY?

library("IO.Pixels")

load.pixel.means()

fpath <- "./Models/Simple-parametric/"
res <- readRDS(paste0(fpath, "Residuals-simple-parametric.rds"))
bp <- readRDS(paste0(fpath, "Bad-pixels.rds"))

####################################################################################################
# pair of dead pixels at c. [250, 150]: worth further investigation
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
xy$shape[xy$count <= 9] <- "Cluster"
xy$shape[xy$width == 1] <- "Vert"
xy$shape[xy$height == 1] <- "Horiz"
xy$shape[xy$width == 1 & xy$height == 1] <- "Singleton"

table(xy$type, xy$shape) 

plot(xy[xy$shape == "N/A", c("x", "y")], pch = 20, cex = 0.1, 
     col = c("blue", "green3", "gold", "red")[xy$type[xy$shape == "N/A"]])

table(round(xy[xy$type == "bright" & xy$shape == "N/A", "ratio"],2))

# RELATIONSHIP BETWEEN TYPE & SHAPE                                                             ####

table(xy[xy$type == "bright" & xy$shape == "Horiz", c("count")])
table(xy[xy$type == "bright" & xy$shape == "Cluster", c("count")])

