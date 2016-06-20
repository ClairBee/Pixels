
# identification of 'root' pixel for each feature

library("IO.Pixels"); library("CB.Misc")

bp <- readRDS("./Notes/Final-classifications/fig/bad-px.rds")

####################################################################################################

# get all roots (lines & clusters) together -> single root for 'cluster with line'
# will need to distinguish between line & adjacent bright pixels

# should change this to assign root status to brightest pixel!
bpx.features <- function(bpx) {
    
    # by default, remove edge pixels & screen spots (if still included) - redundant to cluster these
    px <- bpx[!(bpx$type %in% c("edge", "screen.spot")),]
    
    # get clusters without lines to avoid clustering with incidental adjacent pixels
    # (globally bright pixels that lie on lines will be retained)
    px <- px[!(px$type %in% c("line.b", "line.d")),]
    
    
    # ---------------------------------------------------------------------------
    # IDENTIFY CLUSTERS OF DEFECTIVE PIXELS
    
    # convert to raster & get clumps EXCLUDING LINES
    cc <- clump(m2r(bpx2im(px)), dir = 4)
    
    # coordinates of all clumped cells, including unique ID per clump
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    # get size of each clump, discard 1-pixel clumps
    xy <- merge(xy, count(xy, "id"), all = T)[,c("x", "y", "id", "freq")]
    xy <- xy[xy$freq > 1,]
    
    # ---------------------------------------------------------------------------
    # FIND CLUSTER CENTRES
    
    # iterate over all clusters
    xy$f.type <- "cl.body"
    xy$bval <- pw.m[,,"black", "160430"][as.matrix(xy[,1:2])]
    for (i in unique(xy$id)) {
        tmp <- xy[xy$id == i,]
        
        # find brightest pixel in cluster
        cand <- tmp[tmp$bval == max(tmp$bval),]
        if (nrow(cand) > 1) {
            
            # if multiple brightest pixels in cluster, get horizontal midpoint
            cand <- cand[cand$x - floor(mean(tmp$x)) == min(cand$x - floor(mean(tmp$x))),]
            
            if (nrow(cand) > 1) {
                
                # if multiple brightest pixels in midline, get closest to panel edge
                if (min(cand$y) > 992.5) {
                    cand <- cand[which.max(cand$y),]
                } else {
                    cand <- cand[which.min(cand$y),]
                }
            }
        }
        xy$f.type[xy$x == cand$x & xy$y == cand$y] <- "cl.root"
    }
    
    # ---------------------------------------------------------------------------
    # LABEL INDIVIDUAL PIXELS ACCORDING TO FEATURE TYPE
    
    bpx <- merge(bpx, xy[,c("x", "y", "f.type", "id")], by.x = c("row", "col"), by.y = c("x", "y"), all = T)
    
    bpx$f.type[bpx$type == "line.b" & is.na(bpx$f.type)] <- "line.body"
    bpx$f.type[is.na(bpx$f.type)] <- "singleton"
    
    bpx$f.type <- as.factor(bpx$f.type)
    
    
    # ---------------------------------------------------------------------------
    # REMOVE CLUSTERS CONTAINING *ONLY* LOCALLY BRIGHT/DIM PIXELS
    # these are treated as higher density of singleton pixels, not as single 'bleeding' cluster
    
    lc <- ddply(bpx, .(id), summarise, worst = min(type))

    bpx$f.type[bpx$id %in% lc$id[lc$worst %in% c("l.bright", "l.dim")]] <- "singleton"
    
    
    # ---------------------------------------------------------------------------
    # MANUALLY LABEL BRIGHT LINES TO ASSOCIATE WITH CLUSTERS
    # WILL NEED TO AUTOMATE THIS AT SOME POINT!
    
    # label bright column 427 as part of same cluster as the one at its end
    bpx$id[bpx$f.type == "cl.root" & bpx$row == 427]
    bpx$id[bpx$type == "line.b" & bpx$row == 427 & bpx$col >= bpx$col[bpx$f.type == "cl.root" & bpx$row == 427]] <- bpx$id[bpx$f.type == "cl.root" & bpx$row == 427]
    
    # mark line at middle of column 427 (no line root - assume one per column)
    bpx$id[is.na(bpx$id) & bpx$type == "line.b" & bpx$row == "427"] <- max(bpx$id, na.rm = T) + 1

    # bright column 809 has no cluster at end <- assign ID to each part and set cluster root as end point
    bpx$f.type[bpx$row == 809 & bpx$col == 178] <- "cl.root"
    bpx$id[bpx$row == 809 & bpx$col <= 178] <- max(bpx$id, na.rm = T) + 1
    
    bpx$id[bpx$row == 809 & bpx$col > 190 & bpx$type == "line.b"] <- max(bpx$id, na.rm = T) + 1
    
    return(bpx)
}

load.pixel.means()
bp.f <- lapply(bp, bpx.features)
    
#saveRDS(bp.f, "./Notes/Final-classifications/fig/bad-px-by-feature.rds")
saveRDS(bp.f, "./Notes/Final-classifications/fig/bad-px-by-feature-incl-local.rds")

####################################################################################################

# CLUSTERS ONLY CONTAINING LOCALLY BRIGHT PIXELS                                                ####

px <- bp$"160430"

qq <- ddply(px, .(id), summarise, worst = min(type))
table(qq$worst)

px$f.type[px$id %in% qq$id[qq$worst %in% c("l.bright", "l.dim")]] <- "singleton"

####################################################################################################

# FIRST ATTEMPT AT AUTOMATION                                                                   ####

# iterate over all clusters, finding midpoint, then row closest to edge in mid-column

for (i in unique(xy$id)) {
    tmp <- xy[xy$id == i,]
    xt <- floor(mean(xy[xy$id == i,"x"]))
    tmp <- tmp[tmp$x == xt,]
    if (min(tmp$y) > 992.5) {
        yt <- max(tmp$y)
    } else {
        yt <- min(tmp$y)
    }
    xy$f.type[xy$x == xt & xy$y == yt] <- "cl.root"
}

# IDENTIFY CLUSTER CENTRES

# alternative (automated) approaches that didn't work as well
{
    # centre at midpoint & edge closest to panel edge
    # may pick up points not in cluster!
    {
        cl <- ddply(xy, .(id, freq), summarise, 
                    cl.x = floor(mean(x)),
                    cl.y = ((min(y) > 992.5) * max(y)) + ((max(y) < 992.5) * min(y)),
                    f.type = "cl.root")
        
        qq <- merge(xy, cl, by.x = c("x", "y"), by.y = c("cl.x", "cl.y"), all = T)
    }
    
    # need to fix priorities if working this way
    # identify brightest pixel in black images, dimmest in white?
    {
        # get black and white values
        xy$bv <- round(pw.m[,,"black", "160430"][as.matrix(xy[,1:2])],0)
        xy$wv <- round(pw.m[,,"white", "160430"][as.matrix(xy[,1:2])],0)
        
        # need to decide whether to prioritise brightness or non-response?
        # may need some physical explanation for this
        cl <- ddply(xy, .(id), summarise, size = max(freq), 
                    bright.x = x[which.max(bv)], bright.y = y[which.max(bv)],
                    dim.x = x[which.min(wv)], dim.y = y[which.min(wv)],
                    choose = x[which.min(bv[bv == 65535])])
    }
}

# LOOK FOR CLUSTERS ASSOCIATED WITH LINES
# need to automate this!
{
lapply(unique(bpx[bpx$type == "line.b","row"]), 
       function (l) xy[xy$x == l,c("x", "y", "id")])

xy[xy$x %in% bpx[bpx$type == "line.b","row"],]$id

knnx.dist(bpx[bpx$type == "line.b",c("row", "col")], xy[xy$x == 427,c("x", "y")], k = 1)
}

####################################################################################################

ccols <- c("yellow", "skyblue", "cyan3", "green3", "blue", "slateblue1", "purple", "magenta3", "red", "orange")
plot(xy[,c("x", "y")], pch = 15, col = ccols[xy$freq])

table(xy$freq)

xy[xy$freq == 4,]

focal.plot(pw.m[ , , "grey", "160430"], centre = c(255, 323), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)
focal.plot(pw.m[ , , "grey", "160430"], centre = c(1104, 490), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)
focal.plot(pw.m[ , , "grey", "160430"], centre = c(427, 1200), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)
focal.plot(pw.m[ , , "grey", "160430"], centre = c(1186, 1126), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)
focal.plot(pw.m[ , , "grey", "160430"], centre = c(74, 1743), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)
focal.plot(pw.m[ , , "grey", "160430"], centre = c(74, 1743), bad.px = bpx, bpx.cex = 6, pt.cex = 1.2,cex.main = 1)

