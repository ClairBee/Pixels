
# identification of 'root' pixel for each feature

library("IO.Pixels"); library("CB.Misc")

bp <- readRDS("./Notes/Final-classifications/fig/bad-px.rds")

####################################################################################################

# get all roots (lines & clusters) together -> single root for 'cluster with line'
# will need to distinguish between line & adjacent bright pixels
# bpx <- bp$"160430"

# should change this to assign root status to brightest pixel!
bpx.features <- function(bpx, im.dim = c(1996, 1996)) {
    
    # by default, remove edge pixels & screen spots (if still included) - redundant to cluster these
    px <- bpx[!(bpx$type %in% c("edge", "screen.spot")),]
    
    # get clusters without lines to avoid clustering with incidental adjacent pixels
    # (globally bright pixels that lie on lines will be retained)
    px <- px[!(px$type %in% c("line.b", "line.d")),]
    
    
    # ---------------------------------------------------------------------------
    # IDENTIFY CLUSTERS OF DEFECTIVE PIXELS
    
    # convert to raster & get clumps EXCLUDING LINES
    cc <- clump(m2r(bpx2im(px, im.dim = im.dim)), dir = 4)
    
    # coordinates of all clumped cells, including unique ID per clump
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    # get size of each clump, discard 1-pixel clumps
    xy <- merge(xy, count(xy, "id"), all = T)[,c("x", "y", "id", "freq")]
    xy <- xy[xy$freq > 1,]
    xy$f.type <- "cl.body"
    
    # ---------------------------------------------------------------------------
    # FIND CLUSTER CENTRES
    
    bpx <- merge(bpx, xy[,c("x", "y", "f.type", "id")], by.x = c("row", "col"), by.y = c("x", "y"), all = T)
    bpx$type <- ordered(bpx$type, levels = levels(bpx$type))
    
    for (i in unique(bpx$id[!is.na(bpx$id)])) {
        tmp <- bpx[which(bpx$id == i),]
        
        # if most severe defect is local brightness/dimness, set cluster as singletons
        if (min(tmp$type) %in% c("l.bright", "l.dim")) {
            bpx$f.type[bpx$id == i] <- "singleton"
        } else {
            
            # find brightest pixel in cluster
            cand <- tmp[tmp$type == min(tmp$type),]
            
            if (nrow(cand) > 1) {
                
                # if multiple brightest pixels in cluster, get horizontal midpoint
                cand <- cand[cand$row - floor(mean(tmp$row)) == min(cand$row - floor(mean(tmp$row))),]
                
                if (nrow(cand) > 1) {
                    
                    # if multiple brightest pixels in midline, get closest to panel edge
                    if (min(cand$col) > 992.5) {
                        cand <- cand[which.max(cand$col),]
                    } else {
                        cand <- cand[which.min(cand$col),]
                    }
                }
            }
            bpx$f.type[bpx$row == cand$row & bpx$col == cand$col] <- "cl.root"
        }
    }
    
    # ---------------------------------------------------------------------------
    # LABEL REMAINING FEATURE TYPES
    
    bpx$f.type[bpx$type == "line.b" & is.na(bpx$f.type)] <- "line.body"
    bpx$f.type[is.na(bpx$f.type)] <- "singleton"
    
    bpx$f.type <- as.factor(bpx$f.type)
    
    
    # ---------------------------------------------------------------------------
    # MANUALLY LABEL BRIGHT LINES TO ASSOCIATE WITH CLUSTERS
    # WILL NEED TO AUTOMATE THIS AT SOME POINT!
    
    # label bright column 427 as part of same cluster as the one at its end
    bpx$id[bpx$type == "line.b" & bpx$row == 427 & bpx$col >= bpx$col[bpx$f.type == "cl.root" & bpx$row == 427]] <- bpx$id[bpx$f.type == "cl.root" & bpx$row == 427]
    
    # mark line at middle of column 427 (no line root - assume one per column)
    bpx$id[is.na(bpx$id) & bpx$type == "line.b" & bpx$row == "427"] <- max(bpx$id, na.rm = T) + 1
    
    # bright column 809 has no cluster at end <- assign ID to each part and set cluster root as end point
    bpx$id[bpx$type == "line.b" & bpx$row == 809 & bpx$col >= bpx$col[bpx$f.type == "cl.root" & bpx$row == 809]] <- bpx$id[bpx$f.type == "cl.root" & bpx$row == 809]
    
    bpx$id[is.na(bpx$id) & bpx$type == "line.b" & bpx$row == "809"] <- max(bpx$id, na.rm = T) + 1
    
    # ---------------------------------------------------------------------------
    
    return(bpx)
}

bp.f <- lapply(bp, bpx.features, im.dim = c(2000,2000))
    
#saveRDS(bp.f, "./Notes/Final-classifications/fig/bad-px-by-feature.rds")
saveRDS(bp.f, "./Other-data/Old-data/bad-px-feature-maps.rds")

####################################################################################################

# APPLIED TO OLD DATA                                                                           ####

# single image for now
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")
bpx <- bp$"131122"

# by default, remove edge pixels & screen spots (if still included) - redundant to cluster these
px <- bpx[!(bpx$type %in% c("edge", "screen.spot", "s.bright", "s.dim", "l.bright", "l.dim")),]

#--------------------------------------------------------
# correct line classifications
ll <- clump(m2r(bpx2im(px, im.dim = c(2000,2000))), dir = 4)

xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                 id = getValues(cc)[!is.na(getValues(cc))])
xy$f.type <- "cl.body"

ss <- ddply(xy, .(id, col = x), summarise,
            ymin = min(y), ymax = max(y), length = length(y))
xy$f.type[xy$id %in% ss$id[ss$length > 20]] <- "line.body"

bpx <- merge(bpx, xy, by.x = c("row", "col"), by.y = c("x", "y"), all = T)
bpx$type[bpx$f.type == "line.body" & bpx$type %in% c("v.dim", "no response")] <- "line.d"
bpx <- bpx[,1:3]

#--------------------------------------------------------
# Now, on with the feature classification

# by default, remove edge pixels & screen spots (if still included) - redundant to cluster these
px <- bpx[!(bpx$type %in% c("edge", "screen.spot", "s.bright", "s.dim", "l.bright", "l.dim")),]

# get clusters without lines to avoid clustering with incidental adjacent pixels
# (globally bright pixels that lie on lines will be retained)
px <- px[!(px$type %in% c("line.b", "line.d")),]

#--------------------------------------------------------
# IDENTIFY CLUSTERS OF DEFECTIVE PIXELS

# convert to raster & get clumps EXCLUDING LINES
cc <- clump(m2r(bpx2im(px, im.dim = c(2000,2000))), dir = 4)

# coordinates of all clumped cells, including unique ID per clump
xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                 id = getValues(cc)[!is.na(getValues(cc))])

# get size of each clump, discard 1-pixel clumps
xy <- merge(xy, count(xy, "id"), all = T)[,c("x", "y", "id", "freq")]
xy <- xy[xy$freq > 1,]
xy$f.type <- "cl.body"

#---------------------------------------------------------------------------
# FIND CLUSTER CENTRES

bpx <- merge(bpx, xy[,c("x", "y", "f.type", "id")], by.x = c("row", "col"), by.y = c("x", "y"), all = T)
bpx$type <- ordered(bpx$type, levels = levels(bpx$type))

for (i in unique(bpx$id[!is.na(bpx$id)])) {
    tmp <- bpx[which(bpx$id == i),]
    
    # if most severe defect is local brightness/dimness, set cluster as singletons
    if (min(tmp$type) %in% c("l.bright", "l.dim")) {
        bpx$f.type[bpx$id == i] <- "singleton"
    } else {
        
        # find brightest pixel in cluster
        cand <- tmp[tmp$type == min(tmp$type),]
        
        if (nrow(cand) > 1) {
            
            # if multiple brightest pixels in cluster, get horizontal midpoint
            cand <- cand[cand$row - floor(mean(tmp$row)) == min(cand$row - floor(mean(tmp$row))),]
            
            if (nrow(cand) > 1) {
                
                # if multiple brightest pixels in midline, get closest to panel edge
                if (min(cand$col) > 992.5) {
                    cand <- cand[which.max(cand$col),]
                } else {
                    cand <- cand[which.min(cand$col),]
                }
            }
        }
        bpx$f.type[bpx$row == cand$row & bpx$col == cand$col] <- "cl.root"
    }
}

#---------------------------------------------------------------------------
# LABEL REMAINING FEATURE TYPES

bpx$f.type[bpx$type %in% c("line.b", "line.d") & is.na(bpx$f.type)] <- "line.body"
bpx$f.type[is.na(bpx$f.type)] <- "singleton"

bpx$f.type <- as.factor(bpx$f.type)

#---------------------------------------------------------------------------

saveRDS(bpx, "./Other-data/Old-data/bad-px-features-131122.rds")

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

