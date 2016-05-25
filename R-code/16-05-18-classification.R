
# CLEAN COPY OF CODE USED TO GENERATE DOCUMENT

####################################################################################################

library("IO.Pixels"); library("CB.Misc")
library(moments)    # for skewness
fpath <- "./Notes/Med-diff-classification/fig/"

load.pixel.means()
load.pixel.sds()
bp <- readRDS(paste0(fpath, "bad-px-maps.rds"))

md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")
cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "gold", "yellow", "grey", "violet", NA, "green3", "green", "lightskyblue", "grey")
headerCat <- gsub("[ ]", "", gsub("[.]", "", cat))
fancyCat <- c("No response", "Dead", "Hot", "V. bright", "Bright", "Bright line", "Locally bright", "Slightly bright", "Screen spot", "Dim line", "Edge", "V. dim", "Dim", "Locally dim", "Slightly dim")

####################################################################################################

# GLOBALLY vs LOCALLY EXTREME PIXELS                                                            ####

# hot, v.bright, v.dim, dead, no response, edge, dim spot                   
{
    bp <- lapply(dimnames(pw.m)[[4]],
                 function(x) rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                                                rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                                       ordered("edge", levels = cat)),
                                            c("row", "col", "type")),
                                   data.frame(no.response(x), type = "no response"),
                                   data.frame(which(pw.m[, , "black", x] == 65535, arr.ind = T), type = "hot"),
                                   data.frame(which(pw.m[, , "white", x] == 0, arr.ind = T), type = "dead"),
                                   screen.spots.xy(x),
                                   get.dim.bright.px(pw.m[,,"white", x]),
                                   get.dim.bright.px(pw.m[,,"grey", x]),
                                   get.dim.bright.px(pw.m[,,"black", x])))

bp <- lapply(lapply(bp, 
                    function(x) x[order(x$type),]),
             function(x) x[!duplicated(x[,1:2]),])
names(bp) <- dimnames(pw.m)[[4]]

saveRDS(bp, paste0(fpath, "bad-px-maps-global.rds"))
}       # ~ 6m to run over all 12 images

bp <- readRDS(paste0(fpath, "bad-px-maps-global.rds"))

# load median-filtered differences
{
    md.b <- readRDS("./Other-data/Median-diffs-black.rds")
    md.g <- readRDS("./Other-data/Median-diffs-grey.rds")
    md.w <- readRDS("./Other-data/Median-diffs-white.rds")
}

# get locally extreme pixels
{
    zz <- lapply(dimnames(pw.m)[[4]],
                 function(x) rbind(setNames(data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T),
                                                        type = "l.bright"), c("x", "y", "type")),
                                    setNames(data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * -2) == 0, arr.ind = T),
                                                        type = "l.dim"), c("x", "y", "type")),
                                   setNames(data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * 2) > 0, arr.ind = T),
                                                       type = "l.bright"), c("x", "y", "type")),
                                   setNames(data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * -2) == 0, arr.ind = T),
                                                       type = "l.dim"), c("x", "y", "type"))))
    zz <- lapply(lapply(zz, 
                        function(x) x[order(x$type),]),
                 function(x) x[!duplicated(x[,1:2]),])
    names(zz) <- dimnames(pw.m)[[4]]
}      # ~ 1m to run over all 12 images

# combine & compare results
{
    px <- lapply(dimnames(pw.m)[[4]],
                 function(x) data.frame(global = ordered(c("normal", cat)[c(bpx2im(bp[[x]])) + 1], levels = c("normal", cat)),
                                        local = ordered(c("normal", "l.bright", "l.dim")[c(bpx2im(zz[[x]])) + 1], levels = c("normal", "l.bright", "l.dim"))))
    names(px) <- dimnames(pw.m)[[4]]
}

write.csv(do.call("cbind", lapply(lapply(px[1:4], table), as.matrix)),
          paste0(fpath, "global-vs-local-first-4.csv"), quote = F)
write.csv(do.call("cbind", lapply(lapply(px[5:8], table), as.matrix)),
          paste0(fpath, "global-vs-local-next-4.csv"), quote = F)
write.csv(do.call("cbind", lapply(lapply(px[9:12], table), as.matrix)),
          paste0(fpath, "global-vs-local-last-4.csv"), quote = F)

####################################################################################################

# GET UNIFIED PIXEL CLASSIFICATION                                                              ####

# ~ 7m to run over all images
{
    bp <- lapply(dimnames(pw.m)[[4]],
                 function(x) rbind(data.frame(edge.px(pw.m), type = ordered("edge", levels = cat)),
                                   data.frame(no.response(x), type = "no response"),
                                   data.frame(which(pw.m[, , "black", x] == 65535, arr.ind = T), type = "hot"),
                                   data.frame(which(pw.m[, , "white", x] == 0, arr.ind = T), type = "dead"),
                                   screen.spots.xy(x),
                                   get.dim.bright.px(pw.m[,,"white", x]),
                                   get.dim.bright.px(pw.m[,,"grey", x]),
                                   get.dim.bright.px(pw.m[,,"black", x]),
                                   data.frame(which(find.lines(pw.m[, , "black", x]) > 0, arr.ind = T), type = "line.b"),
                                   data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * -2) == 0, arr.ind = T), type = "l.dim"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * -2) == 0, arr.ind = T), type = "l.dim")))

    bp <- lapply(lapply(bp, 
                        function(x) x[order(x$type),]),
                 function(x) x[!duplicated(x[,1:2]),])
    names(bp) <- dimnames(pw.m)[[4]]
    
    saveRDS(bp, paste0(fpath, "bad-px-maps.rds"))
}
bp <- readRDS(paste0(fpath, "bad-px-maps.rds"))

# quick summary per image
{
    hh <- setNames(as.data.frame(do.call("cbind", lapply(lapply(lapply(bp, "[", "type"), table), as.matrix))),
                   sapply(names(bp), function (x) paste0(substring(x, 1, 2), "-", substring(x, 3,4), "-", substring(x, 5,6))))
    
    rownames(hh) <- fancyCat
    write.csv(hh, paste0(fpath, "bad-px-identified.csv"), quote = F)
}

####################################################################################################

# TRACK ALL STATE TRANSITIONS                                                                   ####
tr <- list()

for (i in 1:(length(bp) - 1)) {
    tr[[i]] <- table("From" = ordered(c("normal", Cat)[c(bpx2im(bp[[i]])) + 1], levels = c("normal", Cat)),
                   "To" = ordered(c("normal", Cat)[c(bpx2im(bp[[i+1]])) + 1], levels = c("normal", Cat)))
    names(tr)[[i]] <- paste(names(bp)[c(i,(i+1))], collapse = "-")
}

# write CSVs of numbers of pixels transitioning
for (i in 1:length(tr)) {
    tmp <- tr[[i]][-11, -11]
    tmp[tmp == 0] <- "-"
    write.csv(tmp, paste0(fpath, "bpx-change-", names(tr)[i], ".csv"), quote = F)
}

# write CSVs of proportion of pixels transitioning
for (i in 1:length(tr)) {
    tmp <- 100 * tr[[i]][-11, -11] / rowSums(tr[[i]][-11, -11])
    tmp[tmp == 0] <- NA
    tmp <- round(tmp, 1)
    tmp[is.na(tmp)] <- "-"
    write.csv(tmp, paste0(fpath, "bpx-change-", names(tr)[i], "-prop.csv"), quote = F)
}

# mean transition rates
tr <- array(unlist(tr), dim = c(16, 16, 11), 
            dimnames = list(c("normal", Cat), c("normal", Cat), names(tr)))

write.csv(prep.csv(apply(tr, 1:2, mean), dp = 0)[-11, -11],
          paste0(fpath, "transition-px-mean.csv"), quote = F)

write.csv(prep.csv(apply(tr, 1:2, sd), dp = 0)[-11, -11],
          paste0(fpath, "transition-px-sd.csv"), quote = F)

write.csv(prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, mean, na.rm = T) * 100)[-11, -11],
          paste0(fpath, "transition-prop-mean.csv"), quote = F)

write.csv(prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, sd, na.rm = T) * 100)[-11, -11],
          paste0(fpath, "transition-prop-sd.csv"), quote = F)


####################################################################################################

# TABLE OF MEDIAN/SD/SKEW PER IMAGE                                                             ####

zz <- data.frame(med.black = apply(pw.m[,,"black",], 3, median),
                 mad.black = apply(pw.m[,,"black",], 3, mad),
                 sd.black = apply(pw.m[,,"black",], 3, sd),
                 skew.black = apply(apply(pw.m[,,"black",], 3, c), 2, skewness),
                 med.grey = apply(pw.m[,,"grey",], 3, median),
                 mad.grey = apply(pw.m[,,"grey",], 3, mad),
                 sd.grey = apply(pw.m[,,"grey",], 3, sd),
                 skew.grey = apply(apply(pw.m[,,"grey",], 3, c), 2, skewness),
                 med.white = apply(pw.m[,,"white",], 3, median),
                 mad.white = apply(pw.m[,,"white",], 3, mad),
                 sd.white = apply(pw.m[,,"white",], 3, sd),
                 skew.white = apply(apply(pw.m[,,"white",], 3, c), 2, skewness))

zz <- round(zz, 1)
zz[,-c(4, 8, 12)] <- round(zz[,-c(4, 8, 12)], 0)
rownames(zz) <- sapply(rownames(zz), function (x) paste(substring(x, 1, 2), substring(x, 3, 4), substring(x, 5, 6), sep = "-"))

write.csv(zz, paste0(fpath, "image-summaries.csv"), quote = F)


# plot each histogram with its predecessor
{
    pdf(paste0(fpath, "img-hists-black"))
    par(mfrow = c(4, 3), mar = c(2,2,1,2))
    
    hist(pw.m[,,"black", 1], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), ylim = c(0, 40000), xlim = c(4000, 7000), main = "")
    legend("topleft", col = "magenta3", pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1])
    
    hist(pw.m[,,"black", 1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), col = adjustcolor("cyan3", alpha = 0.4), ylim = c(0, 40000), xlim = c(4000, 7000), main = "")
    hist(pw.m[,,"black", 2], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), add = T)
    legend("topleft", col = c("cyan3", "magenta3"), pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1:2])    
    
    for (i in 3:12) {
        hist(pw.m[,,"black", i-2], breaks = "fd", border = adjustcolor("gold", alpha = 0.4), ylim = c(0, 40000), xlim = c(4000, 7000),
             main = "", xlab = "", ylab = "")
        hist(pw.m[,,"black", i-1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), add = T)
        hist(pw.m[,,"black", i], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), add = T)
        legend("topleft", col = c("gold", "cyan3", "magenta3"), pch = 15, bty = "n",
               legend = dimnames(pw.m)[[4]][(i-2):i])
    }
    par(mfrow = c(1,1))
    dev.off()
}

{
    pdf(paste0(fpath, "img-hists-grey.pdf"))
    par(mfrow = c(4, 3), mar = c(2,2,1,2))
    
    hist(pw.m[,,"grey", 1], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), ylim = c(0, 40000), xlim = c(15000, 24000), main = "")
    legend("topleft", col = "magenta3", pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1])
    
    hist(pw.m[,,"grey", 1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), col = adjustcolor("cyan3", alpha = 0.4), ylim = c(0, 40000), xlim = c(15000, 24000), main = "")
    hist(pw.m[,,"grey", 2], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), add = T)
    legend("topleft", col = c("cyan3", "magenta3"), pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1:2])    
    
    for (i in 3:12) {
        hist(pw.m[,,"grey", i-2], breaks = "fd", border = adjustcolor("gold", alpha = 0.4), ylim = c(0, 35000), xlim = c(15000, 24000),
             main = "", xlab = "", ylab = "")
        hist(pw.m[,,"grey", i-1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), add = T)
        hist(pw.m[,,"grey", i], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), add = T)
        legend("topleft", col = c("gold", "cyan3", "magenta3"), pch = 15, bty = "n",
               legend = dimnames(pw.m)[[4]][(i-2):i])
    }
    par(mfrow = c(1,1))
    dev.off()
}

{
    pdf(paste0(fpath, "img-hists-white.pdf"))
    par(mfrow = c(4, 3), mar = c(2,2,1,2))
    hist(pw.m[,,"white", 1], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), ylim = c(0, 55000), xlim = c(35000, 55000), main = "")
        legend("topleft", col = "magenta3", pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1])
    
    hist(pw.m[,,"white", 1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), col = adjustcolor("cyan3", alpha = 0.4), ylim = c(0, 55000), xlim = c(35000, 55000), main = "")
    hist(pw.m[,,"white", 2], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), add = T)
    legend("topleft", col = c("cyan3", "magenta3"), pch = 15, bty = "n", legend = dimnames(pw.m)[[4]][1:2])    
    for (i in 3:12) {
        hist(pw.m[,,"white", i-2], breaks = "fd", border = adjustcolor("gold", alpha = 0.4), col = adjustcolor("gold", alpha = 0.4), ylim = c(0, 55000), xlim = c(35000, 55000),
             main = "", xlab = "", ylab = "")
        hist(pw.m[,,"white", i-1], breaks = "fd", border = adjustcolor("cyan3", alpha = 0.4), col = adjustcolor("cyan3", alpha = 0.4), add = T)
        hist(pw.m[,,"white", i], breaks = "fd", border = adjustcolor("magenta3", alpha = 0.4), col = adjustcolor("magenta3", alpha = 0.4), add = T)
        legend("topleft", col = c("gold", "cyan3", "magenta3"), pch = 15, bty = "n",
               legend = dimnames(pw.m)[[4]][(i-2):i])
        }
    par(mfrow = c(1,1))
    dev.off()
}

# STABILITY OF CLASSES WITHIN EACH ACQUISITION                                                  ####

d <- load.daily(160430)

# get median differences for all images (relatively slow process: ~20m per 20 images)
{
    md.b <- apply(d[,,,"black"], 3, med.diffs)
    md.g <- apply(d[,,,"grey"], 3, med.diffs)
}

md.b <- array(readRDS(paste0(fpath, "med-diffs-20-160430-black.rds")), dim = c(1996, 1996, 20))
md.g <- array(readRDS(paste0(fpath, "med-diffs-20-160430-grey.rds")), dim = c(1996, 1996, 20))

# temporary functions to decouple from pw.m
{
    no.resp <- function(m, limit = 0.01) {
        bn <- qJohnson(c(limit, 1 - limit), JohnsonFit(d[,,m,"black"]))
        un <- rbind(which(matrix(findInterval(d[,,m,"grey"], bn), ncol = 1996) == 1, arr.ind = T), 
                    which(matrix(findInterval(d[,,m,"white"], bn), ncol = 1996) == 1, arr.ind = T))
        return(un[!duplicated(un), ])
    }
    
    spots <- function(im, smooth.span = 1/5, min.diam = 5, edge.width = 10, auto.threshold = T) {
    
        sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
        
        # apply lowess smoothing
        smoo <- lowess.per.column(im, span = smooth.span)
        res <- im - smoo
        
        # flatten further by setting brighter pixels to mean value 
        res[res > mean(res)] <- mean(res)
        
        # dilate resulting image
        dilated <- dilate(res, sk)
        
        # erode resulting image (complete morphological closing)
        eroded <- erode(dilated, sk)
        
        # use k-means thresholding to identify spots
        # use 1-thresholded value to assign 1 to spots, 0 to background
        if (auto.threshold) {
            dim <- array(1, dim = c(1996, 1996)) - threshold(eroded, method = "kmeans")
        } else {
            dim <- array(1, dim = c(1996, 1996)) - threshold(eroded, mean(eroded) - 3*sd(eroded))
        }
        
        #------------------------------------------------------------------------------
        # convert to raster & clump values to identify individual spots
        # (values need to be reordered to maintain image orientation)
        blobs <- clump(raster(t(dim[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5), dir = 4)
        
        # summarise & remove any clusters that are too close to an edge
        # check location of cluster (looking for edge clusters)
        sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                               id = getValues(blobs)[!is.na(getValues(blobs))]),
                    .(id), summarise, 
                    xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y),
                    xm = round(mean(x),0), ym = round(mean(y),0))
        
        sc$n.id <- sc$id
        sc$n.id[(sc$ym >= 1996 - edge.width) | (sc$ym <= edge.width) | 
                    (sc$xm >= 1996 - edge.width) | (sc$xm <= edge.width)] <- NA
        
        zz <- subs(blobs, sc[,c("id", "n.id")])
        
        if (all(is.na(getValues(zz)))) {
            # return empty data frame in order to not 
            return(data.frame("row" = double(), "col" = double(), "type" = character()))
        } else {
            # extract coordinates of cells covered by spots identified, return as data frame
            return(setNames(data.frame(xyFromCell(zz, which(!is.na(getValues(zz)))), type = "screen spot"),
                            c("row", "col", "type")))
        }
    }  
}

# get bad pixel map for all 20 images
bp <- list()
for (n in 1:20) {
    bp[[n]] <- rbind(data.frame(edge.px(d), type = ordered("edge", levels = cat)),
                     data.frame(no.resp(n), type = "no response"),
                     data.frame(which(d[ , , n, "black"] == 65535, arr.ind = T), type = "hot"),
                     data.frame(which(d[ , , n, "white"] == 0, arr.ind = T), type = "dead"),
                     spots(d[ , , n, "white"]),
                     get.dim.bright.px(d[ , , n, "black"]),
                     get.dim.bright.px(d[ , , n, "grey"]),
                     get.dim.bright.px(d[ , , n, "white"]),
                     data.frame(which(find.lines(d[ , , n, "black"]) > 0, arr.ind = T), type = "line.b"),
                     data.frame(which(threshold(md.b[,,n], level = mad(d[ , , n, "black"]) * 2) > 0, arr.ind = T), type = "l.bright"),
                     data.frame(which(threshold(md.b[,,n], level = mad(d[ , , n, "black"]) * -2) == 0, arr.ind = T), type = "l.dim"),
                     data.frame(which(threshold(md.g[,,n], level = mad(d[ , , n, "grey"]) * 2) > 0, arr.ind = T), type = "l.bright"),
                     data.frame(which(threshold(md.g[,,n], level = mad(d[ , , n, "grey"]) * -2) == 0, arr.ind = T), type = "l.dim"))
}

bp <- lapply(lapply(bp, 
                    function(x) x[order(x$type),]),
             function(x) x[!duplicated(x[,1:2]),])

zz <- array(unlist(lapply(bp, bpx2im)), dim = c(1996, 1996, 20))

tr <- array(dim = c(15, 15, 19), 
            dimnames = list(c("normal", cat[-10]), c("normal", cat[-10]), NULL))

# get changes at each step
for (i in 1:19) {
    tr[,,i] <- table("from" = zz[,,i], "to" = zz[,,i+1])
}

tr.mean <- prep.csv(apply(tr, 1:2, mean), dp = 0)
tr.sd <- prep.csv(apply(tr, 1:2, sd), dp = 0)
tr.prop <- prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, mean, na.rm = T) * 100)
tr.prop.sd <- prep.csv(apply(array(apply(tr, 3, function(x) x / rowSums(x)), dim = dim(tr), dimnames = dimnames(tr)), 1:2, sd, na.rm = T) * 100)

####################################################################################################

# STANDARD DEVIATIONS WITHIN EACH CLASS                                                         ####

load.pixel.sds()

# suspect that slightly bright/slightly dim may have higher than expected SD?

summ.by.px.type <- function(dt) {
    
    px <- data.frame(bp[[dt]],
                     pw.sd.b = pw.sd[,,"black", dt][as.matrix(bp[[dt]][,1:2])],
                     pw.sd.g = pw.sd[,,"grey", dt][as.matrix(bp[[dt]][,1:2])],
                     pw.sd.w = pw.sd[,,"white", dt][as.matrix(bp[[dt]][,1:2])])
                
    df <- data.frame(pw.sd.b = aggregate(pw.sd.b ~ type, data = px, mean)[,2],
                           pw.sd.g = aggregate(pw.sd.g ~ type, data = px, mean)[,2],
                           pw.sd.w = aggregate(pw.sd.w ~ type, data = px, mean)[,2],
                           b.sd = aggregate(pw.sd.b ~ type, data = px, sd)[,2],
                           g.sd = aggregate(pw.sd.g ~ type, data = px, sd)[,2],
                           w.sd = aggregate(pw.sd.w ~ type, data = px, sd)[,2],
                           b.mad = aggregate(pw.sd.b ~ type, data = px, mad)[,2],
                           g.mad = aggregate(pw.sd.g ~ type, data = px, mad)[,2],
                           w.mad = aggregate(pw.sd.w ~ type, data = px, mad)[,2],
                           n.px = aggregate(pw.sd.w ~ type, data = px, length)[,2])
    rownames(df) <- c(cat[-10])
    
    return(df)
}

summ <- summ.by.px.type("160430")

# plot histograms of SD for each pixel type
{

    px <- data.frame(bp[[dt]],
                     pw.sd.b = pw.sd[,,"black", dt][as.matrix(bp[[dt]][,1:2])],
                     pw.sd.g = pw.sd[,,"grey", dt][as.matrix(bp[[dt]][,1:2])],
                     pw.sd.w = pw.sd[,,"white", dt][as.matrix(bp[[dt]][,1:2])])
        
    sd.hist <- function(class, x.lim, y.lim) {
        
        if (missing(x.lim)) x.lim <- max(pretty(range(px[px$type == class, 4:6])))
        
        hs <- hist(unlist(px[px$type == class, 4:6]), breaks = "fd", plot = F)
        if (missing(y.lim)) y.lim <- max(pretty(range(hs$counts)))
        
        hist(px[px$type == class, 6],  xlab = "SD", ylab = "", breaks = hs$breaks, prob = T, xlim = c(0, x.lim),
             col = "gold", border = NA, ylim = c(0, y.lim / mean(hs$counts / hs$density, na.rm = T)), yaxt = "n",
             main = paste0(fancyCat[which(cat == class)], " - ", length(which(px$type == class)), " pixels"))
        hist(px[px$type == class, 5], add = T, col = adjustcolor("cyan3", alpha = 0.4), border = NA, prob = T, breaks = hs$breaks)
        hist(px[px$type == class, 4], add = T, col = adjustcolor("magenta3", alpha = 0.4), border = NA, prob = T, breaks = hs$breaks)
        
        axis(2, at = pretty(c(0:y.lim)) / mean(hs$counts / hs$density, na.rm = T),
             labels = pretty(c(0:y.lim)))
        
        s <- sample.healthy(px[, 1:2], n = nrow(px[px$type == class, 1:2]))
        
        lines(density(pw.sd[,,"white", dt][s]), col = "orangered")
        lines(density(pw.sd[,,"grey", dt][s]), col = "blue")
        lines(density(pw.sd[,,"black", dt][s]), col = "purple")
        
        legend("topright", col = c("magenta3", "cyan3", "gold"), pch = 15,
               legend = c("black", "grey", "white"), bty = "n")
    }
    
    sd.hist("no response", x.lim = 100)
    sd.hist("hot")
    sd.hist("v.bright", x.lim = 2500, y.lim = 300)
    sd.hist("bright", x.lim = 1500, y.lim = 300)
    sd.hist("l.bright", x.lim = 1500, y.lim = 8000)
    sd.hist("s.bright", x.lim = 1500, y.lim = 8000)
    sd.hist("l.dim", x.lim = 700, y.lim = 600)
    sd.hist("s.dim", x.lim = 700, y.lim = 300)
}

#---------------------------------------------------------------------------------------------
# increasing trend across acquisition?
length(which(d[,,20, "black"] > d[,,1, "black"])) / (1996^2)    # 0.449955
length(which(d[,,20, "grey"] > d[,,1, "grey"])) / (1996^2)      # 0.4054218

length(which(d[,,20, "white"] > d[,,1, "white"])) / (1996^2)    # 0.9324614
length(which(d[,,20, "white"] > d[,,2, "white"])) / (1996^2)    # 0.5839326
length(which(d[,,2, "white"] > d[,,1, "white"])) / (1996^2)     # 0.8995466
#---------------------------------------------------------------------------------------------

# plot transects of interesting points
{
    d <- load.daily(160430)
    
    matplot(d[1880,1872,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    # hits maximum in black & white, but not in grey...?
    # grey @ 13:23, white @ 13:21, black @ 14:43. Possibly some other usage before black images led to residual current?
    # may be able to confirm this after checking .xtx files
    
    matplot(d[9,1991,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    # v. erratic point, prone to big drops in value in both black & white images
    
    matplot(d[488,763,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[32,1000,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[40,1923,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[32,192,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[1945, 1150,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[75, 1742,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))

    matplot(d[8, 6,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[1587, 145,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[323, 178,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    matplot(d[721, 127,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
    
}
matplot(d[592, 1118,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
matplot(d[1496, 193,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
matplot(d[789, 774,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))
matplot(d[116, 48,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))

matplot(d[378, 1,,], type = "o", pch = 20, col = c("magenta3", "cyan3", "gold"))

####################################################################################################

# STANDARD DEVIATION OF 'NORMAL' PIXELS                                                         ####

# does anything useful come of identifying 'locally unusual' pixels?

s.hist(pw.sd[,,"black", "160430"])
pixel.image(pw.sd[,,"black", "160430"])
# no point in using median-differencing: pattern is not so systematic

# try Johnson quantiles instead
# 'official' definition is > 6 sigma from median

pdf(paste0(fpath, "sd-hist-black-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"black", "160430"], xlim = c(0, 100), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"black", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"black", "160430"]) + (6 * sd(pw.sd[,,"black", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"black", "160430"]) + (6 * mad(pw.sd[,,"black", "160430"])), col = "red", lty = 2)
    dev.off()
}
  
pdf(paste0(fpath, "sd-hist-grey-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"grey", "160430"], xlim = c(0, 400), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"grey", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"grey", "160430"]) + (6 * sd(pw.sd[,,"grey", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"grey", "160430"]) + (6 * mad(pw.sd[,,"grey", "160430"])), col = "red", lty = 2)
    dev.off()
}  

pdf(paste0(fpath, "sd-hist-white-160430.pdf"), width = 7, height = 4); {
    par(mar = c(2,2,1,1))
    s.hist(pw.sd[,,"white", "160430"], xlim = c(0, 700), main = "", xlab = "", ylab = "")
    abline(v = qJohnson(c(0.999, 0.9995, 0.99999), JohnsonFit(pw.sd[,,"white", "160430"])), col = "cyan3", lty = c(1:3))
    abline(v = median(pw.sd[,,"white", "160430"]) + (6 * sd(pw.sd[,,"white", "160430"])), col = "red")
    abline(v = median(pw.sd[,,"white", "160430"]) + (6 * mad(pw.sd[,,"white", "160430"])), col = "red", lty = 2)
    dev.off()
}  

    length(which(pw.sd[,,"black", "160430"] > qJohnson(0.9999, JohnsonFit(pw.sd[,,"black", "160430"]))))
    length(which(pw.sd[,,"black", "160430"] > median(pw.sd[,,"black", "160430"]) + (6 * mad(pw.sd[,,"black", "160430"]))))

# Normal & Johnson Q-Q plots of SD at each power setting
{
    QQ.norm <- function(data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99),...) {
        plot(qnorm(quantiles, mean(data), sd(data)),
             quantile(data, quantiles), pch = 20, asp = T, ylab = "Observed quantile", xlab = "Normal quantile", ...)
        abline(0,1,col = "red")
        
        abline(h = quantile(data, grid.quantiles), col = "skyblue", 
               lty = 2)
        abline(v = qnorm(grid.quantiles, mean(data), sd(data)), col = "skyblue", 
               lty = 2)
    }
        
        QQ.norm(pw.sd[,,"black", "160430"], grid.quantiles = c(0.999, 0.99))     # very very not normal
        QQ.norm(pw.sd[,,"grey", "160430"], grid.quantiles = c(0.999, 0.99))      # close to normal
        QQ.norm(pw.sd[,,"white", "160430"], grid.quantiles = c(0.999, 0.99))     # close to normal
        
        Johnson.QQ(pw.sd[,,"black", "160430"], grid.quantiles = c(0.999, 0.99))
        Johnson.QQ(pw.sd[,,"grey", "160430"], grid.quantiles = c(0.999, 0.99))
        Johnson.QQ(pw.sd[,,"white", "160430"], grid.quantiles = c(0.999, 0.99))
}

# summarise standard deviations
{
    sd.summ <- apply(pw.sd, 4, 
                     function(y) data.frame(mean = apply(y, 3, mean),
                                            median = apply(y, 3, median),
                                            sd = apply(y, 3, sd),
                                            mad = apply(y, 3, mad),
                                            th = apply(y, 3, function(x) median(x) + 6 * sd(x)),
                                            noisy1 = apply(y, 3, 
                                                           function(x) length(which(x > median(x) + 6 * sd(x)))),
                                            q999 = apply(y, 3, function(x) qJohnson(0.999, JohnsonFit(x))),
                                            noisy999 = apply(y, 3, 
                                                             function(x) length(which(x > qJohnson(0.999, JohnsonFit(x))))),
                                            q9999 = apply(y, 3, function(x) qJohnson(0.9999, JohnsonFit(x))),
                                            noisy9999 = apply(y, 3, 
                                                              function(x) length(which(x > qJohnson(0.9999, JohnsonFit(x)))))))
    
    sd.summ <- lapply(sd.summ, round, 1)
    
    for (dt in names(sd.summ)) {
        write.csv(sd.summ[[dt]], paste0(fpath, "sd-summary-", dt, ".csv"), quote = F)
    }
}

# ratio of normal grey/white SD to extreme black SD
    unlist(lapply(sd.summ, function(x) x$median[2] / x$th[1]))
    unlist(lapply(sd.summ, function(x) x$median[3] / x$th[1]))
    
# compare 'official' threshold with 
    unlist(lapply(sd.summ, function(x) x$q9999[2] / x$th[2]))
    unlist(lapply(sd.summ, function(x) x$q9999[3] / x$th[3]))
    
# identify noisy px in white & grey images: same px identified in both?
    
zz <- bpx2im(data.frame(which(pw.sd[,,"grey", "160430"] > qJohnson(0.9999, JohnsonFit(pw.sd[,,"grey", "160430"])), arr.ind = T),
            type = "grey"))
    
image(1:1996, 1:1996, zz)

# relationship between black & grey, black & white, grey & white
bpx <- bpx2im(bp$"160430")
thresholds <- list(black = c(median(pw.sd[,,"black", "160430"]) + 6 * sd(pw.sd[,,"black", "160430"]),
                             qJohnson(0.9999, JohnsonFit(pw.sd[,,"black", "160430"]))),
                   grey = c(median(pw.sd[,,"grey", "160430"]) + 6 * sd(pw.sd[,,"grey", "160430"]),
                             qJohnson(0.9999, JohnsonFit(pw.sd[,,"grey", "160430"]))),
                   white = c(median(pw.sd[,,"white", "160430"]) + 6 * sd(pw.sd[,,"white", "160430"]),
                             qJohnson(0.9999, JohnsonFit(pw.sd[,,"white", "160430"]))))

plot.sp.sd <- function(sp, ps = c("grey", "white")) {
    plot(c(subpanels(pw.sd[,,ps[1], "160430"])[,,sp]), 
         c(subpanels(pw.sd[,,ps[2], "160430"])[,,sp]), 
         xlim = c(0,5000), ylim = c(0,5000),
         col = c("black", cat.cols)[subpanels(bpx)[,,sp] + 1],
         pch = 20, xlab = paste0(ps[1], " SD"), ylab = paste0(ps[2], " SD"), 
         main = paste0("Subpanel ", dimnames(subpanels(pw.sd[,,"black", "160430"]))[[3]][sp]))
    
    abline(h = thresholds[[ps[2]]], col = c("red", "cyan3"), lty = 2)
    abline(v = thresholds[[ps[1]]], col = c("red", "cyan3"), lty = 2)
}

pdf(paste0(fpath, "sd-by-subpanel-grey-white.pdf")); {
    plot.sp.sd(1)
    dev.off()
}

####################################################################################################

# PLOT SD PER PIXEL TYPE                                                                        ####

fpath <- "./Notes/Med-diff-classification/fig/"
load.pixel.sds()
bp <- readRDS(paste0(fpath, "bad-px-maps.rds"))

# for pixels identified as noisy in any dimension, what is behaviour in shading correction?

# standard deviation 
range(unlist(alply(subpanels(pw.sd[,,"black", "160430"]), 3, max, na.rm = T)))   # 501:10001 (1px > 5000)
range(unlist(alply(subpanels(pw.sd[,,"grey", "160430"]), 3, max, na.rm = T)))   # 405:4502 (0px > 5000)
range(unlist(alply(subpanels(pw.sd[,,"white", "160430"]), 3, max, na.rm = T)))   # 546:6297 (2px > 5000)



####################################################################################################

# add noisy pixels into bad pixel map
# what is the SD distribution of that group?
# what is the mean value distribution of that group?
# track per-pixel performance over single acquisition (over all acquisitions?)

####################################################################################################

# IDENTIFY LOCALLY & GLOBALLY BRIGHT/DIM PX IN SHADING-CORRECTED IMAGE. ANY NEW ONES?           ####