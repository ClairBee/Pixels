
library("IO.Pixels"); library("CB.Misc")

pw.m <- readRDS("./Other-data/Old-data/Pixelwise-means.rds")
Cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "gold", "yellow", "grey", "violet", NA, "green3", "green", "lightskyblue", "grey")


p.dim <- panel.edges(left.crop = 24, upper.crop = 24, x.dim = 2000, y.dim = 2000)

xml.df <- read.csv("./Other-data/Old-data/xml-summary.csv", as.is = T, row.names = 1)

####################################################################################################

# need new function to obtain subpanels (or update old one to accept flexible image size)

# second version of images is 1000x1000. May be able to shed light on compression used?

####################################################################################################

# READ XML TO CHECK IMAGE CONSISTENCY                                                           ####

{
    folders <- list.dirs("./Other-data/Old-data", full.names = F, recursive = F)
    folders <- folders[!(folders == "Misc")]
    
    xml.df <- rbind.fill(lapply(folders, 
                    function(x) as.data.frame(lapply(do.call("c",
                                                             c(xmlToList(xmlParse(paste0("./Other-data/Old-data/", x, "/CalibrationParameters.xml"))),
                                                                     list(recursive = TRUE))), 
                                          FUN = unlist), stringsAsFactors = F)))
    rownames(xml.df) <- folders
    write.csv(xml.df, "./Other-data/Old-data/xml-summary.csv", quote = F)
}

####################################################################################################

# READ IMAGE DATA INTO ARRAY                                                                    ####

folders <- list.dirs("./Other-data/Old-data", full.names = F, recursive = F)
folders <- folders[!(folders == "Misc")]

load.images <- function(dt, direct = F) {
    dt <- toString(dt)
    fpath <- paste0("./Other-data/Old-data/", dt)
    tmp.b <- readTIFF(paste0(fpath, "/BadPixelMapBlack.tif"), as.is = T)
    tmp.g <- readTIFF(paste0(fpath, "/BadPixelMapGrey.tif"), as.is = T)
    tmp.w <- readTIFF(paste0(fpath, "/BadPixelMapWhite.tif"), as.is = T)
    
    if (direct) {
        pw.m[,,"black", dt] <<- t(tmp.b[nrow(tmp.b):1, , drop = FALSE])
        pw.m[,,"grey", dt] <<- t(tmp.g[nrow(tmp.g):1, , drop = FALSE])
        pw.m[,,"white", dt] <<- t(tmp.w[nrow(tmp.w):1, , drop = FALSE])
    } else {
        tmp.b <<- t(tmp.b[nrow(tmp.b):1, , drop = FALSE])
        tmp.g <<- t(tmp.g[nrow(tmp.g):1, , drop = FALSE])
        tmp.w <<- t(tmp.w[nrow(tmp.w):1, , drop = FALSE])
    }
}

# may need some padding to properly align panels, so will do individually
pw.m <- array(dim = c(2000, 2000, 3, 6),
               dimnames = list(NULL, NULL, c("black", "grey", "white"), folders))

load.images(130613, direct = T)
load.images(130701, direct = T)
load.images(131002, direct = F)

# pad array before assigning to main array
pw.m[,,,"131002"] <- abind(array(dim = c(2000, 200, 3)),
                           abind(tmp.b, tmp.g, tmp.w, along = 3),
                           array(dim = c(2000, 200, 3)),
                           along = 2)
             
load.images(131122, direct = F)

# pad array before assigning to main array
pw.m[,,,"131122"] <- abind(array(dim = c(2000, 200, 3)),
                           abind(tmp.b, tmp.g, tmp.w, along = 3),
                           array(dim = c(2000, 200, 3)),
                           along = 2)

load.images(140128, direct = T)
load.images(140129, direct = T)

saveRDS(pw.m, "./Other-data/Old-data/Pixelwise-means.rds")

# MEDIAN DIFFERENCES

md.b <- array(unlist(apply(pw.m[,,"black", ], 3, med.diffs)), dim = dim(pw.m[,,"black",]), dimnames = dimnames(pw.m[,,"black",]))
md.g <- array(unlist(apply(pw.m[,,"grey", ], 3, med.diffs)), dim = dim(pw.m[,,"grey", ]), dimnames = dimnames(pw.m[,,"grey",]))
md.w <- array(unlist(apply(pw.m[,,"white", ], 3, med.diffs)), dim = dim(pw.m[,,"white", ]), dimnames = dimnames(pw.m[,,"white",]))

saveRDS(md.b, "./Other-data/Old-data/Med-diffs-black.rds")
saveRDS(md.g, "./Other-data/Old-data/Med-diffs-grey.rds")
saveRDS(md.w, "./Other-data/Old-data/Med-diffs-white.rds")

####################################################################################################

# PIXEL CLASSIFICATION                                                                          ####

md.b <- readRDS("./Other-data/Old-data/Med-diffs-black.rds")
md.g <- readRDS("./Other-data/Old-data/Med-diffs-grey.rds")
md.w <- readRDS("./Other-data/Old-data/Med-diffs-white.rds")
md <- abind(md.b, md.g, md.w, along = 2.5, new.names = dimnames(pw.m))

{
    bp <- lapply(dimnames(pw.m)[[4]],
                 function(x) rbind(data.frame(edge.px(pw.m), type = ordered("edge", levels = Cat)),
                                   data.frame(no.response(x), type = "no response"),
                                   hot.px(x),
                                   dead.px(x),
                                   get.dim.bright.px(pw.m[,,"white", x]),
                                   get.dim.bright.px(pw.m[,,"grey", x]),
                                   get.dim.bright.px(pw.m[,,"black", x]),
                                   data.frame(which(find.lines(pw.m[, , "black", x], midline = 1000.5) > 0, arr.ind = T), type = "line.b"),
                                   data.frame(which(find.lines(pw.m[, , "grey", x], dim.lines = T, midline = 1000.5) > 0, arr.ind = T), type = "line.d"),
                                   locally.bright.px(x, "black"),
                                   locally.dim.px(x, "black"),
                                   locally.bright.px(x, "grey"),
                                   locally.dim.px(x, "grey")))

bp <- lapply(lapply(bp, 
                        function(x) x[order(x$type),]),
                 function(x) x[!duplicated(x[,1:2]),])
    names(bp) <- dimnames(pw.m)[[4]]
    
    saveRDS(bp, "./Other-data/Old-data/bad-px-maps.rds")
}
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")

# single-image bad pixel list, while working out details
{
    bp <- list()
    {
        bp$edge <- data.frame(edge.px(pw.m[,,"black", "131002"]), type = ordered("edge", levels = Cat))
        bp$noResp <- data.frame(no.response("131002"), type = "no response")
        bp$hot <- hot.px("131002")
        bp$dead <- dead.px("131002")
        # no dim spots visible in images, so don't run (need to update function to capture this)
        bp$db <- rbind(get.dim.bright.px(pw.m[,,"black", "131002"]),
                       get.dim.bright.px(pw.m[,,"grey", "131002"]),
                       get.dim.bright.px(pw.m[,,"white", "131002"]))
        bp$line.b <- data.frame(which(find.lines(pw.m[, , "black", "131002"], midline = 1000.5) > 0, arr.ind = T), type = "line.b")
        bp$local.b <- locally.bright.px("131002", "black")
        bp$local.g <- locally.bright.px("131002", "grey")
        bp$local.w <- locally.bright.px("131002", "white")
        bp$local.b.d <- locally.dim.px("131002", "black")
        bp$local.g.d <- locally.dim.px("131002", "grey")
        bp$local.w.d <- locally.dim.px("131002", "white")
    }
    
    sapply(bp, nrow)
    
    all.bp <- rbind.fill(bp)
    all.bp <- all.bp[order(all.bp$type),]
    all.bp <- all.bp[!duplicated(all.bp[,1:2]),]
    
    table(all.bp$type)
    
    plot(bp$"140129", pch = 15, cex = 0.4, col = Cat.cols[bp$"140129"$type])
}


####################################################################################################

# STATE TRANSITIONS                                                                             ####

####################################################################################################

# STANDARD DEVIATION ACROSS EACH ACQUISITION                                                    ####

# plot SD (per subpanel?)

####################################################################################################

# PIXEL SHADING CORRECTION                                                                      ####

####################################################################################################