
pw.m <- readRDS("./Other-data/Old-data/Pixelwise-means.rds")

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
    tmp.g <- readTIFF(paste0(fpath, "/BadPixelMapBlack.tif"), as.is = T)
    tmp.w <- readTIFF(paste0(fpath, "/BadPixelMapBlack.tif"), as.is = T)
    
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

saveRDS(pw.m, "./Other-data/Old-data/Pixelwise-means.rds", as.is = T)

# MEDIAN DIFFERENCES

md.b <- array(unlist(apply(pw.m[,,"black", ], 3, med.diffs), dim = dim(pw.m)))
md.g <- array(unlist(apply(pw.m[,,"grey", ], 3, med.diffs), dim = dim(pw.m)))
md.w <- array(unlist(apply(pw.m[,,"white", ], 3, med.diffs), dim = dim(pw.m)))

saveRDS(md.b, "./Other-data/Old-data/Med-diffs-black.rds")
saveRDS(md.g, "./Other-data/Old-data/Med-diffs-grey.rds")
saveRDS(md.w, "./Other-data/Old-data/Med-diffs-white.rds")

####################################################################################################

# PIXEL CLASSIFICATION                                                                          ####

md.b <- readRDS("./Other-data/Old-data/Med-diffs-black.rds")
md.g <- readRDS("./Other-data/Old-data/Med-diffs-grey.rds")
md.w <- readRDS("./Other-data/Old-data/Med-diffs-white.rds")

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
                                   data.frame(which(threshold(md.b[,,x], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * -2) == 0, arr.ind = T), type = "l.dim"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * -2) == 0, arr.ind = T), type = "l.dim")))
    
    bp <- lapply(lapply(bp, 
                        function(x) x[order(x$type),]),
                 function(x) x[!duplicated(x[,1:2]),])
    names(bp) <- dimnames(pw.m)[[4]]
    
    saveRDS(bp, "./Other-data/Old-data/bad-px-maps.rds")
}
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")

####################################################################################################

# STATE TRANSITIONS                                                                             ####

####################################################################################################

# STANDARD DEVIATION ACROSS EACH ACQUISITION                                                    ####

# plot SD (per subpanel?)

####################################################################################################

# PIXEL SHADING CORRECTION                                                                      ####

####################################################################################################