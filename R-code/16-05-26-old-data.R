
library("IO.Pixels"); library("CB.Misc")

pw.m <- readRDS("./Other-data/Old-data/Pixelwise-means.rds")
bp <- readRDS("./Other-data/Old-data/bad-px-maps.rds")
fpath <- "./Notes/Old-data/fig/"

Cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "yellow", NA, "grey", "violet", NA, "green3", "green", "lightskyblue", NA)

p.dims <- list("130613" = panel.edges(left.crop = 24, upper.crop = 24, x.dim = 2000, y.dim = 2000),
               "130701" = panel.edges(left.crop = 24, upper.crop = 24, x.dim = 2000, y.dim = 2000),
               "131002" = panel.edges(left.crop = 24, upper.crop = 24, x.dim = 2000, y.dim = 2000),
               "131122" = panel.edges(left.crop = 24, upper.crop = 24, x.dim = 2000, y.dim = 2000),
               "140128" = panel.edges(left.crop = 0, upper.crop = 20, x.dim = 2000, y.dim = 2000),
               "140129" = panel.edges(left.crop = 0, upper.crop = 20, x.dim = 2000, y.dim = 2000))

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

# EXPLORATORY PLOTS: HISTOGRAMS ETC                                                             ####

hist(pw.m[,,"black", 1], breaks = "fd", ylim = c(0,30), xlim = c(0,65535))
hist(pw.m[,,"grey", 1], breaks = "fd", ylim = c(0,30), xlim = c(0,65535))
hist(pw.m[,,"white", 1], breaks = "fd", ylim = c(0,30), xlim = c(15000,25000))

# add thresholds
zz <- bad.px.limits(pw.m[,,"black", 1])
zz.J <- qJohnson(c(0.01, 0.99), JohnsonFit(pw.m[, , "black",1][!is.na(pw.m[, , "black", 1])]))
rect(0, 0, zz$dv, 50, col = adjustcolor("slateblue1", alpha = 0.3), border = NA)
rect(zz$dv, 0, zz$dm, 50, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
if (zz$ds > zz$dm) rect(zz$dm, 0, zz$ds, 50, col = adjustcolor("green3", alpha = 0.3), border = NA)
if (zz$bs < zz$bm) rect(zz$bm, 0, zz$bs, 50, col = adjustcolor("red", alpha = 0.3), border = NA)
rect(zz$bm, 0, zz$bv, 50, col = adjustcolor("orange", alpha = 0.3), border = NA)
rect(zz$bv, 0, 65535, 50, col = adjustcolor("gold", alpha = 0.3), border = NA)
rect(zz.J[1], 0, zz.J[2], 50, col = adjustcolor("blue", alpha = 0.3), border = NA)

pp <- which(pw.m[,,"white", 1] > 19000 & pw.m[,,"white", 1] < 21000, arr.ind = T)

plot(bp[[1]][,1:2], pch = 15, col = Cat.cols[bp[[1]]$type], cex = 0.5, asp = T)
points(pp, xlim = c(0,2000), ylim = c(0,2000), pch = 15, asp = T, col = "red", cex = 0.1)

bad.px.limits

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

# NON-RESPONSIVE LINES: WHAT HAPPENS TO THEM?                                                   ####

nr <- unique(unlist(lapply(bp, function(x) x$row[x$type %in% c("no response", "line.b", "line.d")])))

# per-column plots of all non-responsive/dim/bright lines
for (cc in nr) {
    pdf(paste0(fpath, "col-plot-", cc, ".pdf"), height = 3, width = 18); {
        par(mfrow = c(1, 6))
        for (i in 1:6) {
            plot(pw.m[cc,,"white", i], ylim = c(0,65535), col = "purple", type = "l", ylab = "", xlab = "")
            lines(pw.m[cc,,"grey", i], col = "blue")
            lines(pw.m[cc,,"black", i])
            points(1:2000, rep(0, 2000), pch = 15, cex = 0.7, col = c(NA, Cat.cols)[bpx2im(bp[[i]], im.dim = c(2000, 2000))[cc,]+1])
            abline(v = 1000.5, lty = 2, col = "red")
            text(0, 65000, names(bp)[i], pos = 4)
        }
        par(mfrow = c(1, 1))
        title(paste0("Column ", cc))
        dev.off()
    }
    }

# plot column artefacts in each image

qq <- lapply(bp, function(x) x[x$type %in% c("no response", "line.b", "line.d", "hot", "v.bright", "v.dim", "dim"),])
Cat.cols <- c("purple", "black", "magenta3", "red", NA, "gold", NA, NA, NA, "violet", NA, "green3", "green", NA, NA)

for (i in 1:5) {
    pdf(paste0(fpath, "extreme-px-", names(bp)[i], ".pdf")); {
        par(mar = c(2, 2, 1, 1))
        plot(qq[[i]][,1:2], pch = 15, col = Cat.cols[qq[[i]]$type], cex = 0.4, asp = T, 
             xlim = c(0,2000), ylim = c(0,2000))
#        draw.panels(p.dims[[i]])
        dev.off()
    }
}

# plot legend
{
    pdf(paste0(fpath, "plot-bpx-legend.pdf")); {
        plot.new()
        legend("center", col = Cat.cols, legend = Cat, pch = 15)
        dev.off()
        crop.pdf(paste0(fpath, "plot-bpx-legend.pdf"))
    }
}

# bright line: column 988
{
    plot(pw.m[988, 1:1000,"black", 1], type = "l", ylim = c(5000,6000))
    lines(pw.m[989, 1:1000,"black", 1], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[987, 1:1000,"black", 1], col = adjustcolor("orange", alpha = 0.4))
    
    plot(pw.m[988, 1:1000,"black", 2], type = "l", ylim = c(5000,6000))
    lines(pw.m[989, 1:1000,"black", 2], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[987, 1:1000,"black", 2], col = adjustcolor("orange", alpha = 0.4))
    
    plot(pw.m[988, 1:1000,"black", 3], type = "l", ylim = c(5000,6000))
    lines(pw.m[989, 1:1000,"black", 3], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[987, 1:1000,"black", 3], col = adjustcolor("orange", alpha = 0.4))
    
    plot(pw.m[988, 1:1000,"black", 4], type = "l", ylim = c(5000,6000))
    lines(pw.m[989, 1:1000,"black", 4], col = adjustcolor("cyan3", alpha = 0.4))
    lines(pw.m[987, 1:1000,"black", 4], col = adjustcolor("orange", alpha = 0.4))
    
    apply(pw.m[988, 1:1000, "black",], 2, mean, na.rm = T) - apply(pw.m[987, 1:1000, "black",], 2, mean, na.rm = T)
    apply(pw.m[988, 1:1000, "black",], 2, mean, na.rm = T) - apply(pw.m[989, 1:1000, "black",], 2, mean, na.rm = T)
}

dims <- lapply(lapply(bp, function(x) ddply(x[x$type %in% c("no response", "line.b", "line.d"),], .(row), summarise, 
              col.max = max(col), col.min = min(col),                       # top and bottom edge
              freq = length(row))),
              function(y) y[y$freq > 5,])

dim.change <- rmerge.df.list(head = unique(rbind.fill(lapply(dims, "[", 1))), lapply(dims, "[", c(1, 4)), by = 1, all = T)
# no change in dimensions of dim lines

pdf(paste0(fpath, "col-1680-progress.pdf"), width = 7, height = 4); {
    par(mar = c(2, 2, 1, 1))
    plot(pw.m[1680,,"white", 1], type = "l", xlim = c(0,1000), ylim = c(0, 60000), col = "purple")
    lines(pw.m[1680,,"grey", 1], col = "blue")
    lines(pw.m[1680,,"black", 1])
    
    for (i in 2:4) {
        lines(pw.m[1680,,"white", i], col = adjustcolor("violet", alpha = 0.2))
        lines(pw.m[1680,,"grey", i], col = adjustcolor("lightskyblue", alpha = 0.2))
    }
    abline(v = 1000.5, col = "red", lty = 2)
    legend("top", lty = 1, col = c("purple", "blue", "black", adjustcolor(c("violet", "lightskyblue"), alpha = 0.2)),
           legend = c("White 13-06-13", "Grey 13-06-13", "Black 13-06-13", "White (later dates)", "Grey (later dates)"), cex = 0.7)
    dev.off()
}

# transects of bad lines after refurbishment
# 512 & 1920: panel edges
{
    plot(pw.m[1019,,"white", 5], xlim = c(1000,2000), type = "l", col = "purple")
    lines(pw.m[1019,,"grey", 5], col = "blue")
    lines(pw.m[1019,,"black", 5])
    abline(v = 997.5, col = "red", lty = 2)
    
    pdf(paste0(fpath, "col-745-after-refurb.pdf"), width = 7, height = 4); {
        par(mar = c(2, 2, 1, 1))
        plot(pw.m[745,,"white", 5], xlim = c(1,1000), type = "l", col = "purple")
        lines(pw.m[745,,"grey", 5], col = "blue")
        lines(pw.m[745,,"black", 5])
        abline(v = 997.5, col = "red", lty = 2)
        dev.off()
    }

}

    
    pdf(paste0(fpath, "col-1374-before-refurb.pdf"), width = 7, height = 4); {
        par(mar = c(2, 2, 1, 1))
        plot(pw.m[1374,,"white", 1], xlim = c(1000, 2000), type = "l", col = "purple")
        lines(pw.m[1374,,"grey", 1], col = "blue")
        lines(pw.m[1374,,"black", 1])
        abline(v = 1000.5, col = "red", lty = 2)
        dev.off()
    }
    
}
####################################################################################################

# STATE TRANSITIONS                                                                             ####

tr <- list()

for (i in 1:(length(bp) - 1)) {
    tr[[i]] <- table("From" = ordered(c("normal", Cat)[c(bpx2im(bp[[i]], im.dim = c(2000, 2000))) + 1], levels = c("normal", Cat)),
                     "To" = ordered(c("normal", Cat)[c(bpx2im(bp[[i+1]], im.dim = c(2000, 2000))) + 1], levels = c("normal", Cat)))
    names(tr)[[i]] <- paste(names(bp)[c(i,(i+1))], collapse = "-")
}

tr <- array(unlist(tr), dim = c(dim(tr[[1]]), length(tr)), 
            dimnames = list(c("normal", Cat), c("normal", Cat), names(tr)))

# write csv of each transition in turn
for (i in 1:dim(tr)[3]) {
    write.csv(prep.csv(tr[,,i]), paste0(fpath, "trans-", dimnames(tr)[[3]][i], ".csv"), quote = F)
}

# plot each bad pixel map
for (i in 1:length(bp)) {
    pdf(paste0(fpath, "plot-bpx-", names(bp)[i], ".pdf")); {
        plot(bp[[i]][,1:2], pch = 15, cex = 0.3, asp = T, col = Cat.cols[bp[[i]]$type], xlab = "", ylab = "")
        draw.panels(p.dims[[i]])
        dev.off()
    }
}

# plot legend
{
    pdf(paste0(fpath, "plot-bpx-legend.pdf")); {
        plot.new()
        legend("center", col = Cat.cols, legend = Cat, pch = 15)
        dev.off()
    }
}

####################################################################################################

# OLD VS NEW CLASSIFICATIONS                                                                    ####

bp.old <- bp.old <- bp[[5]]
bp.new <- readRDS("./Other-data/bad-px-maps.rds")[[1]]

# yup, pretty sure they aren't related.
plot(bp.old[bp.old$type == "no response",1:2], col = "black", asp = T, pch = 20, xlim = c(0,2000), ylim = c(0,2000))
points(bp.new[bp.new$type == "no response", 1:2], col = "red")

plot(bp.old[bp.old$type == "hot",1:2], col = "black", asp = T, pch = 20, xlim = c(0,2000), ylim = c(0,2000))
points(bp.new[bp.new$type == "hot", 1:2], col = "red")

####################################################################################################

# STANDARD DEVIATION ACROSS EACH ACQUISITION                                                    ####

# plot SD (per subpanel?)

####################################################################################################

# PIXEL SHADING CORRECTION                                                                      ####

####################################################################################################