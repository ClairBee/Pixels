
library("IO.Pixels"); library("CB.Misc")

# LINE DETECTION ALGORITHM UNSATISFACTORY

# final code should use single image unless necessary to do otherwise - easier to replicate
im <- readRDS("./02_Objects/images/pwm-160430.rds")             # pixelwise mean images
md7 <- readRDS("./02_Objects/med-diffs/md7-160430.rds")         # residual after median-smoothing

####################################################################################################

# identify lines of dark pixels & replace with smoothed values
# identify single lines
# identify double lines
# identify split lines (either side of a removed dark line)

# use MCT225 & 131122 as exemplars. Also loan panel?
unique.panels <- c("130613", "140128", "160430", "loan", "MCT225")

pw.m <- load.objects("./02_Objects/images/", otype = "pwm", unique.panels)
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7", unique.panels)

####################################################################################################

# IDENTIFY DARK LINES                                                                           ####

# check histogram of all white images - is there an obvious break?
# small peak between 30000 & 35000 is single dim acquisition 130701
# tiny peak visible at 6000 is collection of dark lines
# reasonable threshold is 25000 
hist(pw.m[,,"white",], breaks = "fd")


# very little difference between using 15000 and 25000 on this data set, so capture maximum # points
ll <- data.frame(d.10 = apply(pw.m[,,"white",], 3, function(im) length(which(im < 10000))),
                 d.15 = apply(pw.m[,,"white",], 3, function(im) length(which(im < 15000))),
                 d.25 = apply(pw.m[,,"white",], 3, function(im) length(which(im < 25000))))

# identify dark columns (& dark rows, if there were any)
cc <- dark.lines(im[,,"white"])
rr <- dark.lines(t(im[,,"white"]))

# display dark lines identified
pixel.plot(cc, cex = 0.4); points(rr, pch = 15, cex = 0.4, col = "red")


####################################################################################################


# REPLACE DARK LINES                                                                            ####

# replacing with NA would unbalance the convolution
# replacing dark pixels with 'expected' median-smoothed value would give 0 residual
md7.adj <- array(apply(md7, 3, function(vv) replace(vv, as.matrix(cc[,1:2]), 0)), dim = dim(md7), dimnames = dimnames(md7))

# display resulting smoothed images
pixel.image(md7.adj[,,"black"], title = "Black image, dark lines removed")
pixel.image(md7.adj[,,"grey"], title = "Grey image, dark lines removed")
pixel.image(md7.adj[,,"white"], title = "White image, dark lines removed")

####################################################################################################

# IDENTIFY BRIGHT / DIM LINES IN BLACK IMAGES                                                   ####

# identify lines in black images first because it is far less noisy than the grey image

# STEP-BY-STEP PROCEDURE (for reference)

# convolution with linear kernel
cl <- convolve.lines(md7.adj[,,"black"])
pixel.image(cl)
{
    # cropped histogram of convolved values shows peaks corresponding to identified columns
    hist(cl, breaks = "fd", ylim = c(0,100))
    abline(v = c(-5500, 5500), col = "red")
    
    sum(cl > 5500, na.rm = T)               # 1521 bright pixels identified
    sum(cl < - 5500, na.rm = T)             # 2086 dim pixels identified
    
    # columns are clearly visible when identified pixels are plotted
    pixel.plot(which(cl > 5500, arr.ind = T), col = "red")
    points(which(cl < -5500, arr.ind = T), pch = 15, cex = 0.4)
}

# threshold resulting convolution
th <- cl > 5500
th[which(cl < -5500, arr.ind = T)] <- 2


# identify & summarise each linear segment
# a single pixel will 'smear' to 5px segment
# a cluster of 4 very bright pixels will 'smear' to 8px
# 2 very bright single pixels, separated by 5px -> line of 10px
# set 10px as minimum length for acceptance as possible line segment


min.length <- 10
lines <- clump(m2r(th), dir = 4)
xy <- data.frame(xyFromCell(lines, which(!is.na(getValues(lines)))), 
                 id = getValues(lines)[!is.na(getValues(lines))])
xy$type <- th[as.matrix(xy[,1:2])]
xy$adj.value <- md7.adj[,,"black"][as.matrix(xy[,1:2])]
xy$conv.value <- cl[as.matrix(xy[,1:2])]

dl <- ddply(xy, .(id, x, type), summarise, ymin = min(y), 
            ymax = max(y), length = length(x), density = length / (ymax - ymin + 1),
            med.gv = median(abs(adj.value)), min.gv = min(abs(adj.value)), max.gx = max(abs(adj.value)),
            med.conv = median(conv.value))
dl <- dl[dl$length > min.length & dl$density > 0.5, ]

l.px <- cbind(x = unlist(apply(dl[, ], 1, function(ll) rep(ll["x"], ll["length"]))), 
              y = unlist(apply(dl[, ], 1, function(ll) seq(ll["ymin"], ll["ymax"]))))

####################################################################################################

# FIND SINGLE MOST SIGNIFICANT CHANGEPOINT IN EACH CANDIDATE

# changepoints are identified in each column to determine exactly how line should be defined
library(changepoint)

cp1 <- cpt.mean(md7.adj[232,1025:2048, "black"][!is.na(md7.adj[232,1025:2048, "black"])])
o.plot(md7.adj[232,1025:2048, "black"])
abline(v = cp1@cpts, col = "red")
abline(h = cp1@param.est$mean, col = "blue")
# changepoint found at 248

cp2 <- cpt.mean(md7.adj[1396,1025:2048, "black"][!is.na(md7.adj[1396,1025:2048, "black"])])
o.plot(md7.adj[1396,1025:2048, "black"])
abline(v = cp2@cpts, col = "red")
abline(h = cp2@param.est$mean, col = "blue")
# changepoint found at 71

cp3 <- cpt.mean(md7.adj[1012, 1024:1, "black"][!is.na(md7.adj[1012, 1024:1, "black"])])
o.plot(md7.adj[1012,1024:1, "black"])
abline(v = cp3@cpts, col = "red")
abline(h = cp3@param.est$mean, col = "blue")
# no changepoints found - whole column

####################################################################################################

# REMOVE LINES FROM GREY RESIDUALS, RETEST                                                      ####

n.plot(md7.adj[,,"black"], 1012, xlim = c(0,1024))
n.plot(md7.adj[,,"black"], 1396, xlim = c(1025, 2048))

# are there any lines that appear only in the white/black images?

####################################################################################################

# FUNCTION CONSTRUCTION AREA                                                                    ####

# to be run on median-differenced image with dark columns set to 0 value

ol <- offset.lines(md7.adj[,,1])
ol.r <- offset.lines(t(md7.adj[,,1]), midline = NA, min.length = 20)

o.plot(im[,1025,"white"] - im[,1024, "white"], ylim = c(0,1000))
o.plot(md[1026,])

####################################################################################################

# RUN OVER ALL IMAGES                                                                           ####

# identify dark lines in white images
dl <- apply(pw.m[,,"white", ], 3, dark.lines)

# replace dark lines with 0 in residual dark images
md7.adj <- abind(sapply(names(dl), 
                  function(dt) {
                      b <- replace(md7[,,"black", dt], as.matrix(dl[[dt]][,1:2]), 0)
                      g <- replace(md7[,,"grey", dt], as.matrix(dl[[dt]][,1:2]), 0)
                      w <- replace(md7[,,"white", dt], as.matrix(dl[[dt]][,1:2]), 0)
                      abind(b, g, w, along = 3, new.names = dimnames(md7[,,,dt]))
                      }, simplify = F),
                 along = 4, new.names = dimnames(md7))


# identify offset lines in adjusted dark residuals
ol.c <- apply(md7.adj, 3, offset.lines, min.length = 15)
ol.r <- apply(md7.adj, 3, function(im) offset.lines(t(im), min.length = 15))

# plots of results
{
    pixel.plot(dl$"130613")
    points(ol.c$"130613"[ol.c$"130613"$"type" == "dim",1:2], pch = 15, cex = 0.4, col = "blue")
    points(ol.c$"130613"[ol.c$"130613"$"type" == "bright",1:2], pch = 15, cex = 0.4, col = "red")
    
    table(ol.c$"130613"[,c("x", "type")])
    o.plot(pw.m[232,,"black", "130613"] - pw.m[231,,"black", "130613"], xlim = c(1025, 2048)); abline(h = 0, col = "red", lty = 2)
    o.plot(pw.m[1012,,"black", "130613"] - pw.m[1011,,"black", "130613"], xlim = c(1, 1024)); abline(h = 0, col = "red", lty = 2)
    o.plot(pw.m[1396,,"black", "130613"] - pw.m[1395,,"black", "130613"], xlim = c(1025, 2048)); abline(h = 0, col = "red", lty = 2)
}


pixel.plot(ol.c$"160430")

####################################################################################################

# COINCIDENT NEIGHBOUR-DIFFERENCES                                                              ####

d1 <- cbind(NA, im[,1:2046,"black"] - im[,2:2047,"black"], NA)
d2 <- cbind(NA, im[,3:2048,"black"] - im[,2:2047,"black"], NA)

px <- which(abs(d1) > 300 | abs(d2) > 300, arr.ind = T)

pixel.plot(px, ylim = c(1000,1100), xlim = c(100,200))
abline(h = 1025, col = "red")
points(which(abs(d1) > 300, arr.ind = T), col = "red")      # identifies correct line
points(which(abs(d2) > 300, arr.ind = T), col = "blue")
points(which(abs(d2) > 300 & abs(d1) > 300, arr.ind = T), col = "cyan2", pch = 25, cex = 0.5)

pixel.plot(px, ylim = c(0,100), xlim = c(100,200))
points(which(abs(d1) > 300, arr.ind = T), col = "red")      # identifies correct line
points(which(abs(d2) > 300, arr.ind = T), col = "blue")
points(which(abs(d2) > 300 & abs(d1) > 300, arr.ind = T), col = "cyan2", pch = 15, cex = 0.5)   # appears in both
abline(h = 77, col = "red")

# First preference: get coordinates from both
# Second preference: get coordinates from only d1
o.plot(im[1024,,"black"])

