
TO DO: JULY 14TH

 - LOOK FOR CHANGE-POINT APPROACH TO LINE DETECTION
 - LOOK FOR MEDIAN-DIFFERENCE CONVOLUTION-BASED APPROACH TO LINE DETECTION
 - CREATE PLOTS OF EACH TYPE OF LINE DEFECT FOR PAPER

####################################################################################################

# double columns in MCT225 (for example) are not picked up by current convolution approach
# check behaviour under convolution: are they always picked up as a line of width 5?

library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)
fpath <- "./Notes/MCT225/fig/"

####################################################################################################

# SIMPLE CONVOLUTION                                                                            ####

####################################################################################################

# CHECK ADJACENT COLUMNS                                                                        ####

# get columns of dark pixels by clustering (7 dark columns identified)
dark <- pw.m[,,"grey", "MCT225"]
dark[which(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] > 5000, arr.ind = T)] <- NA

cc <- clump(m2r(dark), dir = 4)
xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                 id = getValues(cc)[!is.na(getValues(cc))])
xy <- merge(xy, count(xy$id), by.x = "id", by.y = "x", all = T)
xy <- xy[xy$freq > 1,]

# check mean difference between adjacent column & its other neighbour
i <- 1
ll <- xy[xy$id == i & xy$x == 256, c("x", "y")]
mean(pw.m[,,"grey", "MCT225"][cbind(ll$x + 1, ll$y)] - pw.m[,,"grey", "MCT225"][cbind(ll$x + 2, ll$y)])
mean(pw.m[,,"grey", "MCT225"][cbind(ll$x + 2, ll$y)] - pw.m[,,"grey", "MCT225"][cbind(ll$x + 3, ll$y)])

mean(pw.m[,,"grey", "MCT225"][cbind(ll$x - 1, ll$y)] - pw.m[,,"grey", "MCT225"][cbind(ll$x - 2, ll$y)])
mean(pw.m[,,"grey", "MCT225"][cbind(ll$x - 2, ll$y)] - pw.m[,,"grey", "MCT225"][cbind(ll$x - 3, ll$y)])

# automate: iterate over each line identified
{
    im <- pw.m[,,"white", "MCT225"]
    line.px <- c()
    for (i in unique(xy$id)) {
        ll <- xy[xy$id == i, c("x", "y")]
        
        # find primary column
        pc <- ll[ll$x == count(ll$x)$x[which.max(count(ll$x)$freq)],]
        
        # get neighbouring columnwise mean differences at each offset
        col.diffs <- sapply(c(-2, -1, 1, 2), 
               function(os) median(im[cbind(pc$x + os, pc$y)] - im[cbind(pc$x + os + sign(os), pc$y)]))

        # identify secondary column 
        ll <- rbind(ll, data.frame(x = pc$x + c(-2, -1, 1, 2)[which.max(abs(col.diffs))],
                               y = pc$y))
        
        line.px <- rbind(line.px, ll)
        # could also update this to include ANY columns with diff above certain threshold 
        # should probably also output columnwise difference somewhere for manual check
    }
}

plot(line.px, pch = 15, xlim = 256 + c(-20,20), ylim = 1880 + c(-20,20))

# check dark line ends
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(256, 1880), bad.px = line.px, bpx.cex = 5)
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(448, 1137), bad.px = line.px, bpx.cex = 5)
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(882, 1119), bad.px = line.px, bpx.cex = 5)
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(1691, 1359), bad.px = line.px, bpx.cex = 5)
focal.plot(pw.m[,,"grey", "MCT225"], centre = c(1926, 1627), bad.px = line.px, bpx.cex = 5)

plot(pw.m[1688,,"white", "MCT225"], col = "green3", xlim = c(1360, 2048), ylim = c(45000, 60000), type = "l")
lines(pw.m[1689,,"white", "MCT225"], col = "forestgreen")
lines(pw.m[1690,,"white", "MCT225"], col = "gold")
lines(pw.m[1691,,"white", "MCT225"], col = "blue")
lines(pw.m[1692,,"white", "MCT225"], col = "red")

####################################################################################################

####################################################################################################

# REPLACE DARK PIXELS WITH LOCAL MEDIAN                                                         ####

# get columns of dark pixels by clustering (7 dark columns identified)
dark <- pw.m[,,"grey", "MCT225"]
dark[which(pw.m[,,"grey", "MCT225"] - pw.m[,,"black", "MCT225"] > 5000, arr.ind = T)] <- NA

cc <- clump(m2r(dark), dir = 4)
xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                 id = getValues(cc)[!is.na(getValues(cc))])
xy <- merge(xy, count(xy$id), by.x = "id", by.y = "x", all = T)
xy <- xy[xy$freq > 1,]

med.replace <- function(im, px, w = 5) {
    get.px <- cbind(px[1] + c(-floor(w/2):floor(w/2)), px[2])
    median(im[get.px], na.rm = T)
}

mr <- apply(xy[,c("x", "y")], 1, med.replace, im = pw.m[,,"grey", "MCT225"])

org <- pw.m[,,"grey", "MCT225"][as.matrix(xy[,c("x", "y")])]

{
    plot(mr[xy$id == 1], org[xy$id == 1], pch = 20, xlim = c(18000, 21000), ylim = c(4500, 7500))
    points(mr[xy$id == 2], org[xy$id == 2], pch = 20, col = "green3")
    points(mr[xy$id == 3], org[xy$id == 3], pch = 20, col = "red")
    points(mr[xy$id == 4], org[xy$id == 4], pch = 20, col = "blue")
    points(mr[xy$id == 5], org[xy$id == 5], pch = 20, col = "gold")
    points(mr[xy$id == 6], org[xy$id == 6], pch = 20, col = "pink")
    points(mr[xy$id == 7], org[xy$id == 7], pch = 20, col = "cyan3")
}

grey.adj <- pw.m[,,"grey", "MCT225"]
grey.adj[as.matrix(xy[,c("x", "y")])] <- mr

focal.plot(pw.m[,,"grey", "MCT225"], centre = c(448, 1140), surround = 7)
focal.plot(grey.adj, centre = c(448, 1140), surround = 7)

# CONVOLUTION
convolve.vertical <- function(im, k.size = 5) {
    
    # define kernel to use in convolution
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    # perform convolution
    r2m(focal(m2r(im), k))
}

conv <- convolve.vertical(grey.adj)
focal.plot(conv, centre = c(448, 1140), surround = 7)
pixel.image(conv)

# THRESHOLDING
hist(conv, breaks = "fd", ylim = c(0,30))
abline(v = 5500 * c(-1,1), col = "red", lty = 2)

th <- array(findInterval(conv, c(floor(min(conv, na.rm = T)), 5000 * c(-1,1), ceiling(max(conv, na.rm = T)))) - 2,
            dim = dim(conv))
table(th)
plot(which(th < 0, arr.ind = T), col = "blue", pch = ".", xlim = c(0,2048), ylim = c(0,2048))
points(which(th > 0, arr.ind = T), col = "orange", pch = ".")
points(xy[,c("x", "y")], pch = ".", col = "black")

# INTERIM FILTERING? - remove any short segment with no vertical neighbours
cc.int <- clump(m2r(th), dir = 4)
xy.int <- data.frame(xyFromCell(cc.int, which(!is.na(getValues(cc.int)))),
                 id = getValues(cc.int)[!is.na(getValues(cc.int))])

xy.int <- merge(xy.int, count(xy.int$id), by.x = "id", by.y = "x", all = T)
xy.int <- xy.int[xy.int$freq > 1,]

int.summ <- ddply(xy.int, .(id, x), summarise,
                  ymin = min(y), ymax = max(y), len = length(x))
table(int.summ$len)

# SMOOTHING & REPAIR

# FILTER SHORT SEGMENTS



# REPEAT FOR BRIGHT SEGMENTS


####################################################################################################

####################################################################################################
    
# PROCESS INTO COLUMNWISE DIFFERENCES FIRST                                                     ####

diff1 <- pw.m[1:2046,,"black", "131122"] - pw.m[2:2047,,"black", "131122"]
diff2 <- pw.m[3:2048,,"black", "131122"] - pw.m[2:2047,,"black", "131122"]

par(mfrow = c(2, 1))
pixel.image(diff1, xlim = c(950, 1050), ylim = c(300, 400))
pixel.image(diff2, xlim = c(950, 1050), ylim = c(300, 400))
par(mfrow = c(1, 1))

o.plot(diff1[1000:1020, 320])
lines(diff2[1000:1020, 320], col = "red")

diff <- apply(abind(diff1, diff2, along = 3), 1:2, mean, na.rm = T)

hist(diff, breaks = "fd", ylim = c(0,30))
lines(diff[1000:1020, 320], col = "blue")

pixel.image(diff)

dl <- find.lines(diff)
table(dl)
points(which(dl > 0, arr.ind = T), pch = ".")
draw.panels(col = "blue")
dld <- find.lines(diff, dim.lines = T)
points(which(dld > 0, arr.ind = T), pch = ".", col = "red")

pixel.image(diff, xlim = c(800,1200), ylim = c(900,1100))

# this approach gives a head-start on convolution, so thresholds will differ
conv <- convolve.lines(diff)
pixel.image(conv, xlim = c(800,1200), ylim = c(900,1100))
draw.panels()

o.plot(conv[,900], xlim = c(950, 1050))
abline(h = 5500 * c(-1,1), col = "red")

# if 300px higher than neighbours...
test <- array(0, dim = c(55,55))
test[23,] <- 300
o.plot(test[,30], ylim = c(-300,300))

t.diff <- array(dim = dim(test))
t.diff[2:54,] <- apply(abind(test[2:54,] - test[1:53,],
                             test[2:54,] - test[3:55,], along = 3), 1:2, mean, na.rm = T)
o.plot(t.diff[,30], add = T, col = "blue")

td.conv <- convolve.lines(t.diff)
pixel.image(td.conv)

o.plot(td.conv[,30], col = "red", add = T)

pw.m[1:2046,,"black", "131122"] - pw.m[2:2047,,"black", "131122"]
diff2 <- pw.m[3:2048,,"black", "131122"] - pw.m[2:2047,,"black", "131122"]


#------------------------------------------------------------
# LET'S DO THIS A LITTLE MORE SYSTEMATICALLY

# line segment of +300 will convolve to 7500. Threshold at 7000 to avoid picking up 
test.convolution <- function(mat) {
    
    nc <- ncol(mat)
    td <- array(dim = dim(mat))
    td[2:(nc-1),] <- apply(abind(mat[2:(nc-1),] - mat[1:(nc-2),],
                                       mat[2:(nc-1),] - mat[3:nc,], along = 3), 1:2, mean, na.rm = T)
    tdc <- convolve.lines(td)
    
    plot(mat[,23], type = "l", ylim = range(tdc, na.rm = T))
    lines(td[,23], col = "blue")
    lines(tdc[,23], col = "red")
    abline(h = 7000, col = "green3")
    
    text(23, tdc[23,23], tdc[23,23], pos = 4)
    
    return(abind(diff = td, conv = tdc, along = 3))
}

#------------------------------------------------------------
# single line +300: 7500
{
    test <- array(0, dim = c(55,55))
    test[23,] <- 300
    
    tmp <- test.convolution(test)
}

#------------------------------------------------------------
# double line +300: 3750
{
    test <- array(0, dim = c(55,55))
    test[23:24,] <- 300
    
    tmp <- test.convolution(test)
}

#------------------------------------------------------------
# double column (one column +300, one +600): main column 11250, secondary 0
{
    test <- array(0, dim = c(55,55))
    test[23,] <- 600
    test[24,] <- 300
    
    tmp <- test.convolution(test)
    tmp[24,23,2]
}

#------------------------------------------------------------
# double column (one column +300, one +500): main column 8750, secondary 1259
{
    test <- array(0, dim = c(55,55))
    test[23,] <- 500
    test[24,] <- 300
    
    tmp <- test.convolution(test)
    tmp[23:24,23,2]
}

#------------------------------------------------------------
# double column (dark column removed): main column NA, secondary NA
{
    test <- array(0, dim = c(55,55))
    test[23,] <- NA
    test[24,] <- 300
    
    tmp <- test.convolution(test)
    tmp[23:24,23,2]
}

#------------------------------------------------------------
# single point +300: 1500
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 300
    
    tmp <- test.convolution(test)
}

#------------------------------------------------------------
# two points +300: 3000
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 300
    test[23,24] <- 300
    
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# three points +300: 4500
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 300
    test[23,24] <- 300
    test[23,22] <- 300
    
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# four points +300: 6000
{
    test <- array(0, dim = c(55,55))
    test[23,22:25] <- 300
    
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# five points +300: 7500 (essentially a line segment)
{
    test <- array(0, dim = c(55,55))
    test[23,21:25] <- 300
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# six points +300: 7500 (essentially a line segment)
{
    test <- array(0, dim = c(55,55))
    test[23,21:26] <- 300
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# three points, corner +300: 2250
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 300
    test[23,24] <- 300
    test[24,23] <- 300
    
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# four points, square +300: spreads effect, 2 columns at 1500
{
    test <- array(0, dim = c(55,55))
    test[23:24,23:24] <- 300
    
    tmp <- test.convolution(test)
}
#------------------------------------------------------------
# panel edge +300 on one side only: 4375
{
    test <- array(0, dim = c(55,55))
    test[23,] <- 300
    test[24:55,] <- 250
    
    tmp <- test.convolution(test)
}

#------------------------------------------------------------
# SMEARING ISSUES

# single point +1500: 5px > 7500
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 1500
    
    tmp <- test.convolution(test)
    pixel.image(tmp[,,2])
    o.plot(tmp[23,,2])
}

# two points +1500, sep = 1: 7px > 7000
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 1500
    test[23,25] <- 1500
    
    tmp <- test.convolution(test)
    o.plot(tmp[23,,2])
}

# two points +1500, sep = 2: 8px > 7000
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 1500
    test[23,26] <- 1500
    
    tmp <- test.convolution(test)
    o.plot(tmp[23,,2])
    sum(tmp[,,2] > 7000, na.rm = T)
}

# two points +1500, sep = 3: 9px > 7000
{
    test <- array(0, dim = c(55,55))
    test[23,23] <- 1500
    test[23,27] <- 1500
    
    tmp <- test.convolution(test)
    o.plot(tmp[23,,2])
    sum(tmp[,,2] > 7000, na.rm = T)
    
}

# two points +1500, sep = 4: 10px > 7000
{
    test <- array(0, dim = c(55,55))
    test[23,22] <- 1500
    test[23,27] <- 1500
    
    tmp <- test.convolution(test)
    o.plot(tmp[23,,2])
    sum(tmp[,,2] > 7000, na.rm = T)
    
    pixel.image(tmp[,,2])
}

# two points +1500, sep = 5: two 5px segments @ 7500.
{
    test <- array(0, dim = c(55,55))
    test[23,21] <- 1500
    test[23,27] <- 1500
    
    tmp <- test.convolution(test)
    o.plot(tmp[23,,2])
    sum(tmp[,,2] > 7000, na.rm = T)
    
    pixel.image(tmp[,,2])
}

####################################################################################################

# TRY DIFFERENCED APPROACH ON REAL DATA                                                         ####

# 160430
{
    diff.160430 <- rbind(NA, apply(abind(pw.m[2:2047,,"black", "160430"] - pw.m[1:2046,,"black", "160430"],
                                         pw.m[2:2047,,"black", "160430"] - pw.m[3:2048,,"black", "160430"],
                                         along = 3), 1:2, mean, na.rm = T), NA)
    
    hist(diff.160430, breaks = "fd", ylim = c(0,30))
    abline(v = 7500 * c(-1,1), lty = 2, col = "red")
    
    o.plot(diff.160430[429,], xlim = c())
    
    diff.lines.160430 <- find.lines(diff.160430, threshold.at = 7000)
    table(diff.lines.160430)
}

# 131122: doesn't make much difference :(
{
    diff.131122 <- rbind(NA, apply(abind(pw.m[2:2047,,"grey", "131122"] - pw.m[1:2046,,"grey", "131122"],
                                         pw.m[2:2047,,"grey", "131122"] - pw.m[3:2048,,"grey", "131122"],
                                         along = 3), 1:2, mean, na.rm = T), NA)

    diff.lines.131122 <- find.lines(diff.131122, threshold.at = 7000)
    table(diff.lines.131122)
    
    diff.lines.131122.2 <- find.lines(diff.131122, threshold.at = 7000, dim.lines = T)
    table(diff.lines.131122.2)
    
    plot(which(diff.lines.131122 > 0, arr.ind = T), pch = 15, xlim = c(0,2048), ylim = c(0,2048))
    points(which(diff.lines.131122.2 > 0, arr.ind = T), pch = 15, col = "blue")
    draw.panels()
    
    c.lines.131122 <- find.lines(pw.m[,,"grey", "131122"], dim.lines = T)
    points(which(c.lines.131122 > 0, arr.ind = T), pch = 15, col = "red")
}

# convolved with & without differencing:
{
    c.160430 <- convolve.lines(pw.m[,,"black", "160430"])
    cd.160430 <- convolve.lines(diff.160430)
    
    difference <- cd.160430 - c.160430
    
    hist(difference, breaks = "fd", ylim = c(0,30))
    pixel.image(difference)
    
    plot(c.160430[429,], xlim = c(1024,2048), ylim = c(-1000,20000), type = "l")
    lines(cd.160430[429,], col = "blue")
    
    plot(c.160430[811,], xlim = c(0,1024), ylim = c(-1000,20000), type = "l")
    lines(cd.160430[811,], col = "blue")
    
    plot(c.160430[800,], xlim = c(0,1024), ylim = c(-5000,20000), type = "l")
    lines(cd.160430[800,], col = "blue")
    mean(cd.160430[800,1:1024], na.rm = T)
    mean(c.160430[800,1:1024], na.rm = T)
}

####################################################################################################

# ALREADY TRIED                                                                                 ####


# try convolution for lines after setting all dark pixels to NA
active <- pw.m; {
    active[,,,"MCT225"] <-  array(apply(pw.m[,,,"MCT225"], 3, 
                                        function(im) {
                                            im[as.matrix(bpx.MCT225[,1:2])] <- NA
                                            return(im)
                                        }),
                                  dim = dim(pw.m[,,,"MCT225"]), dimnames = dimnames(pw.m[,,,"MCT225"]))
    active[,,,"160430"] <-  array(apply(pw.m[,,,"160430"], 3, 
                                        function(im) {
                                            im[as.matrix(bpx.160430[,1:2])] <- NA
                                            return(im)
                                        }),
                                  dim = dim(pw.m[,,,"160430"]), dimnames = dimnames(pw.m[,,,"160430"]))                                       
    active[,,,"131122"] <-  array(apply(pw.m[,,,"131122"], 3, 
                                        function(im) {
                                            im[as.matrix(bpx.131122[,1:2])] <- NA
                                            return(im)
                                        }),
                                  dim = dim(pw.m[,,,"131122"]), dimnames = dimnames(pw.m[,,,"131122"])) 
}

zz <- find.lines(pw.m[,,"black", "160430"])
ll.160430 <- data.frame(which(zz > 0, arr.ind = T), id = zz[which(zz > 0)])
plot(ll.160430[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0,2048),
     col = c( "red", "green3", "cyan3", "blue", "purple")[ll.160430$id])

zz <- find.lines(pw.m[,,"grey", "MCT225"])
ll.MCT225 <- data.frame(which(zz > 0, arr.ind = T), id = zz[which(zz > 0)])
plot(ll.MCT225[,1:2], pch = ".", xlim = c(0,2048), ylim = c(0,2048),
     col = c( "red", "green3", "cyan3", "blue", "purple")[ll.160430$id])

pixel.image(zz, xlim = 256 + c(-20, 20), ylim = 1880 + c(-20, 20))


# standard algorithm finds nothing in MCT225. Needs revision
# can't set to NA, as this will smear & block neighbouring column
# try replacing with local median before convolution

med.replace <- function(im, px, w = 5) {
    get.px <- cbind(px[1] + c(-floor(w/2):floor(w/2)), px[2])
    median(im[get.px], na.rm = T)
}

# convolution of MCT225 with linear kernel
k.size <- 5
k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
conv <- r2m(focal(m2r(active[,,"grey", "MCT225"]), k, na.rm = T))

pixel.image(conv)

pixel.image(active[,,"grey", "MCT225"], xlim = 256 + c(-20, 20), ylim = 1880 + c(-20, 20))

le.plot(256, 1880)
le.plot(448, 1137)
le.plot(882, 1119)
le.plot(1691, 1359)
le.plot(1926, 1627)

hist(conv, breaks = "fd", ylim = c(0,30))

# subtract whole-column offset from each pixel first
c.adj <- sweep(pw.m[,,"grey", "MCT225"], 2, colSums(pw.m[,,"black", "MCT225"], na.rm = T) / 2048, "-")

pixel.image(c.adj, xlim = 256 + c(-20, 20), ylim = 1880 + c(-20, 20))

####################################################################################################

# DIFFERENCING + MA                                                                             ####

col.diffs <- function(im) {
    
    nc <- ncol(im)
    diff1 <- im[1:(nc-2),] - im[2:(nc-1),]
    diff2 <- im[3:(nc),] - im[2:(nc-1),]
    
    rbind(NA, apply(abind(diff1, diff2, along = 3), 1:2, mean, na.rm = T), NA)
}

cd <- col.diffs(fixed[,,"black", "131122"])

plot(cd[232,], xlim = c(1025, 2048), type = "l")
lines(cd[231,], col = "blue")

# what about median-smoothed values
md <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                   function(nm) readRDS(paste0("./02_Objects/med-diffs/md-", nm, ".rds")), 
                   simplify = F),
            along = 4)

plot(md[232,,"black", "131122"], type = "l", xlim = c(1024, 2048))

sd(md[,,"black", "131122"], na.rm = T)

