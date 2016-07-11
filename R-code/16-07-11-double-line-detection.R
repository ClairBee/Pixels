
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