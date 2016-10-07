
dt <- "loan"
px <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))


####################################################################################################
offset <- c(24, 24)

dense.regions <- function(px, offset, th.u = 0.5, area = 11, th.l = 2/area, dilate.by = 21) {
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
        px$f.type <- NA
    }
    
    # define kernel
    k <-  matrix(rep(1 / area^2, area^2), ncol = area)
    
    # convert px to matrix & pad offset edges
    padded <- array(0, dim = apply(px[,1:2], 2, max))
    padded[as.matrix(fpx[,1:2])] <- 1
    
    padded[1:offset[1],] <- padded[offset[1] + 1,] 
    padded[,1:offset[2]] <- padded[,offset[2] + 1] 
    
    # convolve image with kernel
    px.density <- r2m(focal(m2r(padded), k))
    th.density <- (px.density > th.l) * 1
    
    # get clumps of mid-density pixels
    cc <- clump(m2r(th.density))
    
    if(nrow(cc) == 0) {return(px)}
    
    cand <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                       id = getValues(cc)[!is.na(getValues(cc))])
    cand$density <- px.density[as.matrix(cand[,1:2])]
    
    blocks <- ddply(cand, .(id), max.d = max(density), summarise)
    
    if ((dilate.by == 0) | is.na(dilate.by)) {
        qq <- as.matrix(cand[cand$id %in% blocks$id[blocks$max.d >= th.u], c("x", "y")])
    } else {
        th.density[as.matrix(cand[cand$id %in% blocks$id[blocks$max.d < th.u], c("x", "y")])] <- 0
        th.density[is.na(th.density)] <- 0
        
        sk <- shapeKernel(c(dilate.by, dilate.by), type = "disc")
        zz <- dilate(th.density, sk)
        
        qq <- which(zz == 1, arr.ind = T)
    }
    
    if(nrow(qq) == 0) {return(px)}
    
    # label all free pixels within affected regions
    npx <- merge(px, data.frame(qq, dense = T), by = c(1:2), all.x = T)
    npx$f.type[npx$dense & is.na(npx$f.type)] <- "dense.region"
    
    return(npx[, colnames(npx) %in% colnames(px)])
}

dp <- dense.regions(px[,1:3], offset = c(24, 24), th.u = 0.3)

pixel.plot(dp[dp$f.type == "dense.region",])
points(dp[is.na(dp$f.type),], col = "cyan3", pch = 15, cex = 0.4)

####################################################################################################

# TESTING KERNEL SIZE FOR COLUMN IDENTIFICATION                                                 ####
cc <- find.columns(px[,1:3])

table(cc$row[cc$f.type == "line.c"])


tst <- array(0, dim = c(19,19))
k3 <- round(matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
               rep((kernel.size - 1) / 3, kernel.size * 3), 
               rep(-1, kernel.size * floor(kernel.size / 2))),
             nrow = kernel.size),1)
k5 <- matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
               rep((kernel.size - 1) / 5, kernel.size * 5), 
               rep(-1, kernel.size * floor(kernel.size / 2))),
             nrow = kernel.size)

# single column - k works fine with threshold of 10                                 TH: 10
tst[11,] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 20                  x
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 7                   -
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 4                   -

# single column plus neighbour: no real change anywhere
tst[12,8] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 19                  x
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 8                   -
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 5                   -

# double column
tst[12,] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 15                  xx
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 13                  xx
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 8                   --

# double column with neighbour
tst[10,9] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 15                  xx
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 12, 15              xx
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 9                   --

# triple column
tst[10,] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 10                  xxx
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 8/20/8              -x-
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 12                  xxx

# single neighbour
tst[9,10] <- 1
round(r2m(focal(m2r(tst), k, pad = T, padValue = 0)),0)           # 10/9                --x
round(r2m(focal(m2r(tst), k3, pad = T, padValue = 0)),0)          # 10/19/7             xx-
round(r2m(focal(m2r(tst), k5, pad = T, padValue = 0)),0)          # 11/13               xxx
