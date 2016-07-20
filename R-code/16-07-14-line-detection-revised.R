
# NEEDS MORE WORK ON CLASSIFICATION

# could fit a changepoint model manually (see 'refs' bookmarks for details)?

####################################################################################################

# REVISED IDENTIFICATION OF COLUMN DEFECTS

library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)

md7 <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                    function(nm) readRDS(paste0("./02_Objects/med-diffs/md7-", nm, ".rds")), 
                    simplify = F),
             along = 4)

####################################################################################################

# VERBOSE PROCEDURE                                                                             ####

overplot <- function(im.array, column, dt, xlim = c(0, 2048), hline = 0, vline = 1024.5, ...) {
    dt <- toString(dt)
    par(mar = c(2,2,2,1))
    plot(im.array[column, , "white", dt], type = "l", xlim = xlim, col = "green3", 
         main = paste0(dt, " column ", column, " (", as.list(sys.call())[[2]], ")"), ...)
    lines(im.array[column, , "grey", dt], col = "blue")
    lines(im.array[column, , "black", dt])
    
    if (!is.na(vline)) abline(v = vline, col = "red", lty = 2)
    if (!is.na(hline)) abline(h = hline)
}

{
    overplot(md7, 429, "160430", xlim = c(1024, 2048), ylim = c(-1000, 2000))
    overplot(md7, 745, "140128", xlim = c(0,1024), ylim = c(-1000, 2000))
    
    overplot(md7, 736, "131122", xlim = c(0,1024))
    
    overplot(md7, 750, "131122", xlim = c(0,1024))
    
    
    overplot(md7, 447, "MCT225", xlim = c(1024, 2048))
    overplot(md7, 448, "MCT225", xlim = c(1024, 2048))
    overplot(md7, 449, "MCT225", xlim = c(1024, 2048), ylim = c(-50000, -40000))
    overplot(md7, 450, "MCT225", xlim = c(1024, 2048))
}


#----------------------------------------------------------------------------
# ROUTE 0: REPLACE DARK PIXELS WITH MEDIAN VALUES, THEN CONVOLVE

# taking a new median leads to over-smoothing (variance change), so can only check for mean changepoint
{
    med.replace <- function(im, px, w = 5) {
        
        get.px <- cbind(px[1] + c(-floor(w/2):floor(w/2)), px[2])
        get.px <- get.px[get.px[,1] %in% c(1:2048),]
        median(im[get.px], na.rm = T)
    }
    
    mr7 <- md7          # array to hold median-replaced values
    
    for (dt in dimnames(pw.m)[[4]]) {
        for (icol in dimnames(pw.m)[[3]]) {
            px <- which(pw.m[,,"white", dt] - pw.m[,,"black", dt] < 10000 & pw.m[,,"white", dt] < 15000, arr.ind = T)
            mr7[,,icol,dt][px] <- apply(px, 1, med.replace, im = mr7[,,icol, dt], w = 11)
        }
    }
    
    plot(md7[1926, , "grey", "MCT225"], type = "l")
    lines(mr7[1926,,"grey", "MCT225"], col = "blue")
    
    plot(md7[1927, , "grey", "MCT225"], type = "l")
    lines(mr7[1927, , "grey", "MCT225"], col = "red")
    

}

# convolution (kernel splitting now unnecessary)
{
    k.size <- 5
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    conv.mr <- array(apply(mr7, 3:4, function(im) r2m(focal(m2r(im), k))), dim = dim(md7), dimnames = dimnames(md7))
    {
        overplot(conv.mr, 429, "160430", xlim = c(1024, 2048), hline = 3000)
        overplot(conv.mr, 745, "140128", xlim = c(0,1024), hline = 3000)    # smoothed away
        
        overplot(conv.mr, 736, "131122", xlim = c(0,1024), hline = 3000)
        
        overplot(conv.mr, 447, "MCT225", xlim = c(1024, 2048))
        overplot(conv.mr, 448, "MCT225", xlim = c(1024, 2048))
        overplot(conv.mr, 449, "MCT225", xlim = c(1024, 2048))
        overplot(conv.mr, 450, "MCT225", xlim = c(1024, 2048))  
        
        overplot(md7, 447, "MCT225", xlim = c(1024, 2048))
        overplot(md7, 448, "MCT225", xlim = c(1024, 2048))
        overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
        overplot(md7, 450, "MCT225", xlim = c(1024, 2048))  
        }
    
}


#----------------------------------------------------------------------------
# ROUTE 1: CONVOLUTION WITH EDGE DETECTION KERNEL
# convolution generally doesn't work as well with thicker edges though.
{
    # convolution of median-differenced images with 5-square kernel
    k.size <- 5
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    conv <- array(apply(md7, 3:4, function(im) r2m(focal(m2r(im), k))), dim = dim(md7), dimnames = dimnames(md7))
    {
        overplot(conv, 429, "160430", xlim = c(1024, 2048), hline = 3000)
        overplot(conv, 745, "140128", xlim = c(0,1024), hline = 3000)
        
        overplot(conv, 736, "131122", xlim = c(0,1024), hline = 3000)
        
        overplot(conv, 447, "MCT225", hline = 3000)
        overplot(conv, 448, "MCT225", hline = 3000)
        overplot(conv, 449, "MCT225", hline = 3000)
        overplot(conv, 450, "MCT225", hline = 3000)
    }
    
    # thresholding (using low threshold of 3000)
    th <- array(apply(conv, 3:4, function(im) abs(im) > 3000),
                dim = dim(conv), dimnames = dimnames(conv))
    {
        overplot(th, 429, "160430", xlim = c(1024, 2048)
        
        overplot(th, 745, "140128", xlim = c(0,1024))
        
        overplot(th, 736, "131122", xlim = c(0,1024))
        
        overplot(th, 447, "MCT225")
        overplot(th, 448, "MCT225")
        overplot(th, 449, "MCT225")
        overplot(th, 450, "MCT225")
    }
    
    # summarise candidate columns
    tth <- th[,,"black", "160430"]
    
    xt <- ddply(data.frame(which(abs(tth) > 0, arr.ind = T)), .(x = row), summarise,
                ymin = min(col), ymax = max(col), range = ymax - ymin + 1, length = length(row),
                coverage = length / range)
    xt <- xt[xt$length > 6 & xt$coverage > 0.5,] 

}


#----------------------------------------------------------------------------
# ROUTE 1.2: CONVOLUTION WITH 'SPLIT' KERNEL
{
    k.size <- 5
    split.by <- 1
    
    k <- matrix(c(rep(-1/k.size, k.size * floor(k.size / 2)),
                  rep(0, k.size * split.by),
                  rep(k.size - 1, k.size),
                  rep(0, k.size * split.by),
                  rep(-1/k.size, k.size * floor(k.size / 2))),
                nrow = k.size)
           
    conv.sp <- array(apply(md7, 3:4, function(im) r2m(focal(m2r(im), k))),
                     dim = dim(md7), dimnames = dimnames(md7))
    
    pixel.image(conv.sp[,,"black", "MCT225"], xlim = c(400,500), ylim = c(1100, 1200))
    o.plot(conv.sp[,1160, "grey", "MCT225"], xlim = c(400,500), ylim = c())
    o.plot(md7[,1160, "grey", "MCT225"], xlim = c(400,500), col = "blue", add = T)
    abline(h = -10000, col = "red")
    
    pixel.image(md7[,,"grey", "MCT225"], xlim = c(400,500), ylim = c(1100, 1200))
    points(which(abs(conv.sp[,,"grey", "MCT225"]) > 10000, arr.ind = T), pch = 0)
    plot(conv.sp[,1160, "black", "MCT225"], type = "l", xlim = c(400,500))
    # finds pair of defective columns, with 2-px edge
    
    overplot(conv.sp, 449, "MCT225", xlim = c(1024, 2048))
    overplot(conv.sp, 750, "MCT225", xlim = c(1024, 2048))

    plot(conv.sp[449,,"black", "MCT225"], xlim = c(1024, 2048), type = "l")
    

}


#----------------------------------------------------------------------------
# ROUTE 1.3: CONVOLUTION WITH VARIABLE-LINE-WIDTH KERNEL
{
    k.size <- 5
    col.width <- 2
    
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)),
                  rep((k.size - 1) / col.width, k.size * col.width),
                  rep(-1, k.size * floor(k.size / 2))),
                nrow = k.size)
    if (ncol(k) %% 2 == 0) k <- cbind(-1, k)        # pad if necessary
    
    conv.sp <- array(apply(md7, 3:4, function(im) r2m(focal(m2r(im), k))),
                     dim = dim(md7), dimnames = dimnames(md7))
    
    pixel.image(conv.sp[,,"black", "MCT225"], xlim = c(400,500), ylim = c(1100, 1200))
    o.plot(conv.sp[,1160, "grey", "MCT225"], xlim = c(400,500), ylim = c())
    o.plot(md7[,1160, "grey", "MCT225"], xlim = c(400,500), col = "blue", add = T)
    abline(h = c(-1,1) * 100000, col = "red")
    
    pixel.image(md7[,,"grey", "MCT225"], xlim = c(400,500), ylim = c(1100, 1200))
    points(which(abs(conv.sp[,,"grey", "MCT225"]) > 100000, arr.ind = T), pch = 0)
    plot(conv.sp[,1160, "black", "MCT225"], type = "l", xlim = c(400,500))
    # finds pair of defective columns, with 2-px edge
    
    overplot(conv.sp, 449, "MCT225", xlim = c(1024, 2048))
    overplot(conv.sp, 750, "MCT225", xlim = c(1024, 2048))
    
    plot(conv.sp[449,,"black", "MCT225"], xlim = c(1024, 2048), type = "l")
    
    pixel.plot(which(abs(conv.sp[,,"grey","MCT225"]) > 100000, arr.ind = T), cex = 0.4)
    {
        plot(pw.m[256,,"grey", "MCT225"], type = "l")
        lines(pw.m[255,,"grey", "MCT225"], col = "blue")
        lines(pw.m[257,,"grey", "MCT225"], col = "darkred")
        
        points(pw.m[256,,"grey", "MCT225"], pch = ".", cex = 2,
               col = c(NA, "red")[(abs(conv.sp[256,,"grey","MCT225"]) > 100000) + 1])
        plot(pw.m[257,,"grey", "MCT225"], pch = ".", cex = 2,
               col = c(NA, "cyan3")[(abs(conv.sp[257,,"grey","MCT225"]) > 100000) + 1])
        
        o.plot(conv.sp[255,,"grey", "MCT225"])
        
        plot(conv.sp[256,,"grey", "MCT225"], type = "l", ylim = c(-150000, 100000))
        lines(conv.sp[257,,"grey", "MCT225"], col = "blue")
    }
}


#----------------------------------------------------------------------------
# ROUTE 2: DIRECT THRESHOLDING OF MEDIAN DIFFERENCES
# more intuitive threshold setting
{
    th <- array(apply(md7, 4,
                      function(im) abs(im[,,"black"]) > 300 | abs(im[,,"grey"]) > 300),
                dim = dim(md7[,,1,]), dimnames = dimnames(md7[,,1,]))
    
    pixel.image(th[,,"160430"])
    pixel.plot(which(abs(md7[,,"black", "160430"]) > 300, arr.ind = T))
    pixel.plot(which(abs(md7[,,"grey", "160430"]) > 300, arr.ind = T))
    pixel.plot(which(abs(md7[,,"white", "MCT225"]) > 300, arr.ind = T))
    
    xt <- ddply(data.frame(which(abs(conv) > th, arr.ind = T)), .(x = row), summarise,
                ymin = min(col), ymax = max(col), range = ymax - ymin + 1, length = length(row),
                coverage = length / range)
    xt <- xt[xt$length > 6 & xt$coverage > 0.5,] 
}


#----------------------------------------------------------------------------
# ROUTE 3: CHANGEPOINT ANALYSIS OF MEDIAN DIFFERENCES
# threshold using median differences first to reduce # candidates (also to detrend data)
# changepoint analysis to cut columns
# comparison of mean above and below changepoint 

# test bcp on all columns
# NB NEED TO REMOVE NA VALUES FIRST
# struggles when applied to whole column
{
    samp <- md7[c(28:32, 256, 257, 448, 449, 450), 28:2021,,"MCT225"]
    system.time(bb <- apply(samp, 1, bcp))             # ~ 1105 for whole image
    unlist(lapply(bb, function(bp) which.max(bp$posterior.prob)))
    
    plot(bb[[6]]$data[,3], type = "l")
    abline(v = 28 + which.max(bb[[6]]$posterior.prob), col = "red")
 
    plot(bb[[7]]$data[,3], type = "l")
    abline(v = 28 + which.max(bb[[7]]$posterior.prob), col = "red")   
    
    plot(bb[[8]]$data[,3], type = "l")
    abline(v = 28 + which.max(bb[[8]]$posterior.prob), col = "red")   
}

# test ecp on all columns
# very slow (stopped after ~1239 while attempting to do only 10 columns)
# NB NEED TO REMOVE NA VALUES FIRST
# struggles when applied to whole column
{
    samp <- md7[c(28:32, 256, 257, 448, 449, 450), 28:2021,,"MCT225"]
    ee[["256"]] <- e.divisive(md7[256, 28:2021, ,"MCT225"])
    
    tt <- apply(md7[28:2021,,"grey", "MCT225"], 1, function(cc) t.test(cc)$p.value)
    
    t.test(md7[256,,"grey", "MCT225"])
    tt <- t.test(md7[256, , ,"MCT225"])
    t.test(md7[750, , ,"MCT225"])$p.value
    
    overplot(md7, 256, "MCT225")
    abline(v = 28 + ee$"256"$order.found[1:3], col = "red")
    
    plot(bb[[7]]$data[,3], type = "l")
    abline(v = 28 + which.max(bb[[7]]$posterior.prob), col = "red")   
    
    plot(bb[[8]]$data[,3], type = "l")
    abline(v = 28 + which.max(bb[[8]]$posterior.prob), col = "red")   
}

# package comparisons
{
    vv <- md7[449,,,"MCT225"][!is.na(md7[449,,,"MCT225"])]  # 2 fairly clear changepoints
    dd <- md7[750,,,"MCT225"]  # dummy line - no discernable changepionts
    
    # possible package: ECP
    # advantage: multivariate, more selective than most
    # doesn't find changepoints where mean is constant; finds many where mean changes (ie. on damaged lines)
    # quite slow (~ 37s per 2048-column) - don't want to run on all columns if can avoid! (~ 20m per image?)
    {
        library(ecp)
        # gives no measure of confidence (not sure how goodness-of-fit measure should be interpreted)
        
        zz <- e.cp3o(vv, K = 5)
        overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + zz$estimates, col = "cyan3")
        
        zz.d <- e.cp3o(dd, K = 5)
        overplot(md7, 750, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + zz.d$estimates, col = "cyan3") 
        # found an erroneous breakpoint
        
        qq <- e.agglo(md7[449,1024:2021,,"MCT225"])
        qq.d <- e.agglo(dd)
        abline(v = 1024 + qq.d$estimates, col = "magenta3")
        # more flexible: found start, end & changepoint when there is one.
        
        # in flat data, found 800+ changepoints.
        
        system.time(aa <- e.divisive(md7[256,1024:2021,,"MCT225"]))
        overplot(md7, 256, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + aa$estimates, col = "green", lty = 2)
        abline(v = 1024 + aa$order.found[3], col = "red")
        # finds far too many changepoints (although gives order found, so can check sequentially)
        aa.d <- e.divisive(dd)
        overplot(md7, 750, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + aa.d$estimates, col = "cyan3")
        # found no changepoints in column with no changepoints
        vv[114,] - vv[113,]
        vv[144,] - vv[143,]
        vv[176,] - vv[175,]
        # may be useful if check magnitude of change at each point.
    }
    
    # possible package: BCP
    # use to get prob. of most likely 
    {
        library(bcp)
        bb <- bcp(vv)
        
        overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + which(bb$posterior.prob > 0.95), col = "gold")
        abline(v = 1024 + which(bb$posterior.prob > 0.999), col = "magenta3")
        
        abline(v = 24 + which.max(bb$posterior.prob), col = "red")
        
        
        bb.d <- bcp(dd)
        overplot(md7, 750, "MCT225", xlim = c(1024, 2048))
        abline(v = 1024 + which(bb.d$posterior.prob > 0.95), col = "gold")
        plot(bb.d$posterior.prob)
    }
    
    # rejected packages
    {
        # possible package: CPM
        # doesn't give posterior probability
        # also selects a changepoint slightly before the actual change
        {
            library(cpm)
            
            cc <- detectChangePoint(vv, cpmType = "GLR")
            overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
            abline(v = 1024 + cc$changePoint, col = "red")
            
            cc2 <- detectChangePoint(vv[cc$changePoint:998,2], cpmType = "Student")
            abline(v = cc$changePoint + cc2$changePoint, col = "orange")
            
            cc3 <- detectChangePoint(vv[cc2$changePoint:998,2], cpmType = "Student")
            abline(v = cc$changePoint + cc2$changePoint + cc3$changePoint, col = "gold")
            
            
            cm <- processStream(vv[,2], cpmType = "Cramer-von-Mises")
        }
        
        # possible package: changepoint
        # unable to handle multivariate data. Bleh.
        {
            library(changepoint)
            
            cp <- cpt.mean(vv)
            overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
            abline(v = 1024 + cc$changePoint, col = "red")
            
            cc2 <- detectChangePoint(vv[cc$changePoint:998,2], cpmType = "Student")
            abline(v = cc$changePoint + cc2$changePoint, col = "orange")
            
            cc3 <- detectChangePoint(vv[cc2$changePoint:998,2], cpmType = "Student")
            abline(v = cc$changePoint + cc2$changePoint + cc3$changePoint, col = "gold")
            
            
            cm <- processStream(vv[,2], cpmType = "Cramer-von-Mises")
        }
        
        # possible package: strucchange
        # picked an arbitrary breakpoint, missed the big one. Rejected.
        {
            library(strucchange)
            ss <- breakpoints(ts(vv[,1]) ~ 1)
            overplot(md7, 449, "MCT225", xlim = c(1024, 2048))
            abline(v = 1024 + ss$breakpoints, col = "red")
        }
    }
}


#----------------------------------------------------------------------------


# WHAT ABOUT DOING EVERYTHING OVER MEDIAN-DIFFERENCED IMAGE? WHY NOT?

pixel.plot(which(md7[,,"black", "160430"] > 1500, arr.ind = T), cex = 0.4)
pixel.plot(which(md7[,,"white", "160430"] > 1500, arr.ind = T), cex = 0.4)

pixel.plot(which(md7[,,"white", "140128"] > 5000, arr.ind = T), cex = 0.4)
pixel.plot(which(md7[,,"white", "140128"] > 5000, arr.ind = T), cex = 0.4)

pixel.image(md7[,,"white", "140128"])
pixel.image(md7[,,"white", "MCT225"])

hh <- count(round(c(md7[,,"white", "160430"]), 0))

smoothScatter(hh$x, cumsum(hh$freq), nrpoints = Inf, xlab = "Median difference", 
              ylab = "# pixels", xlim = c(-2000, 2000))
abline(h = 1996^2 / 2, col = "darkgrey")
abline(v = 0, col = "darkgrey")
abline(v = c(-1,1) * 1000, col = "darkred", lty = 2)

####################################################################################################

# TEST OVER ALL FILES                                                                           ####

# support functions found in next section
new.lines <- function(im, md, mid.break = 1024.5, k.size = 5, th = 3000) {
    
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    conv <- r2m(focal(m2r(md), k))
    
    # threshold using fairly loose limit (3000 rather than 5500)
    xt <- ddply(data.frame(which(abs(conv) > th, arr.ind = T)), .(x = row), summarise,
                ymin = min(col), ymax = max(col), range = ymax - ymin + 1, length = length(row),
                coverage = length / range)
    xt <- xt[xt$length > 6 & xt$coverage > 0.5,]    
    
    # get candidate lines
    find.changepoints(xt = xt, im = im, midline = mid.break)
}

# thresholding may need to be tied to SD of original image, for easier thresholding
# how can we compare likelihood of changepoint to likelihood of fitted parametric?
# residuals of changepoint vs residuals of smoothed line?

# 131122
{
    cp.b.131122 <- new.lines(pw.m[,,"black", "131122"], md7[,,"black", "131122"])
    cp.g.131122 <- new.lines(pw.m[,,"grey", "131122"], md7[,,"grey", "131122"])
    
    plot.changepoints(cp.g.131122, md7[,,"grey", "131122"])
    
    # still picking up a lot of healthy pixels when run over grey images.
    # needs more work.
}

# 140128
{
    # 745 known to be damaged
    cp.b.140128 <- new.lines(pw.m[,,"black", "140128"], md7[,,"black", "140128"])
    cp.g.140128 <- new.lines(pw.m[,,"grey", "140128"], md7[,,"grey", "140128"])
    
    plot.changepoints(cp.g.140128, pw.m[,,"grey", "140128"])
}
####################################################################################################

# CONVOLUTION OVER SMOOTHED RESIDUALS                                                           ####
im <- pw.m[,,"black", "131122"]

# smoothed images (kernel size 7)       ~ 74 seconds (~ 205 for three-panel image)
{
    md7 <- im - r2m(focal(m2r(im), matrix(rep(1, 49), ncol = 7), fun = median))
}

# convolution over smoothed vs convolution over pw mean image
{
    ll.org <- find.lines(im, dim.lines = T, threshold = 4000)
    ll.res <- find.lines(md7, dim.lines = T, threshold = 4000)
    
    ddply(data.frame(which(ll.org > 0, arr.ind = T)), .(x = row), summarise,
          ymin = min(col), ymax = max(col), length = length(row))
    
    ddply(data.frame(which(ll.res > 0, arr.ind = T)), .(x = row), summarise,
          ymin = min(col), ymax = max(col), length = length(row))
    
    # same lines identified in bright images
}

#-----------------------------------------------------------------

# convolve lines
{
    k.size <- 5
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    conv <- r2m(focal(m2r(md7), k))
    
    pixel.image(conv)
}

# threshold using fairly loose limit (3000 rather than 5500)
{
    xt <- ddply(data.frame(which(abs(conv) > 3000, arr.ind = T)), .(x = row), summarise,
          ymin = min(col), ymax = max(col), range = ymax - ymin + 1, length = length(row),
          coverage = length / range)
    xt <- xt[xt$length > 6 & xt$coverage > 0.5,]
}

# use change-point model along all candidate columns to assess line ends
{
    library(changepoint)
    cp <- cpt.mean(im[232,1025:2048][!is.na(im[232,1025:2048])], class = F)
    
    plot(im[232,], xlim = c(1025, 2048), type = "l")
    lines(im[231,], col = adjustcolor("forestgreen", alpha = 0.4))
    lines(im[233,], col = adjustcolor("dodgerblue4", alpha = 0.4))
    abline(v = 1024 + cp@cpts[1], col = "red")
    
    # check mean value above & below cutpoint
    # quick changepoint examples (from help file)
    {
        x = c(rnorm(100,0,1),rnorm(100,10,1))
        cpt.mean(x,penalty="SIC",method="AMOC",class=FALSE)
        # conf.value 1 means 99% confident
        cpt.mean(x,penalty="Asymptotic",pen.value=0.01,method="AMOC", class = F)
        cpts(ans)
        
        set.seed(1)
        x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
        cpt.mean(x,penalty="Manual",pen.value="2*log(n)",method="BinSeg",Q=5,class=FALSE)
    }
    
    find.changepoints <- function(xt, im, midline = 1024.5) {
        
        # iterate over lower panels, then upper (candidates may cross midline)
        # will need to adjust to account for panels with no midline
        if (is.na(midline)) {
            xt.lower <- xt
        } else {
            xt.lower <- xt[xt$ymin < midline,]
            xt.upper <- xt[xt$ymax > midline,]
        }

        for (i in 1:nrow(xt.lower)) {
            
            cc <- xt.lower[i,"x"]
            vv <- c(1:floor(midline))[!is.na(im[cc, 1:floor(midline)])]
            cp <- cpt.mean(im[cc, vv])
            
            # if changepoint is at end of column, discard
            if (min(cpts(cp), abs(cpts(cp) - length(vv))) > 2) {
                xt.lower$cpoint[i] <- cpts(cp) + min(vv) - 0.5
                xt.lower$mean1[i] <- param.est(cp)$mean[1]
                xt.lower$mean2[i] <- param.est(cp)$mean[2]
            } else {
                xt.lower$cpoint[i] <- NA
                xt.lower$mean1[i] <- NA
                xt.lower$mean2[i] <- NA
            }
        }
        
        if (!is.na(midline)) {
            for (i in 1:nrow(xt.upper)) {
                
                cc <- xt.upper[i,"x"]
                vv <- c(ceiling(midline):nrow(im))[!is.na(im[cc, ceiling(midline):nrow(im)])]
                cp <- cpt.mean(im[cc, vv])
                # if changepoint is at end of column, discard
            if (min(cpts(cp), abs(cpts(cp) - length(vv))) > 2) {
                xt.upper$cpoint[i] <- cpts(cp) + min(vv) - 0.5
                xt.upper$mean1[i] <- param.est(cp)$mean[1]
                xt.upper$mean2[i] <- param.est(cp)$mean[2]
            } else {
                xt.upper$cpoint[i] <- NA
                xt.upper$mean1[i] <- NA
                xt.upper$mean2[i] <- NA
            }
            }
            xt.lower <- rbind(xt.lower, xt.upper)
        }
        xt.lower <- xt.lower[!is.na(xt.lower$cpoint),]
        xt.lower <- xt.lower[abs(xt.lower$mean1 - xt.lower$mean2) > 100,]
        
        data.frame(col = xt.lower$x,
                   change = xt.lower$cpoint,
                   mean1 = xt.lower$mean1,
                   mean2 = xt.lower$mean2)
    }
    
    cand <- find.changepoints(xt, im)
    
    plot.changepoints <- function(cp, im, midline = 1024.5) {
        
        if (is.na(midline)) midline <- nrow(im) 
        
        for (i in 1:nrow(cp)) {
            plot(im[cp$col[i], ], type = "l", xlab = "", ylab = "", 
                 xlim = (cp$change[i] > midline) * floor(midline) + c(1, floor(midline)),
                 main = paste0("Column ", cp$col[i]))
            abline(v = cp$change[i], col = "red")
            
            m1 <- median(im[cp$col[i], 
                            ((cp$change[i] > midline) * floor(midline) + 1):cp$change[i]], na.rm = T)
            m2 <- median(im[cp$col[i], 
                            ceiling(cp$change[i]) : ((cp$change[i] > midline) * floor(midline) + 1024)], na.rm = T)
            
            lines(c(1, cp$change[i], NA, cp$change[i], nrow(im)),
                  c(m1, m1, NA, m2, m2),
                  col = "cyan3")
        }
    }
    plot.changepoints(cand, im)
}

# check mean/median difference from neighbouring pixels in each line segment to classify
# also check segments with no changepoint - may be whole column offset!
{
    classify.segments <- function(ll, md, midline = 1024.5) {
        
        if (is.na(midline)) midline <- nrow(md) 
        ss <- data.frame()
        
        for (i in 1:nrow(ll)) {
            
            # check first segment
            r1 <- ((ll$change[i] > midline) * floor(midline) + 1) : ll$change[i]
            d1 <- median(md[ll$col[i], r1] - md[ll$col[i] - 1, r1], na.rm = T)
            ss <- rbind(ss, c(ll$col[i],  range(r1), d1))
            
            # check second segment
            r2 <- ceiling(ll$change[i]) : ((ll$change[i] > midline) * floor(midline) + 1024)
            d2 <- median(md[ll$col[i], r2] - md[ll$col[i] - 1, r2], na.rm = T)
            ss <- rbind(ss, c(ll$col[i], range(r2), d2))
        }
        colnames(ss) <- c("col", "from", "to", "med.diff")
        
        ss$type <- c("dim", "healthy", "bright")[findInterval(ss$med.diff,
                                                              c(-65535, -100, 100, 65535))]
        ss
    }

    classify.segments(cand, im)
}

# could we skip convolution and go straight to filtering median diffs along columns?
# median-differencing is very like filter convolution....

# BUT convolution allows short line segments to be filtered out

md7.px <- which(md7 > 150, arr.ind = T)
pixel.plot(md7.px)

pixel.plot(which(conv > 3000, arr.ind = T))

####################################################################################################

# EXAMPLES OF PACKAGE 'CHANGEPOINTS'                                                           ####


data("Lai2005fig4", package = "changepoint")
Lai.default <- cpt.mean(Lai2005fig4[, 5], method = "PELT")
plot(Lai.default, pch = 20, col = "grey", cpt.col = "black", type = "p",
     xlab = "Index")
cpts(Lai.default)

coef(Lai.default)

#----------------------------------------

data("wind", package = "gstat")
ts.plot(wind[, 11], xlab = "Index")

wind.pelt <- cpt.var(diff(wind[, 11]), method = "PELT")
plot(wind.pelt, xlab = "Index")
logLik(wind.pelt)

wind.bs <- cpt.var(diff(wind[, 11]), method = "BinSeg")
ncpts(wind.bs)

# ALTERNATIVE PRE-PROCESSING: LOOK FOR LINES VS SMOOTHED NEIGHBOURS                             ####

ms7.160430 <- pw.m[,,,"160430"] - md7[,,,"160430"]
o.plot(pw.m[429,,"black", "160430"] - pw.m[428,,"black", "160430"],
       xlim = c(1024, 2048), ylim = c(-500, 1500))

plot(pw.m[429,,"black", "160430"], xlim = c(1024, 2048), type = "l", ylim = c(4500,6500))
lines(ms7.160430[429,,"black"], col = "blue")

plot(md7[429,,"black", "160430"], md7[429,,"black", "160430"], type = "l", ylim = c(-500,1500))

ma.pre <- movingFun(md7[429,,"black", "160430"], n = 5, fun = mean, type = "to")
ma.post <- movingFun(md7[429,,"black", "160430"], n = 5, fun = mean, type = "from")

plot(md7[429,,"black", "160430"], ylim = c(-500,15000), xlim = c(1024, 2048), type = "l")
lines(ma.pre * ma.post, col = "blue")

# again: smears effect of individual v. bright pixels

####################################################################################################

####################################################################################################

# TEST CASE FOR SPLIT-KERNEL CONVOLUTION                                                        ####

tst <- array(0, dim = c(23, 23))
tst[12,] <- 300
pixel.image(tst)

rr <- r2m(focal(m2r(tst), k))
pixel.image(rr)
    