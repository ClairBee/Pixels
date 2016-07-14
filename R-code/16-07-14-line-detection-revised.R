
# NEEDS MORE WORK ON CLASSIFICATION
# MAYBE LOOK AT MEDIAN NEIGHBOUR DIFFERENCE ALONG COLUMN SEGMENT...

# also try packages CPM and ECP.
# or fit a changepoint model manually (see 'refs' bookmarks for details)

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

# EXAMPLES OF PACKAGES 'CHANGEPOINTS'                                                           ####


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
