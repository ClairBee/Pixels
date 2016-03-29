
library("IO.Pixels")

load.pixel.maps()

bpm.params <- xmlToList(xmlParse("./Other-data/Other-images/BadPixelMap_160314/CalibrationParameters.xml"))

# extract bad pixel map using alternative methods, compare to 'official' map
# 'official' map is from 16-03-14, so apply tests to same data set initially
img.date <- "160314"
    
m.b <- pw.b[,,img.date]
m.g <- pw.g[,,img.date]
m.w <- pw.w[,,img.date]

sd.b <- pw.sd.b[,,img.date]
sd.g <- pw.sd.g[,,img.date]
sd.w <- pw.sd.w[,,img.date]

# start by getting absolute max & min points
reset.bp <- function() {
    bp <<- rbind(data.frame(which(m.b == 0, arr.ind = T), src = "black", type = "dead"),
                data.frame(which(m.b == 65535, arr.ind = T), src = "black", type = "hot"),
                data.frame(which(m.g == 0, arr.ind = T), src = "grey", type = "dead"),
                data.frame(which(m.g == 65535, arr.ind = T), src = "grey", type = "hot"),
                data.frame(which(m.w == 0, arr.ind = T), src = "white", type = "dead"),
                data.frame(which(m.w == 65535, arr.ind = T), src = "white", type = "hot"))
}

    
# tmp function to quickly overlay 'official' & found points        
compare.maps <- function() {
    if (nrow(bp) > 5000) {
        print("Too many rows. Check data manually first")

    } else {
        bp.coords <- bp[!duplicated(bp[,c(1,2)]),c(1,2,4)]
        
        cols <- c("blue", "red", "gold", "green3", "purple")
        
        plot(bpm[,1:2], xlim = c(1,1996), ylim = c(1, 2100), asp = T)
        points(bp.coords[,1:2], pch = 20, col = cols[bp.coords$type])
        
        legend("top", horiz = T, legend = levels(bp.coords$type), pch = 20, col = cols, bty = "n")
    }
    print(table(bp$src, bp$type))
}


###################################################################################################
#                      IDENTIFY BAD PIXELS BY COMPARISON TO BAD PIXEL MAP                         #
###################################################################################################

{
# get 'bright offset corrected image' using flat field correction
    corr <- 60000 * (m.g - m.b) / (m.w - m.b)
    corr[is.na(corr)] <- 0      # otherwise get NA where FF == D
    corr[corr < 0] <- 0         # negative values not possible (result from different image offset)
    
    pixel.image(corr)
    # still systematic variance per panel, but magnitude much reduced
    s.hist(corr)
    
    bp <- rbind(bp, 
                data.frame(which(corr > (median(corr) * 1.5), arr.ind = T), src = "offset", type = "bright"),
                data.frame(which(corr < (median(corr) * 0.45), arr.ind = T), src = "offset", type = "dim"))
    
    compare.maps()
    
    # addition of 'noisy' category picks up large number of extra (unmatched) points. Discard.
    #bp <- rbind(bp,
    #            data.frame(which(sd.b > (median(sd.b) * 6), arr.ind = T), src = "black", type = "noisy"),
    #            data.frame(which(sd.g > (median(sd.g) * 6), arr.ind = T), src = "grey", type = "noisy"),
    #            data.frame(which(sd.w > (median(sd.w) * 6), arr.ind = T), src = "white", type = "noisy"))


    # use abs. thresholds from calibration file?
    # grey thresholds pick up too many points, but black & white seem ok
    # fills in many of the blanks, plus some additional points
    bp <- rbind(bp,
                data.frame(which(m.b > bpm.params$BlackMaxThreshold, arr.ind = T), src = "black", type = "high"),
                data.frame(which(m.w > bpm.params$WhiteMaxThreshold, arr.ind = T), src = "white", type = "high"),
                data.frame(which(m.b < bpm.params$BlackMinThreshold, arr.ind = T), src = "black", type = "low"),
                data.frame(which(m.w < bpm.params$WhiteMinThreshold, arr.ind = T), src = "white", type = "low"))

    table(bp$src, bp$type)
    
    # trying to establish pattern behind thresholds...
    {
        ecdf.b <- ecdf(m.b)
        ecdf.g <- ecdf(m.g)
        ecdf.w <- ecdf(m.w)
        
        # neither symmetrical nor consistent in quantiles
        ecdf.b(c(bpm.params$BlackMinThreshold, bpm.params$BlackMaxThreshold))    
        ecdf.w(c(bpm.params$WhiteMinThreshold, bpm.params$WhiteMaxThreshold))    
        
        as.numeric(bpm.params$BlackMaxThreshold) / median(m.b); as.numeric(bpm.params$WhiteMaxThreshold) / median(m.w)
        as.numeric(bpm.params$BlackMinThreshold) / median(m.b); as.numeric(bpm.params$WhiteMinThreshold) / median(m.w)
        
        # not symmatric about median/mean. Also not consistent in terms of distance of threshold from median/mean.
        as.numeric(bpm.params$BlackMaxThreshold) - median(m.b);  median(m.b) - as.numeric(bpm.params$BlackMinThreshold)
        as.numeric(bpm.params$WhiteMaxThreshold) - median(m.w);  median(m.w) - as.numeric(bpm.params$WhiteMinThreshold)
        
        as.numeric(bpm.params$BlackMaxThreshold) - mean(m.b);  mean(m.b) - as.numeric(bpm.params$BlackMinThreshold)
        as.numeric(bpm.params$WhiteMaxThreshold) - mean(m.w);  mean(m.w) - as.numeric(bpm.params$WhiteMinThreshold)
        
        # not linearly related to SD of mean values
        (as.numeric(bpm.params$BlackMaxThreshold) - as.numeric(bpm.params$BlackMinThreshold)) / sd(m.b)
        (as.numeric(bpm.params$WhiteMaxThreshold) - as.numeric(bpm.params$WhiteMinThreshold)) / sd(m.w)
    }
    
    # and giving up on that. Easier to try to identify new thresholds, which will pick up same/similar values
}


###################################################################################################
#                               IDENTIFY BAD PIXELS BY QUANTILES                                  #
###################################################################################################

# use ECDF of each image, cut quantiles to get most extreme values


###################################################################################################
#                   IDENTIFY BAD PIXELS BY DIFFERENCING & FINDING OUTLIERS                        #
###################################################################################################

# compare differences across row & column, use to identify outliers



###################################################################################################
#                      IDENTIFY BAD PIXELS IN OFFSET CORRECTION WITH PANEL ADJ                    #
###################################################################################################



###################################################################################################
#                    FIT FULL PARAMETRIC MODELS AND USE TO IDENTIFY BAD PIXELS                    #
###################################################################################################