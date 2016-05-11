
# DETECTION OF SLIGHTLY BRIGHT (OR DIM) LINES

library("IO.Pixels"); library("CB.Misc")
if (getwd() != "/home/clair/Documents/Pixels") {setwd("/home/clair/Documents/Pixels")}
fpath <- "./Notes/Line-detection/fig/"

load.pixel.means()

####################################################################################################

# THURSDAY MAY 12TH: TO DO

# run function over all black images: bright & dim
# fun function over all grey images: bright & dim
# run function over all white images: bright & dim

# for all lines found, plot
#   - development over time at each power level
#   - overplot convolutions over time: is difference between line & neighbours constant?
#   - try to associate lines with superclusters

# apply function to as many old images as possible
# repeat the above for any & all lines found

# write all this up & send it to Julia.

#-------------------------------------------------------------------------------------------------

# also have a quick look at matching that diagonal pattern across images
#   - perhaps use differences across rows/columns? Same pattern of increases/decreases?
#   - may need to normalise (but prob. not in black images?)
#   - try to find a particularly egregious area & focus on details
#   - maybe pattern of gradients is more useful than actual change in value?

####################################################################################################

# SETUP                                                                                         ####

bright.l <- alply(pw.m[,,"black", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
                  sm.size = 11, min.length = 6, .dims = T)

image(1:1996, 1:1996, bright.l$"160314", col = c(NA, "green3"), breaks = c(-0.5, 0.5, 999))

summarise.lines(bright.l$"160314")

####################################################################################################

# RUN OVER OLD DATA                                                                             ####

# try running over old data

# MAY NEED TO REPOSITION CENTRE LINE

tmp.black <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapBlack.tif", as.is = T)
tmp.black <- t(tmp.black[nrow(tmp.black):1, , drop = FALSE])

tmp.grey <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapGrey.tif", as.is = T)
tmp.grey <- t(tmp.grey[nrow(tmp.grey):1, , drop = FALSE])

tmp.white <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapWhite.tif", as.is = T)
tmp.white <- t(tmp.white[nrow(tmp.white):1, , drop = FALSE])

zz <- find.lines(tmp.black)

image(1:2000, 1:2000, zz, col = c(NA, "red"), xlim = c(970,1000), ylim = c(0,300))
# got one!

which(zz > 0, arr.ind = T)

pixel.image(tmp.black, xlim = c(970,1000), ylim = c(0,300))

o.plot(tmp.black[970:1000, 100]); o.plot(tmp.black[970:1000, 103], add = T); o.plot(tmp.black[970:1000, 97], add = T)

#-------------------------------------------------------------------------------
# dim lines?

qq <- find.lines(tmp.black, dim.lines = T)

table(qq) # two lines found
image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999))

summarise.lines(qq)

image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999), xlim = c(190, 220), ylim = c(900, 1300))
pixel.image(tmp.black, xlim = c(190, 220), ylim = c(900, 1300))
o.plot(tmp.black[190:230, 1100]); o.plot(tmp.black[190:230, 1101], add = T); o.plot(tmp.black[190:230, 1099], add = T);
o.plot(tmp.black[208, ], xlim = c(990, 1100))

image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999), xlim = c(1300, 1400), ylim = c(1000, 1100))
pixel.image(tmp.black,  xlim = c(1300, 1400), ylim = c(1000, 1100))
o.plot(tmp.black[1372, ], xlim = c(990, 1100))



# plots of line 1
{
    o.plot(tmp[208, ], xlim = c(990, 1300), ylim = c(0,50000))
    o.plot(tmp.grey[208, ], xlim = c(990, 1300), add = T, col = "grey")
    o.plot(tmp.white[208, ], xlim = c(990, 1300), add = T, col = "gold")
    abline(v = 1000.5, lwd = 2)                 # panel midline
    abline(v = 1247.5, col = "red", lty = 2)    # end of dim line
}

# plots of line 2
{
    o.plot(tmp[1372, ], xlim = c(990, 1100), ylim = c(0,50000))
    o.plot(tmp.grey[1372, ], xlim = c(990, 1100), add = T, col = "grey")
    o.plot(tmp.white[1372, ], xlim = c(990, 1100), add = T, col = "gold")
    abline(v = 1000.5, lwd = 2)                 # panel midline
    abline(v = 1071.5, col = "red", lty = 2)    # end of dim line
}



####################################################################################################

# SCRATCH/MISC                                                                                  ####

# plots of problematic areas
{
    # line segment retained in column 3
    {
        o.plot(pw.m[3, 1940:1970, "black", "160314"], ylim = c(0, 25000))
        o.plot(pw.m[4, 1940:1970, "black", "160314"], add = T, col = adjustcolor("green3", alpha = 0.3))
        o.plot(pw.m[2, 1940:1970, "black", "160314"], add = T, col = adjustcolor("gold", alpha = 0.3))
        
        o.plot(conv[3, 1940:1970], add = T, col = "cyan3")
        abline(h = 5500, col = "red", lty = 2)
        o.plot(sm[3, 1940:1970] * 2500, add = T, col = "blue")
        abline(h = 5.5 * 2500, col = "blue", lty = 2)
        
        pixel.image(pw.m[1:126, , "black", "160314"], xlim = c(0,30), ylim = c(1940,1970))
        rect(2.5, 1950.5, 3.5, 1965.5)
    }
    
    # line segment retained in column 234
    {
        o.plot(pw.m[234, 260:300, "black", "160314"], ylim = c(0, 25000))
        o.plot(pw.m[233, 260:300, "black", "160314"], add = T, col = adjustcolor("green3", alpha = 0.3))
        o.plot(pw.m[232, 260:300, "black", "160314"], add = T, col = adjustcolor("gold", alpha = 0.3))
        
        o.plot(conv[234, 260:300], add = T, col = "cyan3")
        abline(h = 5500, col = "red", lty = 2)
        o.plot(sm[234, 260:300] * 2500, add = T, col = "blue")
        abline(h = 5.5 * 2500, col = "blue", lty = 2)
        
        pixel.image(pw.m[200:300, , "black", "160314"], xlim = c(210,250), ylim = c(260,300))
        rect(233.5, 274.5, 234.5, 293.5)
    }
    
    # line identified in column 427
    {
        o.plot(pw.m[427,, "black", "160314"], ylim = c(-1000, 10000), xlim = c(1100,1996))
        o.plot(conv[427,], add = T, col = "green3")
        abline(h = 5500, col = "red")
        o.plot(sm[427,] * 900, add = T, col = "blue")
        abline(h = 5.5*900, col = "blue", lty = 2)
    }
    
    # line identified in column 809
    {
        o.plot(pw.m[809,, "black", "160314"], ylim = c(-1000, 10000), xlim = c(0,992))
        points(pw.m[810,, "black", "160314"], type = "l", col = "orange")
        points(pw.m[808,, "black", "160314"], type = "l", col = "gold")
        o.plot(conv[809,], add = T, col = "green3")
        abline(h = 5500, col = "red")
        o.plot(sm[809,] * 900, add = T, col = "blue")
        abline(h = 5.5*900, col = "blue", lty = 2)
    }
}
