
# DETECTION OF SLIGHTLY BRIGHT (OR DIM) LINES

library("IO.Pixels"); library("CB.Misc")
if (getwd() != "/home/clair/Documents/Pixels") {setwd("/home/clair/Documents/Pixels")}
fpath <- "./Notes/Line-detection/fig/"

load.pixel.means()

####################################################################################################

# SETUP                                                                                         ####

conv.sq5 <- alply(pw.m[,,"black", ], 3, convolve.lines, k = 5, .dims = T)
th.sq5 <- lapply(conv.sq5, threshold, level = 5500)
sm.sq5 <- lapply(th.sq5, smooth.lines, sm.size = 11)

image(sm.sq5$"160314", col = c(NA, "red"))

####################################################################################################

# finish data cleaning

# clean data, then match lines identified after smoothing to original thresholded convolution
# to pick up full length of lines


