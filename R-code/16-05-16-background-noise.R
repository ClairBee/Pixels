
library("IO.Pixels"); library("CB.Misc")

load.pixel.means()
md.w <- readRDS("./Other-data/Median-diffs-white.rds")

####################################################################################################

# same background noise pattern can be seen in...
#   ...per-acquisition differences in black images
#   ...median differences in white images

# IS PATTERN THE SAME IN ALL IMAGES?

####################################################################################################
# the selected images used were chosen because the pattern was easy to see with the naked eye.
# need to be more thorough/robust if we're going to write anything up!
####################################################################################################

im.b <- pw.m[,,"black", "160430"] - pw.m[,,"black", "160314"]
im.w <- md.w[["141009"]]

s.hist(im.w)
th.w <- threshold(im.w, level = mean(im.w, na.rm = T) + 0.55 * sd(im.w, na.rm = T))
image(1:1996, 1:1996, th.w, col = c(NA, "magenta3"), xlim = c(0,500), ylim = c(0,500))

s.hist(im.b)
th.b <- threshold(im.b, level = mean(im.b, na.rm = T) + 1.5 * sd(im.b, na.rm = T))
image(1:1996, 1:1996, th.b, col = c(NA, adjustcolor("cyan3", alpha = 0.5)), xlim = c(0,500), ylim = c(0,500))

cor(c(im.b[1:500,1:500]), c(im.w[1:500,1:500]), use = "complete.obs")

####################################################################################################

# TEST FOR RANDOMNESS                                                                           ####

G <- Gest(as.ppp(which(im.b > 0, arr.ind = T), c(0.5, 1996.5, 0.5, 1996.5)))
plot(G)

K <- Kest(as.ppp(which(im.b > 0, arr.ind = T), c(0.5, 1996.5, 0.5, 1996.5)))
plot(K)

# perform the test itself with a 100 simulations
E <- envelope(as.ppp(which(im.b > 0, arr.ind = T), c(0.5, 1996.5, 0.5, 1996.5)), Kest)
plot(E, main = NULL)