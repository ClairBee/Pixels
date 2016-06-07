
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Drafts/State-space/fig/"

load.pixel.means()
load.pixel.sds()

md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

Cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")

####################################################################################################

# PLOT THRESHOLD SIGMA VS %                                                                     ####

# sigma = SD of per-array mean value, not pixelwise SD
sig <- apply(pw.m, 3:4, sd)

th.prop <- function(im, n.sigma, midpoint = "median", table = F) {
        if (midpoint == "median") {
            m <- median(im)
        } else {
            m <- mean(im)
        }
        
        s <- sd(im)
        if (length(n.sigma) > 1) {
            p <- unlist(lapply(n.sigma, function(x) 100 * sum(im > m + x * s) / length(im)))
        } else {
            p <- 100 * sum(im > m + n.sigma * s) / length(im)
        }
        
        if (table) {
            return(data.frame(n.sigma, 
                              prop = p,
                              GV = m + n.sigma * s))
        } else {
            return(p)
        }
    }

    ccols <- c("black", "blue", "purple", "magenta3", "orange", "gold", "green", "green3", "cyan3", "skyblue", "slateblue1", "orchid")
    ns <- c(2:10)
    
    pdfheight <- 5; pdfwidth = 7
    
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-black.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"black", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"black",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-black.pdf"))
    }
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-grey.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"grey", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"grey",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-grey.pdf"))
    }
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-white.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"white", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"white",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-white.pdf"))
    }

# get % at each setting by sigma
zz <- sapply(dimnames(pw.m)[[4]], function(i) th.prop(pw.m[,,"black", i], n.sigma = ns))
rownames(zz) <- ns
    
matplot(t(zz), type = "o", col = rev(ccols), lty = 1, pch = 20)

# rescale all rows to have same variance & centre. Very similar shape of increase at all levels.
matplot(scale(t(zz), center = T, scale = T), type = "l", col = rev(ccols), lty = 1,
        xaxt = "n")
axis(1, at = c(1:12), labels = sapply(dimnames(pw.m)[[4]], fancy.date), srt = 90)    

####################################################################################################

# CHECK DEVELOPMENT OF PER-PANEL GRADIENT OVER TIME                                             ####

# fit linear model to each image in turn

panel.lm.black <- lapply(apply(pw.m[,,"black",], 3, panel.lm, robust = T), "[", "models")

x.grad <- do.call("rbind", lapply(lapply(panel.lm.black, "[[", "models"), "[",,"x"))
y.grad <- do.call("rbind", lapply(lapply(panel.lm.black, "[[", "models"), "[",, "y"))

colnames(x.grad) <- colnames(y.grad) <- apply(cbind(c(rep("U", 16), rep("L", 16)), rep(c(1:16), 2)), 1, paste, collapse = "")

par(mfrow = c(2, 1))
matplot(x.grad, type = "l")
matplot(y.grad, type = "l")
par(mfrow = c(1,1))

ccols <- c("black", "blue", "purple", "magenta3", "orange", "gold", "green", "green3", "cyan3", "skyblue", "slateblue1", "orchid")

plot(x.grad[,"U1"], y.grad[,"U1"], type = "o", pch = 21, bg = ccols, xlim = c(0,5), ylim = c(0,1.5),
     xlab = "x-gradient fitted", ylab = "y-gradient fitted")
points(x.grad[,"U2"], y.grad[,"U2"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U3"], y.grad[,"U3"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U4"], y.grad[,"U4"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U5"], y.grad[,"U5"], type = "o", pch = 21, bg = ccols)

points(x.grad[,"U14"], y.grad[,"U14"], type = "o", pch = 21, bg = ccols)

plot(x.grad[1,1:16], type = "l", ylim = c(-7,5))
lines(x.grad[1,17:32])
abline(h = 0, lty = 2)
for (i in 2:12) {
    lines(x.grad[i,1:16], col = rev(ccols)[i])
    lines(x.grad[i,17:32], col = rev(ccols)[i])
}

plot(y.grad[1,1:16], type = "l", ylim = c(-1.5,1))
lines(y.grad[1,17:32])
abline(h = 0, lty = 2)
for (i in 2:12) {
    lines(y.grad[i,1:16], col = rev(ccols)[i])
    lines(y.grad[i,17:32], col = rev(ccols)[i])
}


####################################################################################################

# THRESHOLD BEST CASE VS WORST CASE                                                             ####

####################################################################################################

# CLASSIFY ALL PIXELS                                                                           ####

####################################################################################################

# FIT MARKOV MODEL IN R                                                                         ####

library(msm)

# need a dataframe with time of observation, observed state, and subject id.
# observations must be grouped by subject, and ordered by time within each subject.

# convert pixel map of types to data frame of all pixels, with pixel id.
# best way is probably to unlist and melt a bpx image.

qq <- melt(pw.m)
