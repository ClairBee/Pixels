
library("IO.Pixels")

# How does dark current oscillation vary across different acquisition dates? 
#   - does oscillation always go same way?
#   - does amplitude change with date?
#   - does mean value change with date?

#########################################################################################

# import all black images, get pixelwise means (elapsed: 291.958 for 10 dates)
dates <- list.dirs("./Image-data/", full.names = F, recursive = F)
dates <- dates[dates != "150702"]
n <- length(dates)

pw.m.b <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))

pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
    load.images(dates[i], "black")
    pw.m.b[,,i] <- pixelwise.mean(eval(parse(text = paste0("b.",dates[i]))))
    setTxtProgressBar(pb, i)
}
close(pb)
remove(i, pb)

saveRDS(pw.m.b, file = "./Other-data/Pixelwise-means-black.rds")

#########################################################################################
#                               LOAD PIXELWISE MEAN MATRIX                              #
#########################################################################################

pw.m.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
n <- dim(pw.m.b)[3]

#########################################################################################

# define panel start
x <- c(1023, 1150); y <- c(600, 992)
c <- 4
cols <- c("black", "blue", "purple", "red", "orange", 
          "gold", "green", "darkgreen", "seagreen", "blue")

transects <- pw.m.b[x[1]+c, y[1]:y[2], ]

# plot all pixelwise means
o.plot(transects[,1], xaxt = "n",
       ylim = c(floor(min(transects)), ceiling(max(transects))),
       main = paste0("Pixelwise mean along column ", c+x[1]))
axis(1, at = pretty(0:(y[2]-y[1]), h = 100), labels = pretty(y[1]:y[2], h = 100))

for (i in 2:n) {
    o.plot(pw.m.b[x[1]+c, y[1]:y[2], i], add = T, col = adjustcolor(cols[i], alpha = 0.3))
}

#########################################################################################
#                   IMAGE OFFSET: COMPARE SUCCESSIVE ACQUISITIONS                       #
#########################################################################################
# panel offsets are not consistent, neither is change between image sets.
# change within panels isn't even consistent (see histograms)

pw.m.diffs <- pw.m.b[ , , 2:n] - pw.m.b[ , , 1:(n-1)]

# set parameters
p <- panel.edges()
q <- c(0.025,  0.159, 0.841, 0.975)
q.cols <- c("cyan", "green", "darkgreen", "red", "gold")
hist.cols <- c("blue", "red", "darkgreen", "gold", 
               "blueviolet", "chartreuse4", "cyan2", "cornflowerblue")

# image plots of differences in mean pixelwise values between consecutive acquisition dates
for (k in 1:(n-1)) {
    qn <- round(qnorm(q, mean(pw.m.diffs[,,k]), sd(pw.m.diffs[,,k])),0)
    
    pdf(paste0("./Plots/Dark-change-image-", dimnames(pw.m.b)[[3]][k+1], ".pdf"))
    image(c(1:1996), c(1:1996), pw.m.diffs[,,k], asp = T, 
          breaks = c(-65535, qn, 65535),
          ylim = c(0, 2100), ylab = "", xlab = "", 
          col = q.cols,
          main = paste0("PW mean change between ", dimnames(pw.m.b)[[3]][k],
                        " and ", dimnames(pw.m.b)[[3]][k+1]))
    draw.panels()
    legend("top", cex = 0.8, pch = 15, horiz = T,
           legend = apply(cbind(round(c(-65535, qn), 0), 
                                round(c(qn, 65535), 0)),
                          1, paste, collapse = ":"),
           col = q.cols)
    dev.off()
}


# histogram showing per-panel differences and whole-image differences
for (k in 1:(n-1)) {
    xl <- c(floor(qnorm(0.005, mean(pw.m.diffs[,,k]), sd(pw.m.diffs[,,k]))/10)*10,
            ceiling(qnorm(0.995, mean(pw.m.diffs[,,k]), sd(pw.m.diffs[,,k]))/10)*10)
    qn <- round(qnorm(q, mean(pw.m.diffs[,,k]), sd(pw.m.diffs[,,k])),0)
    
    pdf(paste0("./Plots/Dark-change-hist-", dimnames(pw.m.b)[[3]][k+1], ".pdf"))
    hist(c(pw.m.diffs[,,k]), breaks = "fd", xlim = xl,
         main = paste0("PW mean change between ", dimnames(pw.m.b)[[3]][k],
                       " and ", dimnames(pw.m.b)[[3]][k+1]),
         xlab = "Change in pixelwise mean value", col = "black")
    abline(v = qn, col = "grey")
    
    for (i in c(2,1)) {
        for (j in c(1:16)) {
            hist(pw.m.diffs[p$x[j] : (p$x[j+1]-1),
                            p$y[i] : (p$y[i+1]-1), 
                            k], 
                 breaks = "fd", add = T,
                 col = adjustcolor(hist.cols[(j %% length(hist.cols)) + 1], alpha = 0.3),
                 border = adjustcolor(hist.cols[(j %% length(hist.cols)) + 1], alpha = 0.5))
        }
    }
    # add colours to match image
    points(xl[1]:xl[2],
           rep(-600, length(xl[1]:xl[2])), pch = 15,
           col = c(rep(q.cols[1], qn[1] - xl[1]),
                   rep(q.cols[2], qn[2] - qn[1]),
                   rep(q.cols[3], qn[3] - qn[2]),
                   rep(q.cols[4], qn[4] - qn[3]),
                   rep(q.cols[5], xl[2] - qn[4]), q.cols[5]))
    
    dev.off()
}



#########################################################################################
#                   COMPARE DIFFS BETWEEN SUCCESSIVE ACQUISITIONS                       #
#########################################################################################

# get columnwise diffs for all acquisition dates
diffs <- diffs <- pw.m.b[, c(1:1995),] - pw.m.b[, c(2:1996),]

o.plot(diffs[x[1]+c, y[1]:y[2], 1], xaxt = "n", ylim = c(-200,200),
       main = paste0("Pixelwise diffs along column ", c+x[1]))
axis(1, at = pretty(0:(y[2]-y[1]), h = 100), labels = pretty(y[1]:y[2], h = 100))

for (i in 2:n) {
    o.plot(diffs[x[1]+c, y[1]:y[2], i], add = T, col = adjustcolor(cols[i], alpha = 0.3))
}

diff.sd <- apply(diffs, c(1,2), sd)

mean(diff.sd[x[1]+c, y[1]:y[2]])

hist(diff.sd, breaks = "fd", xlim = c(0, 100))
image(c(1:1996), c(1:1995), diff.sd,# xlim = c(1,100), ylim = c(1,100),
      col = c("gold", "violetred", "black"), 
      breaks = c(0, 20, 100, 65535))

draw.panels()

#########################################################################################
#                               LOESS SMOOTHING ALONG TRANSECT                          #
#########################################################################################

# carry out smoothing, residual plotting etc along a single transect to prove out
smoothed <- list()
for (i in 1:n) {
    smoothed[[i]] <- lowess(transects[,i], f = 1/15)$y
}

o.plot(smoothed[[6]], add = T, col = "gold")

#########################################################################################
#                  EXTEND PER-COLUMN SMOOTHING & PLOTTING TO FULL IMAGE                 #
#########################################################################################

# remove loess-smoothed centre line to leave residuals
smoothed <- list()
res <- list()
for (i in 1:n) {
    smoothed[[i]] <- do.call(rbind, lapply(apply(pw.m.b[,,i], 1, lowess, f = 1/15), "[[", 2))
    res[[i]] <- pw.m.b[,,i] - smoothed[[i]]
}   # ideally, adjust this to use apply in 3 dimensions

o.plot(smoothed[[1]][x[1]+c, y[1]:y[2]], xaxt = "n",
       main = paste0("Loess smoothed values along column ", c+x[1]))

# plot residuals (should all have same/similar magnitude?)


# get differences between loess along column. Is there a constant offset in each channel?


# identify per-channel offset. Does this extend to per-panel offset or does it vary across panel?

#########################################################################################
# GET SHAPE OF DARK CURRENT WITHOUT PANEL INTERFERENCE
#########################################################################################

# get per-channel smoothed values and adjust to give smoothed shape of offset
# need to adjust for panels separately; smoothing over centre line makes no sense
mid.edge <- pw.m.b[,991:994,1]

o.plot(mid.edge[,1], ylim = c(4500, 6000))
midline <- apply(mid.edge, 1, mean) 

#########################################################################################

# Use points with low oscillation SD as 'fixed pattern noise':
# do differences between these points give a constant per-channel offset?
# can per-channel offset be broken into image offset & channel offset?