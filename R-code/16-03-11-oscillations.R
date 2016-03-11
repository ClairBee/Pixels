
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
#                               LOAD PXIELWISE MEAN MATRIX                              #
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
#                       EXTEND SMOOTHING & PLOTTING TO FULL IMAGE                       #
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