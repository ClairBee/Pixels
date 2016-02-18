
library("IO.Pixels")
if (getwd() !=  "~Pixels") {setwd("~/Pixels")}

#####################################################################################################

# Add a function to import all channels for one day, as well as all days for one channel?

# would be useful to extract title from object name (or somehow store this information attached to the array?)
# function to split data out into individual panels (compare distribution across panels)
# come up with a colour scheme based on distance from mean (above or below) & assign using function

#####################################################################################################

# Import daily snapshot
b.150828 <- load.acq(150828, "black")
g.150828 <- load.acq(150828, "grey")
w.150828 <- load.acq(150828, "white")

# pixelwise mean of all 20 acquisitions for each channel
pw.b.150828 <- pixelwise.mean(b.150828)
pw.g.150828 <- pixelwise.mean(g.150828)
pw.w.150828 <- pixelwise.mean(w.150828)

# comparison of times required by each plotting method:
system.time(pixel.image(pw.g.150828))        # 12.09 elapsed
system.time(pixel.contour(pw.g.150828))      # 125.94 elapsed
# pixel.image also produced a much smaller image BUT lower resolution can lead to hatched effect in pdf.

# filled contour plot for each channel, with panel edges added
pdf("./Plots/Black-150828.pdf")
    pixel.image(pw.b.150828, "Black channel mean, 15-08-28")
    show.panels()
    dev.off()

pdf("./Plots/Grey-150828.pdf")
    pixel.image(pw.g.150828, "Grey channel mean, 15-08-28")
    dev.off()

pdf("./Plots/White-150828.pdf")
    pixel.image(pw.w.150828, "White channel mean, 15-08-28")
    dev.off()

#####################################################################################################

# rough identification of panel edges: check by zooming in on small area

pixel.image(pw.b.150828, x.range = c(120:400), y.range = c(980:1020)); show.panels()
pixel.image(pw.g.150828, x.range = c(120:400), y.range = c(980:1020)); show.panels()
pixel.image(pw.w.150828, x.range = c(120:400), y.range = c(980:1020)); show.panels()

# check against different day's shots
b.141009 <- load.acq(141009, "black")
g.141009 <- load.daily(141009, "grey")
w.141009 <- load.daily(141009, "white")

pw.b.141009 <- pixelwise.mean(b.141009)
pw.g.141009 <- pixelwise.mean(g.141009)
pw.w.141009 <- pixelwise.mean(w.141009)


# visually inspect all panel margins
pdf("./Plots/Panel check 14-10-09.pdf")
par(mfrow = c(3,1))
for (n in 0:6) {
    pixel.image(pw.b.141009, x.range = c(120:400)+(128*n*2), y.range = c(980:1020),
                title = paste("Black channel, 14.10.09, offset by ",n*2, sep = "")); show.panels()
    pixel.image(pw.g.141009, x.range = c(120:400)+(128*n*2), y.range = c(980:1020), 
                title = paste("Grey channel, 14.10.09, offset by ",n*2, sep = "")); show.panels()
    pixel.image(pw.w.141009, x.range = c(120:400)+(128*n*2), y.range = c(980:1020),
                title = paste("White channel, 14.10.09, offset by ",n*2, sep = "")); show.panels()
}
par(mfrow = c(1,1))
dev.off()

# lines seem to coincide with changes in intensity. Assume panels are OK unless Jay suggests otherwise?

#####################################################################################################

# rough identification of bad pixels: set upper limit as > 5sd from mean value to cut down # to assess
ul.b.150828 <- max(mean(b.150828) + 3 * sd(b.150828), (1 + mean(b.150828))/2)

# get locations of pixels that ever went above the threshold
ind.high <- unique(which(b.150828 > ul.b.150828, arr.ind = T)[,c(1:2)])

# immediately exclude those whose mean is exactly 1.0
ind.stuck <- which(pw.b.150828 == 1.0, arr.ind = T)     # 125 'always-on' pixels
ind.high <- ind.high[!((ind.high[,1] + (ind.high[,2]/10000)) %in% (ind.stuck[,1] + (ind.stuck[,2]/10000))),]



# plot paths of all pixels to have recorded high values
pdf("./Plots/High values (black 150828).pdf")
par(mfrow = c(3,2), mar = c(1,2,0,0))
line.cols <- c("blue", "green", "purple", "pink", "orange", "gold", "turquoise","green4","grey")
for (i in 0:5) {
    plot(b.150828[c(ind.high[i+1,1]),c(ind.high[i+1,2]),], type = "o", lwd = 2,
         pch = 20, ylim = c(0,1), ylab = "", xlab = "", xaxt = "none")
    abline(h = ul.b.150828, col = "red")
    abline(h = mean(b.150828), col = "red", lwd = 2)
    for (j in 1:9) {
        lines(b.150828[c(ind.high[i+1+j,1]),c(ind.high[i+1+j,2]),], type = "o", lwd = 2, 
              col = line.cols[j], pch = 20)
    }
}
dev.off()

#----------------------------------------------------------------------------------------------------

# now do the same for very low values...
ll.b.150828 <- mean(b.150828) - 3 * sd(b.150828)

ind.low <- unique(which(b.150828 < ll.b.150828, arr.ind = T)[,c(1:2)])
ind.stuck.off <- which(pw.b.150828 == 0.0, arr.ind = T)     # 5 'always-off' pixels
ind.low <- ind.low[!((ind.low[,1] + (ind.low[,2]/10000)) %in% (ind.stuck.off[,1] + (ind.stuck.off[,2]/10000))),]

plot(b.150828[c(ind.low[1]),c(ind.low[2]),], type = "o", lwd = 2,
    pch = 20, ylim = c(0,1), ylab = "", xlab = "", xaxt = "none")
abline(h = ll.b.150828, col = "red")

abline(h = mean(b.150828), col = "red", lwd = 2)

#####################################################################################################

# summary of xml profiles: offset parameters are consistent
summ <- summarise.image.profiles()
summ <- summ[order(summ$acq.time),]
summ <- summ[order(summ$acq.date),]
summ

#####################################################################################################

# compare import methods: as.is = T or F

filenm <- "./Image-data/150828/black/black_280815_1.tif"

m.org <- readTIFF(filenm)
m.asis <- readTIFF(filenm, as.is = T)

system.time(summ.org <- batch.summary(m.org))           # 0.72 elapsed
system.time(summ.asis <- batch.summary(m.asis))         # 0.41 elapsed

# summary statistics are all the same, scaling appears to be linear: ok to use the 'as-is' import method
round(summ.org * 65535,9) == round(summ.asis,9)

#####################################################################################################

# get summary statistics for all images

system.time(z <- summarise.all())
