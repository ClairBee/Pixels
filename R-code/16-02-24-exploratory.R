
# investigate intra-panel dependency: pixelwise mean & sd for each panel
load.images(150828, "black")

pw.m.b.150828 <- pixelwise.mean(b.150828)
pw.sd.b.150828 <- pixelwise.sd(b.150828)

# quick visual check
pixel.image(pw.m.b.150828, title = "Pixelwise mean value, black channel 15-08-28")
show.panels()

pixel.image(pw.sd.b.150828, title = "Pixelwise standard deviation, black channel 15-08-28")
show.panels()

pixel.image(pw.m.b.150828[c(1:126), c(993:1996)], break.levels = sd.levels(pw.m.b.150828))

############################################################################################

# spatial plot of mean/variance ratio - may need to adjust levels/colours
pw.ratio.b.150828 <- pw.sd.b.150828^2 / pw.m.b.150828

pw.ratio.b.150828[which(pw.m.b.150828 == 0)] <- 0
sd.levels(pw.ratio.b.150828)

image(c(1:1996), c(1:1996), pw.ratio.b.150828, asp = T,
      breaks = quantile(pw.ratio.b.150828, c(0:10)/10, na.rm = T),
      col = sd.colours())

abline(h = panel.edges()$y[1:2] - 0.5); abline(v = panel.edges()$x[1:16] - 0.5)
# no readily discernible pattern

# per-panel boxplot of ratio
par(mfrow = c(2,16))
for (i in c(2,1)) {
    for (j in 1:(length(panel.coords$x)-1)) {
        boxplot(pw.ratio.b.150828[c(panel.coords$x[j] : (panel.coords$x[j+1]-1)),
                         c(panel.coords$y[i] : (panel.coords$y[i+1]-1))], add = T,
                xlab = paste(panel.coords$x[j], panel.coords$y[i], sep = ","))                
    }
}
par(mfrow = c(1,1))

############################################################################################

# convert pixelwise mean & SD into vector for each panel
panel.mean <- list()
panel.sd <- list()
panel.coords <- panel.edges()
k <- 1

for (i in c(2,1)) {
    for (j in 1:(length(panel.coords$x)-1)) {
        panel.mean[[k]] <- pw.m.b.150828[c(panel.coords$x[j] : (panel.coords$x[j+1]-1)),
                                         c(panel.coords$y[i] : (panel.coords$y[i+1]-1))]
        panel.sd[[k]] <-  pw.sd.b.150828[c(panel.coords$x[j] : (panel.coords$x[j+1]-1)),
                                        c(panel.coords$y[i] : (panel.coords$y[i+1]-1))]
        k <- k+1
    }
}



# panelwise plots of mean vs variance
pdf("./Plots/Panelwise mean vs sd - black 150828.pdf")
par(mfrow = c(2,16), mar = c(0,0,0,0))
x.lim = c(0,6000)
for (i in 1:length(panels)) {
    plot(panel.sd[[i]]^2, panel.mean[[i]], pch = 20, xlim = x.lim, xaxt = "none", yaxt = "none")
}
dev.off()

pdf("./Plots/Pixelwise mean - black 150828.pdf")
pixel.image(pw.m.b.150828)
abline(v = panel.coords$y[1:2] - 0.5); abline(v = panel.coords$x[1:16] - 0.5)
dev.off()

pdf("./Plots/Pixelwise sd - black 150828.pdf")
pixel.image(pw.sd.b.150828)
abline(h = panel.coords$y[1:2] - 0.5); abline(v = panel.coords$x[1:16] - 0.5)
dev.off()

par(mfrow = c(2,16), mar = c(0,0,0,0))
for (i in c(2,1)) {
    for (j in 1:(length(panel.coords$x)-1)) {
    boxplot(b.150828[c(panel.coords$x[j] : (panel.coords$x[j+1]-1)),
                     c(panel.coords$y[i] : (panel.coords$y[i+1]-1)),], add = T)                
    }
}
par(mfrow = c(1,1))

# check coordinates of panels extracted
c <- c()
for (i in c(2,1)) {
    for (j in 1:(length(panel.coords$x)-1)) {
        c[k] <- paste(panel.coords$x[j], panel.coords$y[i], sep = ",")
        k <- k+1 
    }
}

par(mfrow = c(2,16))
for (i in c(2,1)) {
    for (j in 1:(length(panel.coords$x)-1)) {
        pixel.image(pw.m.b.150828[c(panel.coords$x[j] : (panel.coords$x[j+1]-1)),
                                  c(panel.coords$y[i] : (panel.coords$y[i+1]-1))],
                    break.levels = sd.levels(pw.m.b.150828))                
    }
}
