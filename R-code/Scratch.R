
library("IO.Pixels")
plot.params <- "type = 'o', pch = 20, cex = 0.7, main = ''"

# next steps after meeting on March 8th:
#   - model circular spot (assume bivariate Gaussian)
#   - model panelwise variation
#   - model 'standing wave' effect in consecutive readings along columns
#
#   - create new set of plotting levels, based on IQR and more even distribution of values
#
#   - do signs of differences generally go same way? Even if start point isn't necessarily same?
#   - correlation of adjacent differences? (column-wise & row-wise)

# check if 'standing wave' is common across whole of black image

load.images(150828, "black")
pw.m.b.150828 <- pixelwise.mean(b.150828)
pw.sd.b.150828 <- pixelwise.sd(b.150828)

# get differences in pixelwise mean
diffs <- pw.m.b.150828[, c(1:1995)] - pw.m.b.150828[, c(2:1996)]    # column diffs
h.diffs <- pw.m.b.150828[c(1:1995),] - pw.m.b.150828[c(2:1996),]    # row diffs

hist(diffs, breaks = "fd", xlim = c(-500, 500))
hist(h.diffs, breaks = "fd", xlim = c(-500, 500))
# standing wave not seen across rows, only across columns

# would be much easier to examine this over a single panel: 
# then could be sure of which dimension I'm looking at!
pixel.image(pw.m.b.150828)
x <- c(127, 254); y <- c(1, 992)
lines(as.matrix(cbind(x = c(x[1], x[1], x[2], x[2], x[1]),
                      y = c(y[1], y[2], y[2], y[1], y[1]))))

panel <- pw.m.b.150828[x[1]:x[2], y[1]:y[2]]
pixel.image(panel)

panel.diffs <- diffs[x[1]:x[2], y[1]:y[2]]
pixel.image(panel.diffs)

plot(panel[1,892:992], ylim = c(5250, 6750), type = 'o', pch = 20, cex = 0.7, main = '')  
points(panel[2,892:992], col = adjustcolor("red", alpha = 0.5), type = 'o', pch = 20, cex = 0.7, main = '')

plot(panel[2,1:100], ylim = c(5250, 6750), type = 'o', pch = 20, cex = 0.7, main = '')  
points(panel[3,1:100], col = adjustcolor("red", alpha = 0.5), type = 'o', pch = 20, cex = 0.7, main = '')

plot(panel[3,1:100], ylim = c(5250, 6750), type = 'o', pch = 20, cex = 0.7, main = '')  
points(panel[4,1:100], col = adjustcolor("red", alpha = 0.5), type = 'o', pch = 20, cex = 0.7, main = '')



plot(panel.diffs[1,1:100], type = "o", pch = 20, cex = 0.7,
     main = "Columnwise plot")
plot(panel.diffs[,1], type = "o", pch = 20, cex = 0.7,
     main = "Row-wise plot")

plot(diffs[1,], type = "o", pch = 20, cex = 0.7, main = "Column-wise diffs", ylim = c(-3500, 3500))
plot(diffs[,1], type = "o", pch = 20, cex = 0.7, main = "Row-wise diffs", ylim = c(-3500, 3500))

# 'standing wave' phenomenon running across rows: row-wise mean of column differences
plot(apply(diffs, 2, mean)[1:993],  type = "o", pch = 20, cex = 0.7)
plot(apply(diffs, 2, mean)[994:1996],  type = "o", pch = 20, cex = 0.7)



diffs.1 <- b.150828[, c(1:1995),1] - b.150828[, c(2:1996),1]

marks <- c(rep(20, 5), rep(18, 5), rep(17, 5), rep(15, 5))
cols <- rep(c("orange", "blue", "red", "darkgreen", "purple"), 5) 

plot(0, ylim = c(4800, 9300), xlim = c(0,1995), type = "n")
c <- 2
for (i in 1:20) {
    points(b.150828[c,,i], type = "o", cex = 0.5, 
           pch = marks[i], col = adjustcolor(cols[i], alpha = 0.5))
}
points(pwm.b.150828[c,], type = "o", pch = 20, cex = 0.7)

plot(pw.m.b.150828[2,], pw.sd.b.150828[2,], pch = 20, cex = 0.7)
abline(coef(line(pw.m.b.150828[2,], pw.sd.b.150828[2,])), col = "red")

# how does slope of fitted line change across rows? Across columns? 
# Compare to panel-wise results.


################################################################################
# trying to remember/work out how to fit a bivariate normal (circular Gaussian) model...
{
# lets first simulate a bivariate normal sample
library(MASS)
bivn <- mvrnorm(100000, mu = c(0, 0), Sigma = matrix(c(1996, 0, 0, 1996), 2))

# now we do a kernel density estimate
bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)

# now plot your results
contour(bivn.kde, asp = T)
image(bivn.kde)
persp(bivn.kde, phi = 45, theta = 30)

# fancy contour with image
image(bivn.kde); contour(bivn.kde, add = T)

# fancy perspective
persp(bivn.kde, phi = 45, theta = 30, shade = .1, border = NA)
}
################################################################################