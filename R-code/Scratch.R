
library("IO.Pixels")
library(xts)

load.images(150828, "black")
pw.m.b.150828 <- pixelwise.mean(b.150828)
pw.sd.b.150828 <- pixelwise.sd(b.150828)
lvls <- sd.levels(pw.m.b.150828)

iqr.levels <- function(data) {
    lq <- quantile(data, 0.25)
    med <- median(data)
    uq <- quantile(data, 0.75)
    iqr <- IQR(data)
    
    c(0,
      max(lq - 1.5 * iqr,         # low outliers
          lq/2),
      lq,
      lq + (med - lq) / 3,
      med - (med - lq) / 3,
      med,
      med + (uq - med) / 3,
      uq - (uq - med) / 3,
      uq,
      min(uq + 1.5 * iqr,    # high outliers
          uq + ((65535-uq) / 2)),
      65535)
}
a.lvls <- abs.levels(panel)

panel.lvls <- c(0,
                lq - 1.5 * iqr, # low outliers
                lq,
                lq + (med - lq) / 3,
                med - (med - lq) / 3,
                med,
                med + (uq - med) / 3,
                uq - (uq - med) / 3,
                uq,
                uq + 1.5 * iqr,    # high outliers
                65535)

# looking at common periodic behaviour across panels
x <- c(1023, 1150); y <- c(1,992)
pixel.image(pw.m.b.150828)
lines(as.matrix(cbind(x = c(x[1], x[1], x[2], x[2], x[1]),
                      y = c(y[1], y[2], y[2], y[1], y[1]))))

panel <- pw.m.b.150828[x[1]:x[2], y[1]:y[2]]
pixel.image(panel, break.levels = a.lvls)

plot(panel[4,], ylim = c(4500, 6000), type = 'o', pch = 20, cex = 0.7, main = '')  
points(panel[5,], col = adjustcolor("red", alpha = 0.5), type = 'o', pch = 20, cex = 0.7, main = '')
points(lowess(panel[4,], f = 1/5), col = "red", type = "l", lwd = 2)
plot(panel[4,] - lowess(panel[4,], f = 1/5)$y, type = "o", pch = 20)
mean(panel[4,] - lowess(panel[4,], f = 1/5)$y)
sd(panel[4,] - lowess(panel[4,], f = 1/5)$y)
abline(h = mean(panel[4,] - lowess(panel[4,], f = 1/5)$y), col = "red", lwd = 2)
abline(h = (sd(panel[4,] - lowess(panel[4,], f = 1/5)$y)) * c(-1,1), col = "purple", lwd = 2)
abline(h = (sd(panel[4,] - lowess(panel[4,], f = 1/5)$y)) * 1.5 * c(-1,1), col = "blue", lwd = 2)

res <- panel[4,] - lowess(panel[4,], f = 1/5)$y
median(abs(res))

med.res <- movingFun(abs(res), 5, median, 'around')

# now to apply this approach across whole panel


plot(panel[4,400:450], type = "o", pch = 20)
points(mean(panel[4,400:450]) + rep(c(1,-1) * 50, 51), col = adjustcolor("red", alpha = 0.5), type = "o")

osc <- rep(c(1,-1) * 50, 51)
pixel.image(panel[, 200:600])

acf.4 <- acf(panel[4,])$acf[-1]
acf.5 <- acf(panel[5,])$acf[-1]

# measure correlation between autocorrelation functions
cor(acf.4, acf.5)
cor(pacf(panel[4,])$acf, pacf(panel[5,])$acf)

acf <- list()
pacf <- list()
ccf <- list()

for (i in 1:128) {
    acf[[i]] <- acf(panel[i,], plot = F)$acf[-1]
    pacf[[i]] <- pacf(panel[i,], plot = F)$acf
    ccf[[i]] <- ccf(panel[i,], panel[i+1,], plot = F)$acf
}

acf.mat <- do.call(rbind, acf)
pacf.mat <- do.call(rbind, pacf)
ccf.mat <- do.call(rbind, ccf)

df <- cbind(acf.2 = acf.mat[,2],
            pacf.2 = pacf.mat[,2], 
            ccf.1 = c(ccf.mat[,28], NA),
            ccf.2 = c(ccf.mat[,29], NA),
            c(acf.cor, NA), c(pacf.cor, NA))

plot(df[,1], type = "o", pch = 20, ylab = "ACF", main = "Per-column ACF, lag 2")
plot(df[,2], type = "o", pch = 20, ylab = "PACF", main = "Per-column PACF, lag 2")

plot(df[,3], type = "o", pch = 20, ylab = "CCF", main = "Between-column CCF, lag 1")
plot(df[,4], type = "o", pch = 20, ylab = "CCF", main = "Between-column CCF, lag 2")

plot(df[,5], type = "o", pch = 20, ylab = "Cor", main = "ACF correlation across columns")
plot(df[,6], type = "o", pch = 20, ylab = "Cor", main = "PACF correlation across columns")

plot(0 ,type = "n", ylim = c(-1,1), xlim = c(0,128), ylab = "")
points(df[,1], type = "o", pch = 20)
points(df[,2], type = "o", pch = 20, col = adjustcolor("blue", alpha = 0.5))

points(df[,5], type = "o", pch = 20, col = adjustcolor("red", alpha = 0.5))
points(df[,6], type = "o", pch = 20, col = adjustcolor("orange", alpha = 0.5))



acf.cor <- c(); pacf.cor <- c()
for (i in 1:127) {
    acf.cor[i] <- cor(acf.mat[i,], acf.mat[i+1,])
    pacf.cor[i] <- cor(pacf.mat[i,], pacf.mat[i+1,])
}

    
    df[i,] <- c(acf(panel[i,])$acf[3],
                cor(acf(panel[i,])$acf, acf(panel[i+1,])$acf),
                pacf(panel[i,])$acf[2],
                cor(pacf(panel[i,])$acf, pacf(panel[i+1,])$acf),
                ccf(panel[i,], panel[i+1,]))
}

pacf <- pacf(panel[4,])$acf
ccf <- ccf(panel[4,], panel[5,])[2]

z <- arima(panel[4,], c(2,1,2))
plot(z$residuals)

plot(arima(panel[4,], c(2,1,2))$residuals)

#####################################################################################

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



# get differences in pixelwise mean
diffs <- pw.m.b.150828[, c(1:1995)] - pw.m.b.150828[, c(2:1996)]    # column diffs
h.diffs <- pw.m.b.150828[c(1:1995),] - pw.m.b.150828[c(2:1996),]    # row diffs

hist(diffs, breaks = "fd", xlim = c(-500, 500))
hist(h.diffs, breaks = "fd", xlim = c(-500, 500))
# standing wave not seen across rows, only across columns

# would be much easier to examine this over a single panel: 
# then could be sure of which dimension I'm looking at!
pixel.image(pw.m.b.150828)
lines(as.matrix(cbind(x = c(x[1], x[1], x[2], x[2], x[1]),
                      y = c(y[1], y[2], y[2], y[1], y[1]))))

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