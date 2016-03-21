
library("IO.Pixels")
library(reshape)

###########################################################################################
# FIT LINEAR MODELS ACROSS ALL SUBPANELS
pw.m.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
pw.m.g <- readRDS("./Other-data/Pixelwise-means-grey.rds")
pw.m.w <- readRDS("./Other-data/Pixelwise-means-white.rds")

pw.sd.b <- readRDS("./Other-data/Pixelwise-sds-black.rds")

{
    n <- dim(pw.m.b)[[3]]
    
    # get linear gradients across panels - dark images
    {
        wedges.b <- array(dim = c(32, 4, n), 
                          dimnames = list(apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                                      rep(c(1:16), 2)), 1, paste, collapse = ""),
                                          c("offset", "x", "y", "ratio"),
                                          dimnames(pw.m.b)[[3]]))
        smoothed.panels.tmp <- array(dim = c(128, 1024, 32))
        smoothed.panels.b <- array(dim = dim(pw.m.b))
        
        # takes less than a minute to run for 11 images
        for (l in 1:n) {
            tmp <- subpanels(pw.m.b[ , , l])
            for (s in 1:32) {
                lm <- lm(value ~ X1 + X2, data = melt(tmp[ , , s]))
                wedges.b[s, 1:3 , l] <- coef(lm)
                smoothed.panels.tmp[,,s] <- predict(lm, melt(tmp[ , , s]))
            }
            wedges.b[ , 4, l] <- wedges.b[ , 2, l] / wedges.b[ , 3, l]
            smoothed.panels.b[,,l] <- join.panels(smoothed.panels.tmp)
        }
        saveRDS(wedges.b, file = "./Other-data/Panel-gradients-black.rds")
        saveRDS(smoothed.panels.b, file = "./Other-data/Smoothed-panels-black.rds")
    }
    
    # get linear gradients across panels - grey images
    {
        wedges.g <- array(dim = c(32, 4, n), 
                          dimnames = list(apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                                      rep(c(1:16), 2)), 1, paste, collapse = ""),
                                          c("offset", "x", "y", "ratio"),
                                          dimnames(pw.m.g)[[3]]))
        
        # takes less than a minute to run for 11 images
        for (l in 1:n) {
            tmp <- subpanels(pw.m.g[ , , l])
            for (s in 1:32) {
                wedges.g[s, 1:3 , l] <- coef(lm(value ~ X1 + X2, data = melt(tmp[ , , s])))
            }
            wedges.g[ , 4, l] <- wedges.g[ , 2, l] / wedges.g[ , 3, l]
        }
        saveRDS(wedges.g, file = "./Other-data/Panel-gradients-grey.rds")
    }
    
    # get linear gradients across panels - white images
    {
        wedges.w <- array(dim = c(32, 4, n), 
                          dimnames = list(apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                                      rep(c(1:16), 2)), 1, paste, collapse = ""),
                                          c("offset", "x", "y", "ratio"),
                                          dimnames(pw.m.w)[[3]]))
        
        # takes less than a minute to run for 11 images
        for (l in 1:n) {
            tmp <- subpanels(pw.m.w[ , , l])
            for (s in 1:32) {
                wedges.w[s, 1:3 , l] <- coef(lm(value ~ X1 + X2, data = melt(tmp[ , , s])))
            }
            wedges.w[ , 4, l] <- wedges.w[ , 2, l] / wedges.w[ , 3, l]
        }
        saveRDS(wedges.w, file = "./Other-data/Panel-gradients-white.rds")
    }
}

###########################################################################################
wedges.b <- readRDS("./Other-data/Panel-gradients-black.rds")
smoothed.panels.b <- readRDS("./Other-data/Smoothed-panels-black.rds")
wedges.g <- readRDS("./Other-data/Panel-gradients-grey.rds")
wedges.w <- readRDS("./Other-data/Panel-gradients-white.rds")

n <- dim(wedges.b)[[3]]
###########################################################################################
# try linear regression with panel as factor

zz <- subpanels(pw.m.b[,,"150828"])

df <- data.frame(panel = character(), x = double(), y = double(), value = numeric())
for (i in 1:16) {   # upper panels
    df <- rbind(df, cbind(dimnames(zz)[[3]][i], melt(zz[,,i])))
}
for (i in 17:32) {  # rotate lower panels to get same direction of gradient
    df <- rbind(df, cbind(dimnames(zz)[[3]][i], melt(zz[(128:1),(1024:1),1])))
}
colnames(df) <- c("panel", "x", "y", "value")
df <- df[!is.na(df$value),]

array.lm <- lm(value ~ x + y, data = df)
summary(array.lm)


###########################################################################################
# linear regression over distance from centre  at 1023.5, 992.5
z.dist <- merge(x = c(1:1996), y = c(1:1996))
z.dist$z <- sqrt((z.dist$x - 1023.5)^2 + (z.dist$y - 992.5)^2)

# convert to square array to be plotted
circ.dist <- as.matrix(cast(z.dist, x ~ y))

pw.m <- pw.m.b[,,"150828"]

zz <- cbind(melt(pw.m), z.dist)[,c(4,5,6,3)]

# what is relationship between value and distance?
plot(zz$z, zz$value)

circ.lm <- lm(value ~ z, zz)
coef(circ.lm)

# may need to adjust to convert by row - try moving centre spot to see
circ.fitted <- matrix(circ.lm$fitted.values, ncol = 1996)
circ.res <- matrix(circ.lm$residuals, ncol = 1996)

pixel.image(circ.fitted); draw.panels()
pixel.image(circ.res); draw.panels()

# now fit per-panel linear gradient to remainder
circ.res.panels <- subpanels(circ.res)

smoo <- array(dim = c(128, 1024, 32))
for (i in 1:32) {
    smoo[,,i] <- matrix(predict(lm(value ~ X1 + X2, 
                                   data = melt(circ.res.panels[ , , i]),
                                   na.action = na.omit),
                                melt(circ.res.panels[ , , i])),
                        nrow = 128)
}

smoothed <- join.panels(smoo)
smoothed.res <- circ.res - smoothed
pixel.image(smoothed.res, break.levels = sd.levels(pw.m))
hist(pixel.image, breaks = "fd")

pixel.image(pw.m)

scale.hist(smoothed.res, main = "Residuals after smoothing")
sd(smoothed.res)
mad(smoothed.res)

scale.hist(pw.m)
sd(pw.m)
mad(pw.m)

###########################################################################################
# convert coefficients into cross-panel gain
zz[3,1,1]       # BL = 5178.6
zz[3,1004,1]    # TL = 5644.9   - sensor
zz[128,1,1]     # BR = 5657 
zz[128,1004,1]  # TR = 5632.9

zz[3,1004,1] - zz[3,1,1] 
wedges.b[1,"x",1] * 128     # 432.1894
wedges.b[1,"y",1] * 1004    # 337.7277

tmp <- melt(zz[,,1])
tmp <- tmp[!is.na(tmp[,3]),]

lm <- lm(value ~ X1 * X2, data = tmp)
h.grad <- apply(tmp, 2, lm, value ~ X1)

zz[128,1,1] - zz[3,1,1]         # BR to BL: 478.4 
zz[128,1004,1] - zz[3,1004,1]   # TR to TL: -12

zz[3,1004,1] - zz[128,1,1]

y.diffs <- zz[128,,1] - zz[3,,1]

pixel.image(zz[,,1], break.levels = sd.levels(pw.m.b[,,1]))
points(matrix(c(128,1), ncol = 2, byrow = T))

acq.cols <- c("darkblue", "blue", "purple", "violetred", "red", "orangered",
          "orange", "gold", "green", "chartreuse4", "cyan4")

plot(abs(wedges.b[1:16, "x", 1]), type = "o", xaxt = "none", ylim = c(0, 5.5), xlab = "", ylab = "",
     col = acq.cols[1], pch = 20, main = "X-gradient change over panels: upper panels")
axis(1, at = c(1:16), labels = dimnames(wedges.b)[[1]][1:16])
for (i in 2:n) {
    points(abs(wedges.b[1:16, "x", i]), type = "o", col = acq.cols[i], pch = 20)
}
