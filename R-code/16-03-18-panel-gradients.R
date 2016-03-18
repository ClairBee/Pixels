
library("IO.Pixels")
library(reshape)

pw.m.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
pw.sd.b <- readRDS("./Other-data/Pixelwise-sds-black.rds")
n <- dim(pw.m.b)[[3]]

wedges.b <- array(dim = c(32, 4, n), 
                dimnames = list(apply(cbind(c(rep("U", 16), rep("L", 16)), 
                                            rep(c(1:16), 2)), 1, paste, collapse = ""),
                                c("offset", "x", "y", "ratio"),
                                dimnames(pw.m.b)[[3]]))

# takes less than a minute to run for 11 images
for (l in 1:n) {
    tmp <- subpanels(pw.m.b[ , , l])
    for (s in 1:32) {
        wedges.b[s, 1:3 , l] <- coef(lm(value ~ X1 + X2, data = melt(tmp[ , , s])))
    }
    wedges.b[ , 4, l] <- wedges.b[ , 2, l] / wedges.b[ , 3, l]
}

saveRDS(wedges.b, file = "./Other-data/Panel-gradients-black.rds")

# gradient change over time
{
# x-gradient: use absolute value (sign switches for upper & lower panels)
o.plot(abs(wedges.b.b[1,"x",]), xaxt = "none", ylim = c(0, 5.5),
       main = "LM gradient change over time: X")
axis(1, at = c(1:n), labels = dimnames(wedges.b)[[3]])
for (i in 2:32) {
    o.plot(abs(wedges.b[i, "x", ]), add = T, col = c(rep("black", 16), rep("blue", 16))[i])
}

o.plot(abs(wedges.b[1,"y",]), xaxt = "none", ylim = c(0, 2),
       main = "LM gradient change over time: Y")
axis(1, at = c(1:n), labels = dimnames(wedges.b)[[3]])
for (i in 2:32) {
    o.plot(abs(wedges.b[i, "y", ]), add = T, col = c(rep("black", 16), rep("blue", 16))[i])
}

o.plot(wedges.b[1,4,], xaxt = "none", main = "LM gradient change over time: diagonal",
       ylim = c(1, 100))
axis(1, at = c(1:n), labels = dimnames(wedges.b)[[3]])
for (i in 2:32) {
    o.plot(wedges.b[i, 4, ], add = T, col = c(rep("black", 16), rep("blue", 16))[i])
}
}


# gradient change over subpanels
o.plot(wedges.b[, 4, 1], xaxt = "none", ylim = c(1, 100),
       main = "LM gradient across panels: diagonal")
abline(v = 16.5, col = "red", lwd = 2)
axis(1, at = c(1:32), labels = dimnames(wedges.b)[[1]])
for (i in 2:n) {
    o.plot(wedges.b[, 4, i], add = T)
}

o.plot(abs(wedges.b[, "x", 1]), xaxt = "none", ylim = c(0, 5.5),
       main = "LM gradient across panels: X")
abline(v = 16.5, col = "red", lwd = 2)
axis(1, at = c(1:32), labels = dimnames(wedges.b)[[1]])
for (i in 2:n) {
    o.plot(abs(wedges.b[, "x", i]), add = T)
}

o.plot(abs(wedges.b[, "y", 1]), xaxt = "none", ylim = c(0, 2),
       main = "LM gradient across panels: Y")
abline(v = 16.5, col = "red", lwd = 2)
axis(1, at = c(1:32), labels = dimnames(wedges.b)[[1]])
for (i in 2:n) {
    o.plot(abs(wedges.b[, "y", i]), add = T)
}
