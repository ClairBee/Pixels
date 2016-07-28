
load.pixel.means()
load.pixel.sds()

dim(pw.m)
dim(pw.sd)

pw.snr <- pw.m / pw.sd
pw.snr[is.infinite(pw.snr)] <- NA

hist(pw.snr[,,"grey", "160430"], breaks = "fd", xlim = c(0,1000), prob = T, ylim = c(0,0.0005))
mean(pw.snr[,,"grey", "160430"], na.rm = T)
lines(dnorm(0:1000, mean(pw.snr[,,"grey", "160430"], na.rm = T), sd(pw.snr[,,"grey", "160430"], na.rm = T)),
            col = "red", lwd = 2)
lines(dJohnson(0:1000, JohnsonFit(pw.snr[,,"grey", "160430"][!is.na(pw.snr[,,"grey", "160430"])])), col = "cyan3", lwd = 2)
abline(v = mean(pw.snr[,,"grey", "160430"], na.rm = T), col = "orange")

plot(which(pw.snr[,,"grey", "160430"] > qJohnson(0.99995, JF.g), arr.ind = T), cex = 0.4, pch = 15, col = "red", xlim = c(0,1996), ylim = c(0,1996))
points(which(pw.snr[,,"grey", "160430"] < qJohnson(0.00005, JF.g), arr.ind = T), cex = 0.4, pch = 15, col = "blue", xlim = c(0,1996), ylim = c(0,1996))

JF.g <- JohnsonFit(pw.snr[,,"grey", "160430"][!is.na(pw.snr[,,"grey", "160430"])])
JF.g
abline(v = qJohnson(c(0.0005, 0.9995), JF.g), col = "gold")

