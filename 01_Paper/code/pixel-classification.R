
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

unique.panels <- c("130613", "140128", "160430", "loan", "MCT225")
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7")
pw.m <- load.objects("./02_Objects/images/", otype = "pwm")

sc <- array(apply(pw.m, 4, shading.corrected), dim = dim(pw.m)[-3], dimnames = dimnames(pw.m)[-3])

# need to get remaining pixelwise SDs (for old data)
# check SDs of dead pixels

####################################################################################################

# METHOD FROM MANUAL                                                                            ####

# get offset-corrected images
os <- pw.m[,,c("grey", "white"), ] - pw.m[,,c("black", "black"), ] # offset-corrected images

# get median-differenced shading corrected images
{
    zz <- apply(sc, 3, function(im) r2m(focal(m2r(im), matrix(rep(1, 81), nrow = 9), fun = median)))
    saveRDS(scmd.99, "./02_Objects/sc-med-diffs/sc-md-9x9.rds")
}
scmd.99 <- readRDS("./02_Objects/sc-med-diffs/sc-md-9x9.rds")


# signal sensitivity
{
    ss <- list(up.b.g = which(os[,,"grey", "160430"] > 1.5 * median(os[,,"grey", "160430"], na.rm = T), arr.ind = T),
               up.b.w = which(os[,,"white", "160430"] > 1.5 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T),
               n.g = which(os[,,"grey", "160430"] < 10000 | os[,,"white", "160430"] < 15000, arr.ind = T),
               up.d.g = which(os[,,"grey", "160430"] < 0.45 * median(os[,,"grey", "160430"], na.rm = T), arr.ind = T),
               up.d.w = which(os[,,"white", "160430"] < 0.45 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T))
               
    spx <- rbind.fill(lapply(names(ss[sapply(ss, nrow) > 0]), function(tp) data.frame(ss[[tp]], type = tp)))
}

# noise
{
    np <- list(bn.g = which(pw.sd[,,"grey", "160430"] > 6 * median(pw.sd[,,"grey", "160430"], na.rm = T), arr.ind = T),
               bn.w = which(pw.sd[,,"white", "160430"] > 6 * median(pw.sd[,,"white", "160430"], na.rm = T), arr.ind = T),
               dn.b = which(pw.sd[,,"black", "160430"] > 6 * median(pw.sd[,,"black", "160430"], na.rm = T), arr.ind = T))
    npx <- rbind.fill(lapply(names(np[sapply(np, nrow) > 0]), function(tp) data.frame(np[[tp]], type = tp)))
}

# uniformity
{
    uu <- list(gl.u = which(sc[,,"160430"] > 1.02 * median(sc[,,"160430"], na.rm = T) | sc[,,"160430"] < 0.98 * median(sc[,,"160430"], na.rm = T), arr.ind = T),
               loc.u = which(sc[,,"160430"] > 1.01 * scmd.99[,,"160430"] | sc[,,"160430"] < 0.99 * scmd.99[,,"160430"], arr.ind = T))
    upx <- rbind.fill(lapply(names(uu[sapply(uu, nrow) > 0]), function(tp) data.frame(uu[[tp]], type = tp)))
}

# combine into single map
{
    # merge & remove duplicates
    px <- rbind(spx, npx, upx)
    px$type <- ordered(px$type, levels = c("n.g", "up.b.g", "up.d.g", "up.d.w", "dn.b", "gl.u", "loc.u"))
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    
    pixel.plot(px, cex = 0.5, col = c("black", "red", "blue", "blue", "green3", "orange", NA)[px$type])
}


# and now, to apply over all images
{
    bpx <- sapply(dimnames(pw.m)[[4]],
                 function (dt) {
                     
                     # only run noise check if pw sd exists 
                     if (dt %in% dimnames(pw.sd)[[4]]) {
                         nn <- list(bright.noise = which(pw.sd[,,"grey", dt] > 6 * median(pw.sd[,,"grey", dt], na.rm = T) | 
                                                             pw.sd[,,"white", dt] > 6 * median(pw.sd[,,"white", dt], na.rm = T), arr.ind = T),
                                    dark.noise = which(pw.sd[,,"black", dt] > 6 * median(pw.sd[,,"black", dt], na.rm = T), arr.ind = T))
                         nn <- rbind.fill(lapply(names(nn[sapply(nn, nrow) > 0]), function(tp) data.frame(nn[[tp]], type = tp)))
                     } else {
                         nn <- NULL
                     }
                     
                     ss <- list(no.gain = which(os[,,"grey", dt] < 10000 | os[,,"white", dt] < 15000, arr.ind = T),
                                up.bright = which(os[,,"grey", dt] > 1.5 * median(os[,,"grey", dt], na.rm = T) |
                                                   os[,,"white", dt] > 1.5 * median(os[,,"white", dt], na.rm = T), arr.ind = T),
                                up.dim = which(os[,,"grey", dt] < 0.45 * median(os[,,"grey", dt], na.rm = T) |
                                                   os[,,"white", dt] < 0.45 * median(os[,,"white", dt], na.rm = T), arr.ind = T),
                                gl.unif = which(sc[,,dt] > 1.02 * median(sc[,,dt], na.rm = T) | sc[,,dt] < 0.98 * median(sc[,,dt], na.rm = T), arr.ind = T),
                                local.unif = which(sc[,,dt] > 1.01 * scmd.99[,,dt] | sc[,,dt] < 0.99 * scmd.99[,,dt], arr.ind = T))
                    
                     ss <- rbind.fill(lapply(names(ss[sapply(ss, nrow) > 0]), function(tp) data.frame(ss[[tp]], type = tp)))
                     px <- rbind(nn, ss)
                     px$type <- ordered(px$type, levels = c("no.gain", "up.bright", "up.dim", "dark.noise", "bright.noise", "gl.unif", "local.unif"))
                     px <- px[order(px$type),]
                     px <- px[!duplicated(px[,1:2]),]
                     }, simplify = F)
}

saveRDS(bpx, paste0(fpath, "official-bad-px.rds"))

px.summ <- rbind.fill.matrix(sapply(bpx, function(px) t(as.matrix(table(px$type))), simplify = F))
row.names(px.summ) <- dimnames(pw.m)[[4]]

bpx <- readRDS(paste0(fpath, "official-bad-px.rds"))

table(bpx$"141009"$type)

ssp <- screen.spots(pw.m[,,"white", "141009"], enlarge = T)

####################################################################################################

# MEDIAN-DIFFERENCING OVER INDIVIDUAL IMAGES                                                    ####

JF.threshold <- function(im, cut.probs = pnorm(c(-6,6), 0, 1)) {
    
    dat <- im[!is.na(im)]
    JF <- JohnsonFit(dat)
    
    cut.points <-  qJohnson(cut.probs, JF)
    
    rbind(which(im < min(cut.points), arr.ind = T),
          which(im > max(cut.points), arr.ind = T))
}

# write JF-threshold function, compare to 'new categories' by SD thresholding,
# make final decision, write up.

# histograms of all observed dark residuals after median-differencing
hist(md7[,,"black",], breaks = "fd", ylim = c(0,100), xlab = "Observed GV", ylab = "Frequency", main = "Residuals after median-smoothing in black images")

JF <- JohnsonFit(md7[,,"black", ][!is.na(md7[,,"black", ])])
abline(v = qJohnson(pnorm(c(-6,6), 0, 1), JF), col = "red")

sum(md7[,,"black",] > qJohnson(pnorm(6, 0, 1), JF), na.rm = T)
sum(md7[,,"black",] < qJohnson(pnorm(-6, 0, 1), JF), na.rm = T)    

# replace dark lines with 0 in residual dark images
md7.adj <- abind(sapply(names(dl), 
                        function(dt) {
                            b <- replace(md7[,,"black", dt], as.matrix(dl[[dt]][,1:2]), 0)
                            g <- replace(md7[,,"grey", dt], as.matrix(dl[[dt]][,1:2]), 0)
                            w <- replace(md7[,,"white", dt], as.matrix(dl[[dt]][,1:2]), 0)
                            abind(b, g, w, along = 3, new.names = dimnames(md7[,,,dt]))
                        }, simplify = F),
                 along = 4, new.names = dimnames(md7))
# alternatively: try setting all dark pixels to 0?

# compare results with different cutpoints
JF.px.b.5 <- apply(md7.adj[,,"black",], 3, JF.threshold, cut.probs = pnorm(c(-5,5), 0, 1))
JF.px.b.6 <- apply(md7.adj[,,"black",], 3, JF.threshold, cut.probs = pnorm(c(-6,6), 0, 1))
JF.px.g <- apply(md7.adj[,,"grey",], 3, JF.threshold, cut.probs = pnorm(c(-6,6), 0, 1))

sapply(JF.px.b.5, nrow)
sapply(JF.px.b.6, nrow)
sapply(JF.px.g, nrow)

hist(md7.adj[,,"black", "160430"], breaks = "fd", ylim = c(0,100))

six.sig <- sapply(names(JF.px.g),
                  function(dt) {
                      px <- rbind(data.frame(JF.px.g[[dt]], type = "g"),
                                  data.frame(JF.px.b.6[[dt]], type = "b"))
                      px[!duplicated(px[,1:2]),]
                  })


# check fit of Johnson distribution
Johnson.QQ(md7[,,"black", "140128"])

# histogram of adjusted residuals with thresholds marked - black
hist(md7[,,"black",-3], breaks = "fd", ylim = c(0, 100), xlab = "Observed GV", ylab = "Frequency", main = "Residuals after median-smoothing in black images")
abline(v = qJohnson(pnorm(c(-5,5), 0, 1), JohnsonFit(md7[,,"black", ][!is.na(md7[,,"black", ])])), col = "red")

# histogram of adjusted residuals with thresholds marked - grey
hist(md7[,,"grey",], breaks = "fd", ylim = c(0, 100), xlab = "Observed GV", ylab = "Frequency", main = "Residuals after median-smoothing in grey images")
abline(v = qJohnson(pnorm(c(-5,5), 0, 1), JohnsonFit(md7[,,"grey", ][!is.na(md7[,,"grey", ])])), col = "red")

JF.px.b <- apply(md7[,,"black",], 3, JF.threshold, cut.probs = pnorm(c(-5,5), 0, 1))

sapply(JF.px.b, nrow)

ol <- offset.lines(md7[,,"grey", "MCT225"])
table(ol[,c("x", "type")])


Johnson.QQ(md7[,,"black", 5])
lines(-30000:60000, dJohnson(-30000:60000, JF), col = "cyan3")
abline(v = qJohnson(pnorm(c(-6,6), 0, 1), JF), col = "red")
           
sum(md7[,,"grey",] > qJohnson(pnorm(5, 0, 1), JF), na.rm = T)
sum(md7[,,"grey",] < qJohnson(pnorm(-5, 0, 1), JF), na.rm = T)

####################################################################################################

# NONLINEARITY                                                                                  ####

####################################################################################################

# MEDIAN-DIFFERENCING OVER SHADING CORRECTION                                                   ####

####################################################################################################

# COMPARE RESULTS                                                                               ####

# check for constant offset in any bad pixels
{
    
}

# check correlation between pixelwise SD and median difference (per image/power)
{
    
}

