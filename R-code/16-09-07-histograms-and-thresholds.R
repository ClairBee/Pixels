
library("IO.Pixels"); library("CB.Misc")

# LINE DETECTION ALGORITHM UNSATISFACTORY

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7")

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

####################################################################################################

# FUNCTIONS                                                                                     ####

# find MAD thresholds based on two half-normal densities, conjoined at modal density
asymm.bounds <- function(dat, n = 6) {
    
    zz <- density(dat, n = 65536, na.rm = T)
    
    # break at point of maximum density
    mu <- zz$x[which.max(zz$y)]
    
    # calculate MAD on either side of breakpoint, centred at breakpoint
    c(mu - n * mad(dat[which(dat <= mu)], center = mu),
      mu + n * mad(dat[which(dat > mu)], center = mu))
}

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "", JF = F, ...) {
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "", ...)
    abline(v = asymmetric.mad(dat, n = 6), lty = 2, col = "red")
    abline(v = asymmetric.mad(dat, n = 5), lty = 3, col = "red")
    
    if (JF) {
        abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JohnsonFit(dat[!is.na(dat)])), col = "cyan3", lty = c(2,3,3,2))
    }
}

fit.lm <- function(im, terms = "g ~ b * w") {
    
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    
    # fit linear model to central part of image only (excludes edge effects)
    w.lm <- lm(as.formula(terms), 
               data = df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                             findInterval(df$y, c(40.5, 2008.5)) == 1, ])
    
    df$fv <- predict(w.lm, df)
    
    target <- gsub(" ~.*$", "", terms)
    
    df$res <- eval(parse(text = paste0("df$", target))) - df$fv
    list(df = df, r2 = round(summary(w.lm)$adj.r.squared, 3), rmse = round(summary(w.lm)$sigma, 2))
}

####################################################################################################

# THRESHOLDING BY ASYMMETRIC MAD                                                                ####

fpath <- "./Image-plots/thresholds/"

px <- list()

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "thresholds-", dt, ".pdf"), width = 14)
                     par(mfrow = c(2, 4), mar = c(3,3,3,1))
                     
                     im.b <- pw.m[,,"black", dt]; im.g <- pw.m[,,"grey", dt]
                     res.b <- md7[,,"black", dt]; res.g <- md7[,,"grey", dt]
                     
                     hist.with.boundaries(im.b, title = paste0(dt, " - black values"))
                     hist.with.boundaries(im.g, title = paste0(dt, " - grey values"), xlim = c(0,25000))
                     
                     hist.with.boundaries(res.b, title = paste0(dt, " - black residuals"), xlim = c(-1000,1000))
                     hist.with.boundaries(res.g, title = paste0(dt, " - grey residuals"), xlim = c(-1000,1000))
                     
                     n.val.b <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2)))
                     n.val.g <- length(which(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2)))
                     n.res.b <- length(which(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2)))
                     n.res.g <- length(which(findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     n.all.b <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2) | findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2)))
                     n.all.g <- length(which(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2) | findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     n.all <- length(which(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2) | findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2) |
                                               findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2) | findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2)))
                     
                     px[[dt]] <<- data.frame(n.val.b, n.val.g, n.res.b, n.res.g, n.all.b, n.all.g, n.all)
                     
                     pixel.plot(which(array(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2),
                                            dim = dim(im.b)), arr.ind = T), 
                                main = paste0(dt, " - black values (", n.val.b, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2),
                                            dim = dim(im.g)), arr.ind = T),
                                main = paste0(dt, " - grey values  (", n.val.g, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2),
                                            dim = dim(res.b)), arr.ind = T), 
                                main = paste0(dt, " - black residuals (", n.res.b, ")"))
                     draw.panels(col = "grey")
                     
                     pixel.plot(which(array(findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2),
                                            dim = dim(res.g)), arr.ind = T), 
                                main = paste0(dt, " - grey residuals (", n.res.g, ")"))
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))
                     

df <- rbind.fill(px)
rownames(df) <- dimnames(pw.m)[[4]]

write.csv(df, paste0(fpath, "px-identified.csv"), quote = F)

####################################################################################################

# CLASSIFICATION BY ASYMMETRIC MAD                                                              ####

dt <- "160314"

# get all abnormal pixels, by mean value & median-smoothed residual
bpx <- which(array(findInterval(pw.m[,,"black", dt], asymm.bounds(pw.m[,,"black", dt])) %in% c(0,2) |
                       findInterval(pw.m[,,"grey", dt], asymm.bounds(pw.m[,,"grey", dt])) %in% c(0,2) |
                       findInterval(md7[,,"black", dt], asymm.bounds(md7[,,"black", dt])) %in% c(0,2) |
                       findInterval(md7[,,"grey", dt], asymm.bounds(md7[,,"grey", dt])) %in% c(0,2), dim = dim(pw.m[,,"black", dt])),
             arr.ind = T)

os <- data.frame(x = bpx[,1],
                 y = bpx[,2],
                 b = pw.m[,,"black", dt][bpx],
                 g = pw.m[,,"grey", dt][bpx],
                 res.b = md7[,,"black", dt][bpx],
                 res.g = md7[,,"grey", dt][bpx])


abline(h = asymm.bounds(md7[,,"grey", dt]), col = "green3", lty = 2)
abline(v = asymm.bounds(md7[,,"black", dt]), col = "green3", lty = 2)


# cut data midway between bounds & max values
class.boundaries <- function(dat) {
    
    rng <- c(floor(min(dat, na.rm = T)), ceiling(max(dat, na.rm = T)) + 1)
    inner.bounds <- asymm.bounds(dat)
    
    vb <- inner.bounds[2] + ((rng[2] - inner.bounds[2]) * c(0.25, 0.5, 0.75))

    sort(c(rng, inner.bounds, vb))
}

.smoothScatter(os$res.b, os$res.g, nrpoints = Inf, xlab = "Black residuals", ylab = "Grey residuals", main = paste0(dt))
abline(0,1, col = "red")

abline(v = class.boundaries(md7[,,"black", dt]), col = "skyblue", lty = 2)
abline(h = class.boundaries(md7[,,"grey", dt]), col = "skyblue", lty = 2)

table(findInterval(md7[,,"black", dt], class.boundaries(md7[,,"black", dt])), useNA = "ifany")

os$class.b <- findInterval(os$b, class.boundaries(os$b))
os$class.g <- findInterval(os$g, class.boundaries(os$g))

table("black" = os$class.b, "grey" = os$class.g)


cor(os[!is.na(os$res.b), c("res.b", "g")])
cor(os[!is.na(os$res.g), c("res.g", "b")])





cor(os[!is.na(os$res.b) & os$g < 65535 & os$b > 0, c("res.b", "res.g")])

smoothScatter(os$res.b, os$res.g, nrpoints = Inf)
abline(0,1, col = "red")

smoothScatter(os$b, os$res.b, nrpoints = Inf, xlab = "Values", ylab = "Residuals", main = paste0(dt, " - black"))
abline(a = -5000, b = 1, col = "red")

smoothScatter(os$g, os$res.g, nrpoints = Inf, xlab = "Values", ylab = "Residuals", main = paste0(dt, " - grey"))
abline(a = -17500, b = 1, col = "red")

plot(os, lower.panel = panel.cor)

####################################################################################################

# CHECK LINEARITY ISSUES VS OFFSET ERRORS                                                       ####

# decide whether to include as separate category
# are most nonlinear pixels already identified using residual/extreme-value approach?

wlm <- fit.w.lm(pw.m[,,,dt])

.smoothScatter(wlm$df$g.x, wlm$df$g.y)

hist.with.boundaries(wlm$df$res, xlim = c(-2000,2000))

nlpx <- wlm$df[findInterval(wlm$df$res, asymm.bounds(wlm$df$res)) %in% c(0,2),]
             
.smoothScatter(nlpx$res.b, nlpx$res.g)
abline(0,1)

pixel.plot(nlpx, col = "cyan3")
pixel.plot(bpx, col = "orange")

px <- merge(data.frame(nlpx, type = "nl"), data.frame(os, type = "os"), by = c("x", "y"), all = T)
px$b <- apply(px[,c("b.x", "b.y")], 1, mean, na.rm = T)
px$g <- apply(px[,c("g.x", "g.y")], 1, mean, na.rm = T)

px$res.b <- md7[,,"black", dt][as.matrix(px[,c("x", "y")])]
px$res.g <- md7[,,"grey", dt][as.matrix(px[,c("x", "y")])]

.smoothScatter(px[is.na(px$type.y) & !is.na(px$type.x), c("res.b", "res.g")])
abline(v = class.boundaries(md7[,,"black", dt]), lty = 2)
abline(h = class.boundaries(md7[,,"grey", dt]), lty = 2)

# check number of nonlinear pixels NOT identified as offset errors in all images



####################################################################################################

# LINEARITY VS OFFSET IN ALL IMAGES                                                             ####


os.px <- function(dt) {
    
    dt <- toString(dt)
    which(array(findInterval(pw.m[,,"black", dt], asymm.bounds(pw.m[,,"black", dt])) %in% c(0,2) |
                    findInterval(pw.m[,,"grey", dt], asymm.bounds(pw.m[,,"grey", dt])) %in% c(0,2) |
                    findInterval(md7[,,"black", dt], asymm.bounds(md7[,,"black", dt])) %in% c(0,2) |
                    findInterval(md7[,,"grey", dt], asymm.bounds(md7[,,"grey", dt])) %in% c(0,2), dim = dim(pw.m[,,"black", dt])),
          arr.ind = T)
}

#=========================================================================

bpx <- os.px("loan")

wlm <- fit.lm(pw.m[,,, "loan"], "g ~ b * w")

nlpx <- wlm$df[findInterval(wlm$df$res, asymm.bounds(wlm$df$res)) %in% c(0,2),]

hist.with.boundaries(wlm$df$res, xlim = c(-10000,10000))

.smoothScatter(wlm$df$fv, wlm$df$res, xlab = "fitted", ylab = "residual")

wlm.w <- fit.lm(pw.m[,,, "loan"], "w ~ b * g")
.smoothScatter(wlm.w$df$fv, wlm.w$df$res, xlab = "fitted", ylab = "residual")
hist.with.boundaries(wlm.w$df$res, xlim = c(-10000,10000))

pixel.image(array(wlm$df$res, dim = c(2048, 2048)))
pixel.image(array(wlm.w$df$res, dim = c(2048, 2048)), title = "white fitting")

nlpx.g <- wlm$df[findInterval(wlm$df$res, asymm.bounds(wlm$df$res)) %in% c(0,2),]
nlpx.w <- wlm.w$df[findInterval(wlm.w$df$res, asymm.bounds(wlm.w$df$res)) %in% c(0,2),]

all(wlm$df$y == wlm.w$df$y)
.smoothScatter(wlm$df$res, wlm.w$df$res, xlab = "grey", ylab = "white")

abline(v = asymm.bounds(wlm$df$res), col = "grey")
abline(h = asymm.bounds(wlm.w$df$res), col = "grey")

pdf("./tmp-offset-corrected--loan.pdf")
pixel.image(shading.corrected(pw.m[,,,"loan"]))
draw.outlines(sc.spots)
dev.off()

pixel.plot(nlpx.g)
pixel.plot(nlpx.w)

hist.with.boundaries(wlm$df$res, xlim = c(-2000,2000))

pixel.plot(which(array(wlm$df$res, dim = c(2048, 2048)) < -700, arr.ind = T))
pixel.plot(which(array(wlm$df$res, dim = c(2048, 2048)) > 700, arr.ind = T))

sc.spots <- screen.spots(pw.m[,,"white","loan"])
pixel.image(shading.corrected(pw.m[,,,"loan"]))
draw.outlines(sc.spots)

sc <- shading.corrected(pw.m[,,, "loan"])
sc.px <- which(array(findInterval(sc, asymm.bounds(sc)) %in% c(0,2), dim = c(2048, 2048)), arr.ind = T)

hist.with.boundaries(sc, xlim = c(10000, 30000)) 
pixel.plot(which(array(findInterval(sc, asymm.bounds(sc, n = 5)) %in% c(0,2), dim = c(2048, 2048)), arr.ind = T))

####################################################################################################

# ALL BAD PIXELS IN ALL IMAGES                                                                  ####

fpath <- "./Image-plots/thresholds/"

all.bpx <- list(); jf.bpx <- list()

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     bmp(paste0(fpath, "thresholds-", dt, ".bmp"), height = 1920, width = 2880)
                     par(mfrow = c(4, 6), mar = c(3,3,3,1))
                     bpx <- list()
                     
                     # select images of interest
                     {
                         im.b <- pw.m[,,"black", dt]; im.g <- pw.m[,,"grey", dt]
                         res.b <- md7[,,"black", dt]; res.g <- md7[,,"grey", dt]
                         linear <- fit.lm(pw.m[,,, dt], "g ~ b * w")
                         sc <- shading.corrected(pw.m[,,,dt])
                     }
                     
                     # fit Johnson distriutions
                     {
                         JF <- list(im.b = JohnsonFit(im.b[!is.na(im.b)]),
                                    im.g = JohnsonFit(im.g[!is.na(im.g)]),
                                    res.b = JohnsonFit(res.b[!is.na(res.b)]),
                                    res.g = JohnsonFit(res.g[!is.na(res.g)]),
                                    linear = JohnsonFit(linear$df$res[!is.na(linear$df$res)]),
                                    sc = JohnsonFit(sc[!is.na(sc)]))
                     }
                     
                     # images
                     {
                         pixel.image(im.b, title = paste0(dt, " - black image"))
                         pixel.image(im.g, title = paste0(dt, " - grey image"))
                         pixel.image(res.b, title = paste0(dt, " - black ms residuals"))
                         pixel.image(res.g, title = paste0(dt, " - grey ms residuals"))
                         pixel.image(array(linear$df$res, dim = dim(im.b)), title = paste0(dt, " - linear residuals"))
                         pixel.image(sc, title = paste0(dt, " - shading-corrected"))
                     }
                     
                     # draw histograms with boundaries
                     {
                         hist.with.boundaries(im.b, JF = T, title = paste0(dt, " - black values"))
                         hist.with.boundaries(im.g, JF = T, title = paste0(dt, " - grey values"), xlim = c(0,25000))
                         
                         hist.with.boundaries(res.b, JF = T, title = paste0(dt, " - black residuals"), xlim = c(-1000,1000))
                         hist.with.boundaries(res.g, JF = T, title = paste0(dt, " - grey residuals"), xlim = c(-1000,1000))
                         
                         hist.with.boundaries(linear$df$res, JF = T, title = paste0(dt, " - linear residuals"), xlim = c(-2000,2000))
                         hist.with.boundaries(sc, JF = T, title = paste0(dt, " - shading-corrected"), xlim = median(sc, na.rm = T) + c(-1000,1000))
                     }

                # plot pixels identified by MAD
                     {
                     # black values
                     {
                         bpx$val.b <- which(array(findInterval(im.b, asymm.bounds(im.b)) %in% c(0,2),
                                                       dim = dim(im.b)), arr.ind = T)
                         
                         pixel.plot(bpx$val.b, main = paste0(dt, " - black values (", nrow(bpx$val.b), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # grey values
                     {
                         bpx$val.g <- which(array(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2),
                                                  dim = dim(im.g)), arr.ind = T)
                         pixel.plot(bpx$val.g, main = paste0(dt, " - grey values  (", nrow(bpx$val.g), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # black residuals
                     {
                         bpx$res.b <- which(array(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2),
                                                  dim = dim(res.b)), arr.ind = T)
                         pixel.plot(bpx$res.b, main = paste0(dt, " - black residuals (", nrow(bpx$res.b), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # grey residuals
                     {
                         bpx$res.g <- which(array(findInterval(res.g, asymm.bounds(res.g)) %in% c(0,2),
                                                  dim = dim(res.g)), arr.ind = T)
                         pixel.plot(bpx$res.g, main = paste0(dt, " - grey residuals (", nrow(bpx$res.g), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # linear model residuals
                     {
                         bpx$linear <- which(array(findInterval(linear$df$res, asymm.bounds(linear$df$res)) %in% c(0,2),
                                               dim = dim(im.b)), arr.ind = T)
                         pixel.plot(bpx$linear, main = paste0(dt, " - linear residuals (", nrow(bpx$linear), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # shading-corrected residuals
                     {
                         bpx$sc <- which(array(findInterval(sc, asymm.bounds(sc)) %in% c(0,2),
                                               dim = dim(sc)), arr.ind = T)
                         pixel.plot(bpx$sc, main = paste0(dt, " - shading-corrected residuals (", nrow(bpx$sc), ")"))
                         draw.panels(col = "grey")
                     }}
                     
                # plot pixels identified by Johnson distribution
                     {
                         # black values
                         {
                             jfpx$val.b <- which(array(findInterval(im.b, qJohnson(pnorm(c( -5, 5), 0, 1), JF$im.b)) %in% c(0,2),
                                                      dim = dim(im.b)), arr.ind = T)
                             pixel.plot(jfpx$val.b, main = paste0(dt, " - black values (", nrow(jfpx$val.b), ")"))
                             draw.panels(col = "grey")
                         }
                         
                         # grey values
                         {
                             jfpx$val.g <- which(array(findInterval(im.g, qJohnson(pnorm(c( -5, 5), 0, 1), JF$im.g)) %in% c(0,2),
                                                      dim = dim(im.g)), arr.ind = T)
                             pixel.plot(jfpx$val.g, main = paste0(dt, " - grey values  (", nrow(jfpx$val.g), ")"))
                             draw.panels(col = "grey")
                         }
                         
                         # black residuals
                         {
                             jfpx$res.b <- which(array(findInterval(res.b, qJohnson(pnorm(c( -5, 5), 0, 1), JF$res.b)) %in% c(0,2),
                                                      dim = dim(res.b)), arr.ind = T)
                             pixel.plot(jfpx$res.b, main = paste0(dt, " - black residuals (", nrow(jfpx$res.b), ")"))
                             draw.panels(col = "grey")
                         }
                         
                         # grey residuals
                         {
                             jfpx$res.g <- which(array(findInterval(res.g, qJohnson(pnorm(c( -5, 5), 0, 1), JF$res.g)) %in% c(0,2),
                                                      dim = dim(res.g)), arr.ind = T)
                             pixel.plot(jfpx$res.g, main = paste0(dt, " - grey residuals (", nrow(jfpx$res.g), ")"))
                             draw.panels(col = "grey")
                         }
                         
                         # linear model residuals
                         {
                             jfpx$linear <- which(array(findInterval(linear$df$res, qJohnson(pnorm(c( -5, 5), 0, 1), JF$linear)) %in% c(0,2),
                                                       dim = dim(im.b)), arr.ind = T)
                             pixel.plot(jfpx$linear, main = paste0(dt, " - linear residuals (", nrow(jfpx$linear), ")"))
                             draw.panels(col = "grey")
                         }
                         
                         # shading-corrected residuals
                         {
                             jfpx$sc <- which(array(findInterval(sc, qJohnson(pnorm(c( -5, 5), 0, 1), JF$sc)) %in% c(0,2),
                                                   dim = dim(sc)), arr.ind = T)
                             pixel.plot(jfpx$sc, main = paste0(dt, " - shading-corrected residuals (", nrow(jfpx$sc), ")"))
                             draw.panels(col = "grey")
                         }
                     }
                         
                     all.bpx[[dt]] <<- bpx
                     jf.bpx[[dt]] <<- jfpx
                     dev.off()
                     }))

saveRDS(all.bpx, paste0(fpath, "px-identified.rds"))
saveRDS(jf.bpx, paste0(fpath, "px-identified-jf.rds"))


####################################################################################################

# CATEGORICAL CLASSIFICATION OF BAD PIXELS                                                      ####

# dark, dim, offset, stuck, nonlinear. Distnguish between high-density regions and otherwise

# ratio of black:grey value, black:grey offset
dt <- "130613"

px <- do.call("rbind", all.bpx[[dt]])
df <- data.frame(b = pw.m[,,"black", dt][px],
                 g = pw.m[,,"grey", dt][px],
                 res.b = md7[,,"black", dt][px],
                 res.g = md7[,,"grey", dt][px])
df$ratio.val <- df$b / df$g
df$ratio.res <- df$res.b / df$res.g

.smoothScatter(df$ratio.val, df$ratio.res, xlab = "value ratio", ylab = "residual ratio")

hist(df$ratio.val, breaks = "fd")

# use k-means clustering to divide the two groups
th <- kmeans(df$ratio.val, centers = 2)

abline(v = median(min(df$ratio.val[th$cluster == which.max(th$centers)]),
                  max(df$ratio.val[th$cluster == which.min(th$centers)])),
       col = "red")

# check: does this occur in all images?
fpath <- "./Image-plots/stuck-vs-offset-px/"

invisible(lapply(names(all.bpx),
                 function(dt) {
                     bmp(paste0(fpath, "stuck-vs-offset-px-", dt, ".bmp"), width = 1440, height = 720)
                     par(mfrow = c(1,2))
                     
                     px <- do.call("rbind", all.bpx[[dt]])
                     df <- data.frame(b = pw.m[,,"black", dt][px],
                                      g = pw.m[,,"grey", dt][px],
                                      res.b = md7[,,"black", dt][px],
                                      res.g = md7[,,"grey", dt][px])
                     df$diff <- df$g - df$b
                     rv <- df$b / df$g
                     df$ratio.res <- df$res.b / df$res.g
                     
                     rv <- rv[!is.na(rv)]
                     
                     th <- kmeans(rv, centers = 2)
                     th.cut <- median(min(rv[th$cluster == which.max(th$centers)]),
                                      max(rv[th$cluster == which.min(th$centers)]))
                     
                     hist(rv, breaks = "fd", xlim = c(0,1.2), xlab = "", ylab = "", main = paste0(dt, " - black:grey ratio"))
                     abline(v = th.cut, col = "red")
                     abline(v = th$centers, col = "green3", lty = 2)
                     
                     pixel.plot(which(pw.m[,,"black",dt] / pw.m[,,"grey",dt] > th.cut, arr.ind = T),
                                main = paste0("Ratio above ", round(th.cut, 2)))
                     points(which(pw.m[,,"black",dt] / pw.m[,,"grey",dt] > 0.95, arr.ind = T), col = "red", pch = 15, cex = 0.3)
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))

.smoothScatter(df[abs(df$ratio.res) > 10, c("res.b", "res.g")])

####################################################################################################



####################################################################################################

# AREAS OF EXTREMELY HIGH DENSITY                                                               ####

fpath <- "./Image-plots/thresholds/"
mad.bpx <- readRDS(paste0(fpath, "px-identified.rds"))
jf.bpx <- readRDS(paste0(fpath, "px-identified-JF.rds"))

dt <- "130613"

# convert px to binary image, calculate density in small area (11x11?)
dense.px <- function(px, th = 0.9, area = 11) {
    px.im <- array(0, dim = c(2048, 2048))
    px.im[px] <- 1
    
    # define kernel
    k <-  matrix(rep(1 / area^2, area^2), ncol = area)
    
    # convolve image with kernel
    px.density <- r2m(focal(m2r(px.im), k))
    
    # return coordinates exceeding threshold
    which(px.density > th, arr.ind = T)
}

dp <- dense.px(mad.bpx[[dt]]$val.b)

# consider picking up all dense regions in which density > 0.1, as long as associated with region > 0.5
invisible(lapply(names(mad.bpx),
                 function(dt) {
                     bmp(paste0(fpath, "spatial/density-", dt, ".bmp"), width = 2400)
                     par(mfrow = c(1,5), mar = c(2,2,1,1))
                     
                     pixel.plot(mad.bpx[[dt]]$val.b, main = paste0(dt, " - black"))
                     points(dense.px(mad.bpx[[dt]]$val.b, th = 0.1), col = adjustcolor("gold", alpha = 0.4), pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$val.b, th = 0.5), col = adjustcolor("red", alpha = 0.4), pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$val.b, th = 0.9), col = adjustcolor("magenta3", alpha = 0.4), pch = 15, cex = 0.4)
                     
                     pixel.plot(mad.bpx[[dt]]$val.g, main = paste0(dt, " - grey"))
                     points(dense.px(mad.bpx[[dt]]$val.g, th = 0.1), col = "gold", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$val.g, th = 0.5), col = "red", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$val.g, th = 0.9), col = adjustcolor("magenta3", alpha = 0.4), pch = 15, cex = 0.4)
                     
                     pixel.plot(mad.bpx[[dt]]$res.b, main = paste0(dt, " - black residuals"))
                     points(dense.px(mad.bpx[[dt]]$res.b, th = 0.1), col = "gold", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$res.b, th = 0.5), col = "red", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$res.b, th = 0.9), col = adjustcolor("magenta3", alpha = 0.4), pch = 15, cex = 0.4)
                     
                     pixel.plot(mad.bpx[[dt]]$res.g, main = paste0(dt, " - grey residuals"))
                     points(dense.px(mad.bpx[[dt]]$res.g, th = 0.1), col = "gold", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$res.g, th = 0.5), col = "red", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$res.g, th = 0.9), col = adjustcolor("magenta3", alpha = 0.4), pch = 15, cex = 0.4)
                     
                     pixel.plot(mad.bpx[[dt]]$linear, main = paste0(dt, " - linear residuals"))
                     points(dense.px(mad.bpx[[dt]]$linear, th = 0.1), col = "gold", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$linear, th = 0.5), col = "red", pch = 15, cex = 0.4)
                     points(dense.px(mad.bpx[[dt]]$linear, th = 0.9), col = adjustcolor("magenta3", alpha = 0.4), pch = 15, cex = 0.4)
                     
                     dev.off()
                 }))

pixel.plot(dense.px(mad.bpx[[dt]]$val.b, th = 0.1), col = "gold")
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.6), col = "orange", pch = 15, cex = 0.4)
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.7), col = "red", pch = 15, cex = 0.4)
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.8), col = "darkred", pch = 15, cex = 0.4)
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.9), col = "black", pch = 15, cex = 0.4)

pixel.plot(mad.bpx[[dt]]$val.b, xlim = c(1800, 2048), ylim = c(0,200))
zz <- dense.px(mad.bpx[[dt]]$val.b, th = 0.05)

points(dense.px(mad.bpx[[dt]]$val.b, th = 0.5), col = adjustcolor("gold", alpha = 0.3), pch = 15, cex = 0.4)
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.05), col = adjustcolor("green3", alpha = 0.3), pch = 15, cex = 0.4)

draw.outlines(dense.px(mad.bpx[[dt]]$val.b, th = 0.6), col = "orange", lwd = 2)



hh <- density(px.ppp$"val.b")
plot(hh)
plot(px.ppp$"val.b")

# look at values of dense pixels
d.im.b <- array(dim = dim(pw.m[,,"black", dt]))
d.im.b[dense.px(mad.bpx[[dt]]$val.b, th = 0.1)] <- pw.m[,,"black", dt][dense.px(mad.bpx[[dt]]$val.b, th = 0.1)]

# identify dense areas for exclusion:
# areas of density > 0.1, contiguous with

dense.regions <- function(px, th.l = 0.05, th.u = 0.5, area = 11, dilate.by = 21) {
    
    px.im <- array(0, dim = c(2048, 2048))
    px.im[px] <- 1
    
    # define kernel
    k <-  matrix(rep(1 / area^2, area^2), ncol = area)
    
    # convolve image with kernel
    px.density <- r2m(focal(m2r(px.im), k))
    th.density <- (px.density > th.l) * 1
    
    # get clumps of higher-density pixels
    cc <- clump(m2r(th.density))
    
    cand <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                       id = getValues(cc)[!is.na(getValues(cc))])
    cand$density <- px.density[as.matrix(cand[,1:2])]
    
    blocks <- ddply(cand, .(id), max.d = max(density), summarise)

    if ((dilate.by == 0) | is.na(dilate.by)) {
        return(as.matrix(cand[cand$id %in% blocks$id[blocks$max.d >= th.u], c("x", "y")]))
    } else {
        th.density[as.matrix(cand[cand$id %in% blocks$id[blocks$max.d < th.u], c("x", "y")])] <- 0
        th.density[is.na(th.density)] <- 0
        
        sk <- shapeKernel(c(dilate.by, dilate.by), type = "disc")
        zz <- dilate(th.density, sk)
        
        return(which(zz == 1, arr.ind = T))
    }
}

pixel.plot(dense.regions(mad.bpx[[dt]]$val.b))
points(dense.px(mad.bpx[[dt]]$val.b, th = 0.1), col = "gold", pch = 15, cex = 0.3)

invisible(lapply(names(mad.bpx),
                 function(dt) {
                     bmp(paste0(fpath, "spatial/density-", dt, ".bmp"), width = 2400)
                     par(mfrow = c(1,5), mar = c(2,2,1,1))
                     
                     pixel.plot(mad.bpx[[dt]]$val.b, main = paste0(dt, " - black"), cex = 0.3)
                     points(dense.regions(mad.bpx[[dt]]$val.b), col = "cyan3", pch = 15, cex = 0.3)
                     
                     pixel.plot(mad.bpx[[dt]]$val.g, main = paste0(dt, " - grey"))
                     points(dense.regions(mad.bpx[[dt]]$val.g), col = "cyan3", pch = 15, cex = 0.3)
                    
                     pixel.plot(mad.bpx[[dt]]$res.b, main = paste0(dt, " - black residuals"))
                     points(dense.regions(mad.bpx[[dt]]$res.b), col = "cyan3", pch = 15, cex = 0.3)
                     
                     pixel.plot(mad.bpx[[dt]]$res.g, main = paste0(dt, " - grey residuals"))
                     points(dense.regions(mad.bpx[[dt]]$res.g), col = "cyan3", pch = 15, cex = 0.3)
                     
                     pixel.plot(mad.bpx[[dt]]$linear, main = paste0(dt, " - linear residuals"))
                     points(dense.regions(mad.bpx[[dt]]$linear), col = "cyan3", pch = 15, cex = 0.3)
                     
                     dev.off()
                 }))

####################################################################################################

# CLASSIFICATION                                                                                ####

fpath <- "./Image-plots/thresholds/"
mad.bpx <- readRDS(paste0(fpath, "px-identified.rds"))

classify.all.px <- function(dt) {
    bpx <- unique(do.call("rbind", mad.bpx[[dt]][1:5]))
    dp <- dense.regions(bpx)
    wlm <- fit.lm(pw.m[,,, dt], "g ~ b * w")
    
    # create df with all necessary pixel attributes
    
    px <- data.frame(bpx,
                     b.v = pw.m[,,"black", dt][bpx],
                     g.v = pw.m[,,"grey", dt][bpx],
                     w.v = pw.m[,,"white", dt][bpx],
                     b.res = md7[,,"black", dt][bpx],
                     g.res = md7[,,"grey", dt][bpx])
    
    px <- merge(px, setNames(wlm$df[, c("x", "y", "res")], nm = c("x", "y", "lm.res")), by = c(1:2), all.x = T)
    
    #---------------------------------------------------------------------------------
    
    # absolute dark & hot pixels
    px$type[px$w.v < 15000] <- "dark"
    px$type[px$b.v == 65535] <- "hot"
    
    # ordinal classification by grey values
    # hist(pw.m[,,"grey", dt], breaks = "fd")
    # abline(v = class.boundaries(pw.m[,,"grey", dt])[-6], col = "red")
    
    px$type[is.na(px$type)] <- c("dim", NA, "s.bright", "bright", "v.bright", "v.bright")[findInterval(px$g.v, class.boundaries(pw.m[,,"grey", dt]))][is.na(px$type)]
    
    # further ordinal classification by grey residuals after median-smoothing
    px$type[is.na(px$type)] <- c("l.dim", NA, "l.bright", "bright", "v.bright", "v.bright")[findInterval(px$g.res, class.boundaries(md7[,,"grey", dt]))][is.na(px$type)]
    
    # yet further ordinal classification by black values & smoothed residuals
    px$type[is.na(px$type)] <- c("dim", NA, "s.bright", "bright", "v.bright", "v.bright")[findInterval(px$b.v, class.boundaries(pw.m[,,"black", dt]))][is.na(px$type)]
    px$type[is.na(px$type)] <- c("l.dim", NA, "l.bright", "bright", "v.bright", "v.bright")[findInterval(px$b.res, class.boundaries(md7[,,"black", dt]))][is.na(px$type)]
    
    # finally, pixels with nonlinear response
    px$type[is.na(px$type) & px$lm.res > asymm.bounds(wlm$df$res)[2]] <- "nl.bright"
    px$type[is.na(px$type) & px$lm.res < asymm.bounds(wlm$df$res)[1]] <- "nl.dim"
    
    px$type <- factor(px$type, levels = c("hot", "dark", "v.bright", "bright", "s.bright", "l.bright",
                                          "nl.bright", "nl.dim", "l.dim", "dim"))
    
    if (nrow(dp) > 0) px <- merge(px, data.frame(dp, dense = T), by = c(1,2), all.x = T)
    
    return(px)
}

zz <- invisible(sapply(names(mad.bpx), classify.all.px, simplify = F))

saveRDS(zz, "./02_Objects/px-classified.rds")

#---------------------------------------------------------------------------------
pixel.plot(px[!is.na(px$type),], col = "lightgrey", cex = 0.3)
points(px[is.na(px$type),], pch = 15, cex = 0.3)
points(px[px$dense,], pch = 15, cex = 0.3, col = adjustcolor("cornflowerblue", alpha = 0.3))

table(px$type, useNA = "ifany")
#---------------------------------------------------------------------------------

do.call("rbind", lapply(lapply(zz, "[", "type"), table, useNA = "ifany"))

# plot black, grey & white values of each class against one another. Patterns?

####################################################################################################

# SPATIAL DISTRIBUTION OF BAD PIXELS                                                            ####

fpath <- "./Image-plots/thresholds/"
mad.bpx <- readRDS(paste0(fpath, "px-identified.rds"))
jf.bpx <- readRDS(paste0(fpath, "px-identified-JF.rds"))

dt <- "130613"
px.ppp <- lapply(mad.bpx[[dt]], function(px) ppp(px[,1], px[,2], c(1,2048), c(1,2048)))


bmp(paste0(fpath, "spatial/spatial-dist-", dt, ".bmp"))
par(mfrow = c(3,5), mar = c(2,2,2,1))

plot(envelope(ppp(mad.bpx[[dt]]$val.b[,1], mad.bpx[[dt]]$val.b[,2], c(1,2048), c(1,2048)),
              Gest, nsim = 99, nrank = 2), main = "Black values")
plot(envelope(ppp(mad.bpx[[dt]]$val.g[,1], mad.bpx[[dt]]$val.g[,2], c(1,2048), c(1,2048)),
              Gest, nsim = 99, nrank = 2), main = "Grey values")
plot(envelope(ppp(mad.bpx[[dt]]$res.b[,1], mad.bpx[[dt]]$res.b[,2], c(1,2048), c(1,2048)),
              Gest, nsim = 99, nrank = 2), main = "Black residuals")
plot(envelope(ppp(mad.bpx[[dt]]$res.g[,1], mad.bpx[[dt]]$res.g[,2], c(1,2048), c(1,2048)),
              Gest, nsim = 99, nrank = 2), main = "Grey residuals")
plot(envelope(ppp(mad.bpx[[dt]]$linear[,1], mad.bpx[[dt]]$linear[,2], c(1,2048), c(1,2048)),
              Gest, nsim = 99, nrank = 2), main = "Nonlinear response")
par(mfrow = c(1,1))



####################################################################################################

# GET BAD PIXELS BY JOHNSON THRESHOLDING FOR COMPARISON                                         ####

jf.bpx <- list()

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     bmp(paste0(fpath, "thresholds-JF-", dt, ".bmp"), height = 1440, width = 2880)
                     par(mfrow = c(3, 6), mar = c(3,3,3,1))
                     bpx <- list()
                     
                     # select images of interest
                     {
                         im.b <- pw.m[,,"black", dt]; im.g <- pw.m[,,"grey", dt]
                         res.b <- md7[,,"black", dt]; res.g <- md7[,,"grey", dt]
                         linear <- fit.lm(pw.m[,,, dt], "g ~ b * w")
                         sc <- shading.corrected(pw.m[,,,dt])
                     }
                     

                     
                     # images
                     {
                         pixel.image(im.b, title = paste0(dt, " - black image"))
                         pixel.image(im.g, title = paste0(dt, " - grey image"))
                         pixel.image(res.b, title = paste0(dt, " - black ms residuals"))
                         pixel.image(res.g, title = paste0(dt, " - grey ms residuals"))
                         pixel.image(array(linear$df$res, dim = dim(im.b)), title = paste0(dt, " - linear residuals"))
                         pixel.image(sc, title = paste0(dt, " - shading-corrected"))
                     }
                     
                     # draw histograms with boundaries
                     {
                         hist.with.boundaries(im.b, title = paste0(dt, " - black values"), JF = T)
                         hist.with.boundaries(im.g, title = paste0(dt, " - grey values"), JF = T, xlim = c(0,25000))
                         
                         hist.with.boundaries(res.b, title = paste0(dt, " - black residuals"), JF = T, xlim = c(-1000,1000))
                         hist.with.boundaries(res.g, title = paste0(dt, " - grey residuals"), JF = T, xlim = c(-1000,1000))
                         
                         hist.with.boundaries(linear$df$res, title = paste0(dt, " - linear residuals"), JF = T, xlim = c(-2000,2000))
                         hist.with.boundaries(sc, title = paste0(dt, " - shading-corrected"), JF = T, xlim = median(sc, na.rm = T) + c(-1000,1000))
                     }
                     
                     # plot pixels identified by each approach
                     # black values
                     {
                         bpx$val.b <- which(array(findInterval(im.b, qJohnson(pnorm(c( -5, 5), 0, 1), JF$im.b)) %in% c(0,2),
                                                  dim = dim(im.b)), arr.ind = T)
                         pixel.plot(bpx$val.b, main = paste0(dt, " - black values (", nrow(bpx$val.b), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # grey values
                     {
                         bpx$val.g <- which(array(findInterval(im.g, qJohnson(pnorm(c( -5, 5), 0, 1), JF$im.g)) %in% c(0,2),
                                                  dim = dim(im.g)), arr.ind = T)
                         pixel.plot(bpx$val.g, main = paste0(dt, " - grey values  (", nrow(bpx$val.g), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # black residuals
                     {
                         bpx$res.b <- which(array(findInterval(res.b, qJohnson(pnorm(c( -5, 5), 0, 1), JF$res.b)) %in% c(0,2),
                                                  dim = dim(res.b)), arr.ind = T)
                         pixel.plot(bpx$res.b, main = paste0(dt, " - black residuals (", nrow(bpx$res.b), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # grey residuals
                     {
                         bpx$res.g <- which(array(findInterval(res.g, qJohnson(pnorm(c( -5, 5), 0, 1), JF$res.g)) %in% c(0,2),
                                                  dim = dim(res.g)), arr.ind = T)
                         pixel.plot(bpx$res.g, main = paste0(dt, " - grey residuals (", nrow(bpx$res.g), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # linear model residuals
                     {
                         bpx$linear <- which(array(findInterval(linear$df$res, qJohnson(pnorm(c( -5, 5), 0, 1), JF$linear)) %in% c(0,2),
                                                   dim = dim(im.b)), arr.ind = T)
                         pixel.plot(bpx$linear, main = paste0(dt, " - linear residuals (", nrow(bpx$linear), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     # shading-corrected residuals
                     {
                         bpx$sc <- which(array(findInterval(sc, qJohnson(pnorm(c( -5, 5), 0, 1), JF$sc)) %in% c(0,2),
                                               dim = dim(sc)), arr.ind = T)
                         pixel.plot(bpx$sc, main = paste0(dt, " - shading-corrected residuals (", nrow(bpx$sc), ")"))
                         draw.panels(col = "grey")
                     }
                     
                     jf.bpx[[dt]] <<- bpx
                     dev.off()
                 }))

saveRDS(jf.bpx, paste0(fpath, "px-identified-JF.rds"))

####################################################################################################

# PLOTS FOR EXPLANATION OF THRESHOLD CHANGE FOR JULIA                                           ####

fpath <- "./Notes/00_thresholds/fig/"

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "", JF = F, ...) {
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "", ...)
    abline(v = asymm.bounds(dat, n = 6), lty = 2, col = "red")

    if (JF) {
        abline(v = qJohnson(pnorm(c(-5, 5), 0, 1), JohnsonFit(dat[!is.na(dat)])), col = "cyan3", lty = 2)
    }
}


# 141009 black residuals - far too many values picked up
res.b <- md7[,,"black", "141009"]
pdf(paste0(fpath, "hist-black-res-141009.pdf")); {
    par(mar = c(2,2,1,1))
    hist.with.boundaries(res.b, xlim = c(-1000,1000), JF = T)
    dev.off()
}  
pdf(paste0(fpath, "plot-JF-black-res-141009.pdf")); {
    par(mar = c(2,2,1,1))
    pixel.plot(which(array(findInterval(res.b, qJohnson(pnorm(c(-5, 5), 0, 1), JohnsonFit(res.b[!is.na(res.b)]))) %in% c(0,2),
                           dim = dim(res.b)), arr.ind = T), col = "cyan3")
    dev.off()
}
pdf(paste0(fpath, "plot-MAD-black-res-141009.pdf")); {
    par(mar = c(2,2,1,1))
    pixel.plot(which(array(findInterval(res.b, asymm.bounds(res.b)) %in% c(0,2),
                           dim = dim(res.b)), arr.ind = T), col = "red")
    dev.off()
}  

# 130613 grey values - dark lines not picked up
im.g <- pw.m[,,"grey", "130613"]
pdf(paste0(fpath, "hist-grey-values-130613.pdf")); {
    par(mar = c(2,2,1,1))
    hist.with.boundaries(im.g, xlim = c(0,30000), JF = T)
    dev.off()
}  
pdf(paste0(fpath, "plot-JF-grey-values-130613.pdf")); {
    par(mar = c(2,2,1,1))
    pixel.plot(which(array(findInterval(im.g, qJohnson(pnorm(c(-5, 5), 0, 1), JohnsonFit(im.g[!is.na(im.g)]))) %in% c(0,2),
                           dim = dim(im.g)), arr.ind = T), col = "cyan3")
    dev.off()
}
pdf(paste0(fpath, "plot-MAD-grey-values-130613.pdf")); {
    par(mar = c(2,2,1,1))
    pixel.plot(which(array(findInterval(im.g, asymm.bounds(im.g)) %in% c(0,2),
                           dim = dim(im.g)), arr.ind = T), col = "red")
    dev.off()
}  

####################################################################################################

# PERSISTENCE THROUGH REFURBISHMENT?                                                            ####

pre <- data.frame(unique(do.call("rbind", all.bpx$"131122")), pre = T)
post <- data.frame(unique(do.call("rbind", all.bpx$"140128")), post = T)

px.match <- merge(pre, post, by = c(1:2), all = T)
px.match[is.na(px.match)] <- F

table(px.match[,c("pre", "post")])
stayed1 <- as.matrix(px.match[px.match$pre & px.match$post,1:2])

pw.m[,,"black", "131122"][stayed1]
pw.m[,,"grey", "131122"][stayed1]

md7[,,"black", "131122"][stayed1]
md7[,,"grey", "131122"][stayed1]

pw.m[,,"black", "140128"][stayed1]
pw.m[,,"grey", "140128"][stayed1]



#==============================================================================

pre <- data.frame(unique(do.call("rbind", all.bpx$"140129")), pre = T)
post <- data.frame(unique(do.call("rbind", all.bpx$"141009")), post = T)

px.match <- merge(pre, post, by = c(1:2), all = T)
px.match[is.na(px.match)] <- F

table(px.match[,c("pre", "post")])
stayed <- as.matrix(px.match[px.match$pre & px.match$post,1:2])

pw.m[,,"black", "140129"][stayed]
pw.m[,,"grey", "140129"][stayed]

pw.m[,,"black", "141009"][stayed]
pw.m[,,"grey", "141009"][stayed]

md7[,,"black", "140129"][stayed]
md7[,,"black", "141009"][stayed]

md7[,,"grey", "140129"][stayed]
md7[,,"grey", "141009"][stayed]

pixel.plot(stayed)
pixel.plot(stayed1)

####################################################################################################

# BEHAVIOUR OF ROW DEFECTS IN LOAN PANEL                                                        ####

o.plot(pw.m[, 1025, "black", "loan"] - pw.m[, 1026, "black", "loan"])
abline(h = 0, col = "red")

o.plot(pw.m[, 77, "black", "loan"] - pw.m[, 78, "black", "loan"])
abline(h = 0, col = "red")

####################################################################################################

####################################################################################################

# DISCONTINUED                                                                                  ####

# JOHNSON DISTRIBUTIONS OF ALL BLACK & GREY IMAGES

{fpath <- "./Image-plots/histograms/"

normal.Johnson.QQ <- function (data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99), title = "Normal Q-Q plot", ...) {
    
    data <- data[!is.na(data)]
    m <- mean(data)
    s <- mad(data)
    jf <- JohnsonFit(data, moment = "quant")
    
    plot(qnorm(quantiles, m, s), quantile(data, quantiles), 
         pch = 20, asp = T, ylab = "Observed quantile", xlab = "Fitted quantile", 
         main = title, col = adjustcolor("seagreen", alpha = 0.3), ...)
    abline(h = quantile(data, grid.quantiles), col = "lightseagreen", lty = 2)
    abline(v = qnorm(grid.quantiles, m, s), col = "lightseagreen", lty = 2)
    
    points(qJohnson(quantiles, jf), quantile(data, quantiles), pch = 20)
    abline(v = qJohnson(grid.quantiles, jf), col = "skyblue", lty = 2)
    abline(v = qnorm(grid.quantiles, m, s), col = "skyblue", lty = 2)
    
    abline(0, 1, col = "red", lty = 2)
    
    legend("topleft", cex = 0.7, title = "Gridlines: ", legend = paste0("Q", grid.quantiles * 100), lty = 2, col = "skyblue", box.col = "white")
    legend("bottomright", cex = 0.7, title = "Points", col = c("seagreen", "black"), pch = 20, box.col = "white",
           legend = c("Normal (MAD)", "Johnson"))
}

# histograms with fit, Q-Q plot & identified defects: black images
invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "hist-black-", dt, ".pdf"))
                     par(mfrow = c(2,2), mar = c(3,3,3,1))
                     
                     JF <- JohnsonFit(pw.m[,,"black", dt][!is.na(pw.m[,,"black", dt])])
                     
                     hist(pw.m[,,"black", dt], breaks = "fd", prob = T, xlim = c(0,10000), main = paste0(dt, " - black"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), sd(pw.m[,,"black", dt], na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), mad(pw.m[,,"black", dt], na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = median(pw.m[,,"black", dt], na.rm = T) + c(-6,-5,5,6) * mad(pw.m[,,"black", dt], na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topright", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     hist(pw.m[,,"black", dt], breaks = "fd", prob = T, ylim = c(0,0.0005), xlim = c(0,10000), main = paste0(dt, " - black (cropped)"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), sd(pw.m[,,"black", dt], na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(pw.m[,,"black", dt], na.rm = T), mad(pw.m[,,"black", dt], na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = median(pw.m[,,"black", dt], na.rm = T) + c(-6,-5,5,6) * mad(pw.m[,,"black", dt], na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     legend("topright", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     normal.Johnson.QQ(pw.m[,,"black", dt], title = paste0("Q-Q plot - ", dt, " black"))
                     
                     med <- median(pw.m[,,"black", dt], na.rm = T)
                     sp <- mad(pw.m[,,"black", dt], na.rm = T)
                     
                     pixel.plot(which(pw.m[,,"black", dt] > med + 6 * sp | pw.m[,,"black", dt] < med - 6 * sp, arr.ind = T), col = "gold", cex = 0.3, main = paste0("Defects - ", dt))
                     points(which(pw.m[,,"black", dt] > qJohnson(pnorm(5, 0, 1), JF) | pw.m[,,"black", dt] < qJohnson(pnorm(-5, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3, col = "cyan3")
                     points(which(pw.m[,,"black", dt] > qJohnson(pnorm(6, 0, 1), JF) | pw.m[,,"black", dt] < qJohnson(pnorm(-6, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3)
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))


# histograms with fit, Q-Q plot & identified defects: grey images
invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     pdf(paste0(fpath, "hist-grey-", dt, ".pdf"))
                     par(mfrow = c(2,2), mar = c(3,3,3,1))
                     
                     f.im <- im <- pw.m[,,"grey", dt]
                     
                     zz <- density(im, n = 65536, na.rm = T)
                     mu <- zz$x[which.max(zz$y)]
                     
                     # filter image to remove dark pixels
                     f.im <- f.im[which(f.im > 10000)]
                     JF <- JohnsonFit(f.im)
                     
                     hist(im, breaks = "fd", prob = T, xlim = c(0, 26000), main = paste0(dt, " - grey"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), sd(f.im, na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), mad(f.im, na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = mu + c(-6,-5,5,6) * mad(im, na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topleft", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     hist(im, breaks = "fd", prob = T, ylim = c(0,0.0005), xlim = c(0, 26000), main = paste0(dt, " - grey (cropped)"), xlab = "", ylab = "")
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), sd(f.im, na.rm = T)), col = "red", lwd = 1)
                     lines(0:65535, dnorm(0:65535, mean(f.im, na.rm = T), mad(f.im, na.rm = T)), col = "orange", lwd = 1)
                     lines(0:65535, dJohnson(0:65535, JF), col = "cyan3", lwd = 1)
                     abline(v = mu + c(-6,-5,5,6) * mad(im, na.rm = T), col = "orange", lty = c(3,2,2,3))
                     abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JF), col = "cyan3", lty = c(3,2,2,3))
                     abline(h = 0.0005, col = "grey", lty = 3)
                     legend("topleft", col = c("red", "orange", "cyan3"), lwd = 1, legend = c("Normal - SD", "Normal - MAD", "Johnson"), bty = "n")
                     
                     normal.Johnson.QQ(f.im, title = paste0("Q-Q plot - ", dt, " grey"))
                     
                     med <- median(f.im, na.rm = T)
                     sp <- mad(f.im, na.rm = T)
                     
                     pixel.plot(which(im > mu + 6 * sp | im < mu - 6 * sp, arr.ind = T), col = "gold", cex = 0.3, main = paste0("Defects - ", dt))
                     points(which(im > qJohnson(pnorm(5, 0, 1), JF) | im < qJohnson(pnorm(-5, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3, col = "cyan3")
                     points(which(im > qJohnson(pnorm(6, 0, 1), JF) | im < qJohnson(pnorm(-6, 0, 1), JF), arr.ind = T), pch = 15, cex = 0.3)
                     draw.panels(col = "grey")
                     
                     dev.off()
                 }))
}
####################################################################################################

# HALF-NORMAL DISTRIBUTION FOR HISTOGRAMS 
{
library(VGAM)

dat <- pw.m[,,"black", "160430"]
hist(dat, breaks = c(0:65535), xlim = c(4000,6000), prob = T)
zz <- density(dat, n = 65536, na.rm = T)

mu <- zz$x[which.max(zz$y)]
abline(v = mu, col = "red")

# break at point of maximum density
dat.l <- abs(dat[which(dat <= mu)] - mu)
dat.u <- dat[which(dat > mu)] - mu

sig.l <- sqrt(mean(dat.l^2))
sig.u <- sqrt(mean(dat.u^2))

hist(dat.l, breaks = "fd", xlim = c(0,5000), prob = T)
lines(0:5000, dhalfnorm(0:5000, sd2theta(sig.l)), col = "red", lwd = 2)

hist(dat.u, breaks = "fd", xlim = c(0,5000), prob = T)
lines(0:5000, dhalfnorm(0:5000, sd2theta(sig.u)), col = "red", lwd = 2)}

# images with strong upward drift are poorly fitted by half-normal. Hey ho.

####################################################################################################
