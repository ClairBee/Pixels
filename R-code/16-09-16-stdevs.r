
# standard deviation vs value for all images
library("IO.Pixels"); library("CB.Misc")

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
pw.sd <- load.objects("./02_Objects/stdevs/", otype = "pwsd")

fpath <- "./Image-plots/SDs/"

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

####################################################################################################

# FUNCTIONS                                                                                     ####

hist.with.boundaries <- function(dat, xlim = c(0,10000), title = "", JF = F, ...) {
    
    hist(dat, breaks = "fd", xlim = xlim, main = title, xlab = "", ylab = "", ...)
    abline(v = asymmetric.mad(dat, n = 6), lty = 2, col = "red")
    abline(v = asymmetric.mad(dat, n = 5), lty = 3, col = "red")
    
    if (JF) {
        abline(v = qJohnson(pnorm(c(-6, -5, 5, 6), 0, 1), JohnsonFit(dat[!is.na(dat)])), col = "cyan3", lty = c(2,3,3,2))
    }
}

extreme.px <- function(im, n = 6) {
    which(array(findInterval(im, asymmetric.mad(im, n = n)) %in% c(0,2),
                dim = dim(im)), arr.ind = T)
}

####################################################################################################

# PLOTS                                                                                         ####

px <- readRDS("./Image-plots/thresholds/px-identified.rds")

invisible(lapply(dimnames(pw.sd)[[4]],
                 function(dt) {
                     
                     bmp(paste0(fpath, "sd-vs-mean-", dt, ".bmp"), height = 960 * 3, width = 960 * 3); {
                         par(mfrow = c(3,3))
                         
                         .smoothScatter(pw.m[,,"black", dt], pw.sd[,,"black", dt], nbin = 256, 
                                        xlab = "mean value", ylab = "sd", main = paste0(dt, " - black pixel behaviour"))
                         .smoothScatter(pw.m[,,"grey", dt], pw.sd[,,"grey", dt], nbin = 256, 
                                        xlab = "mean value", ylab = "sd", main = paste0(dt, " - grey pixel behaviour"))
                         .smoothScatter(pw.m[,,"white", dt], pw.sd[,,"white", dt], nbin = 256,
                                        xlab = "mean value", ylab = "sd", main = paste0(dt, " - white pixel behaviour"))
                         
                         hist.with.boundaries(pw.sd[,,"black", dt], xlim = c(0,1000), title = "black SDs")
                         hist.with.boundaries(pw.sd[,,"grey", dt], xlim = c(0,1000), title = "grey SDs")
                         hist.with.boundaries(pw.sd[,,"white", dt], xlim = c(0,1000), title = "white SDs")
                         
                         pixel.plot(extreme.px(pw.sd[,,"black", dt]), col = "red", main = "px with high SD - black")
                            points(unique(do.call("rbind", px[[dt]])), pch = 15, cex = 0.4)
                         
                         pixel.plot(extreme.px(pw.sd[,,"grey", dt]), col = "red", main = "px with high SD - grey")
                            points(unique(do.call("rbind", px[[dt]])), pch = 15, cex = 0.4)
                         
                         pixel.plot(extreme.px(pw.sd[,,"white", dt]), col = "red", main = "px with high SD - white")
                            points(unique(do.call("rbind", px[[dt]])), pch = 15, cex = 0.4)
                            
                        dev.off()
                   }
                 }))

####################################################################################################

# THRESHOLDING                                                                                  ####

invisible(lapply(dimnames(pw.sd)[[4]],
                 function(dt) {
                     
                     px[[dt]]$sd.b <<- extreme.px(pw.sd[,,"black", dt])
                     px[[dt]]$sd.gw <<-  unique(rbind(extreme.px(pw.sd[,,"grey", dt]),
                                                     extreme.px(pw.sd[,,"white", dt])))
                 }))

saveRDS(px, "./Image-plots/thresholds/px-identified.rds")

####################################################################################################

# COMPARE SD PX WITH VALUE PX                                                                   ####

pw.px <- lapply(lapply(px[7:21], "[", 1:6), 
                function(lst) data.frame(unique(do.call("rbind", lst)), value = T))
sd.b.px <- lapply(lapply(px[7:21], "[", 7),
                function(lst) data.frame(lst, sd.b = T))
sd.gw.px <- lapply(lapply(px[7:21], "[", 8),
                   function(lst) data.frame(lst, sd.gw = T))

matched.px <- invisible(sapply(names(sd.gw.px),
                               function(dt) {
                                   merge(merge(pw.px[[dt]], sd.b.px[[dt]], by = c(1:2), all = T),
                                         sd.gw.px[[dt]], by = c(1:2), all = T)
                               }, simplify = F))
matched.px <- lapply(matched.px,
                      function(df)  {
                          df[is.na(df)] <- F
                          df
                      })

lapply(matched.px, function(pp) table(pp[,c("value", "sd.b")], useNA = "ifany"))

# very small numbers of NA, NA (high grey/white SD, not unusual black SD or value)
# numbers of high-black-SD pixels seem to increase in line with bright pixels.

# plot SD vs value for high-black-SD pixels

dt <- "160314"

pixel.plot(sd.b.px[[dt]], main = paste0("high-SD black pixels - ", dt))
points(matched.px[[dt]][matched.px[[dt]]$sd.gw & !matched.px[[dt]]$sd.b & !matched.px[[dt]]$value,1:2],
       pch = 15, cex = 0.4, col = "red")

sd.b <- lapply(matched.px, function(df) {as.matrix(df[df$sd.b & !df$value & !df$sd.gw,1:2])})
sd.gw <- lapply(matched.px, function(df) {as.matrix(df[!df$sd.b & !df$value & df$sd.gw,1:2])})


.smoothScatter(pw.m[,,"black", dt], pw.m[,,"grey", dt])
points(pw.m[,,"black", dt][sd.b[[dt]]], pw.m[,,"grey", dt][sd.b[[dt]]], col = "cyan3", pch = 15, cex = 0.3)

sc <- array(apply(pw.m, 4, shading.corrected), dim = c(2048, 2048, 21), dimnames = dimnames(pw.m[,,"black", ]))

.smoothScatter(pw.sd[,,"black", dt], sc[,,dt])
points(pw.sd[,,"black", dt][sd.b[[dt]]], sc[,,dt][sd.b[[dt]]], col = "cyan3", pch = 15, cex = 0.3,
       xlab = "Pixelwise SD", ylab = "Shading-corrected value")

invisible(lapply(names(sd.b),
                 function(dt) {
                     
                     bmp(paste0(fpath, "black-sd-vs-SC-", dt, ".bmp"),
                         height = 960, width = 960)
                     .smoothScatter(pw.sd[,,"black", dt], sc[,,dt], main = paste0(dt, " - black SD vs shading correction"),
                                    xlab = "Pixelwise SD", ylab = "Shading-corrected value")
                     points(pw.sd[,,"black", dt][sd.b[[dt]]], sc[,,dt][sd.b[[dt]]], col = "cyan3", pch = 15, cex = 0.3)
                 
                     dev.off()
                 }))

invisible(lapply(names(sd.gw),
                 function(dt) {
                     
                     bmp(paste0(fpath, "white-sd-vs-SC-", dt, ".bmp"),
                         height = 960, width = 960)
                     .smoothScatter(pw.sd[,,"white", dt], sc[,,dt], main = paste0(dt, " - white SD vs shading correction"),
                                    xlab = "Pixelwise SD", ylab = "Shading-corrected value")
                     points(pw.sd[,,"white", dt][sd.gw[[dt]]], sc[,,dt][sd.gw[[dt]]], col = "cyan3", pch = 15, cex = 0.3)
                     
                     dev.off()
                 }))

####################################################################################################

# HIGH WHITE/GREY STANDARD DEVIATION IN LOAN PANEL                                              ####
dt <- "loan"

.smoothScatter(pw.sd[,,"grey", dt], pw.sd[,,"white", dt], xlab = "grey", ylab = "white", xlim = c(0,2500))

pixel.plot(extreme.px(pw.sd[,,"black", dt]))
pixel.plot(extreme.px(pw.sd[,,"grey", dt]))
pixel.plot(extreme.px(pw.sd[,,"white", dt]))

pixel.plot(extreme.px(pw.sd[,,"grey", dt]))
rect(1500, 0, 2048, 200, border = "red")

# compare SD in the damaged corner region
.smoothScatter(pw.sd[1500:2048, 0:200, "grey", dt], pw.sd[1500:2048, 0:200, "white", dt], xlab = "grey", ylab = "white")

# compare SD of area in other corners (assume should be roughly same)
plot(pw.sd[1:548, 1848:2048, "grey", dt], pw.sd[1:548, 1848:2048, "white", dt], xlab = "grey", ylab = "white",
     pch = 20)
points(pw.sd[1500:2048, 0:200, "grey", dt], pw.sd[1500:2048, 0:200, "white", dt], pch = 20, col = "blue")

.smoothScatter(pw.sd[,, "grey", dt], pw.sd[,, "white", dt], xlab = "grey", ylab = "white", xlim = c(0,1000), ylim = c(0,1000))
points(pw.sd[1500:2048, 1:200, "grey", dt], pw.sd[1500:2048, 1:200, "white", dt], col = "cyan3", pch = ".")

# white SD seems normal - grey is high

pixel.plot(p)

hist(pw.m[,,"black", "160314"], breaks = "fd", xlim = c(0,10000))
hist(pw.m[,,"black", "160314"][extreme.px(pw.sd[,,"black", "160314"])], breaks= "fd", add = T, border = "red", col = "red")
abline(v = asymmetric.mad(pw.m[,,"black", "160314"]))

invisible(lapply(dimnames(pw.sd)[[4]],
                 function(dt) {
                     bmp(paste0(fpath, "value-hists-by-SD-", dt, ".bmp"),
                         height = 480 * 2, width = 480 * 3)
                     par(mfrow = c(2,3))
                     
                     sdpx.b <- extreme.px(pw.sd[,,"black", dt])
                     sdpx.g <- extreme.px(pw.sd[,,"grey", dt])
                     sdpx.w <- extreme.px(pw.sd[,,"white", dt])
                     
                     hist(pw.m[,,"black", dt], breaks = "fd", ylim = c(0,10000), col = "black", xlab = "", ylab = "",
                          main = paste0(dt, " - black values, split by SD: ", nrow(sdpx.b)))
                     hist(pw.m[,,"black", dt][sdpx.b], breaks= "fd", add = T, border = "red", col = "red")
                     abline(v = asymmetric.mad(pw.m[,,"black", dt]), col = "orange")
                     
                     hist(pw.m[,,"grey", dt], breaks = "fd", ylim = c(0,10000), col = "black", main = paste0(dt, " - grey values, split by SD: ", nrow(sdpx.g)), xlab = "", ylab = "")
                     hist(pw.m[,,"grey", dt][sdpx.g], breaks= "fd", add = T, border = "red", col = "red")
                     abline(v = asymmetric.mad(pw.m[,,"grey", dt]), col = "orange")
                     
                     hist(pw.m[,,"white", dt], breaks = "fd", ylim = c(0,10000), col = "black",  main = paste0(dt, " - white values, split by SD: ", nrow(sdpx.w)), xlab = "", ylab = "")
                     hist(pw.m[,,"white", dt][sdpx.w], breaks= "fd", add = T, border = "red", col = "red")
                     abline(v = asymmetric.mad(pw.m[,,"white", dt]), col = "orange")
                     
                     pixel.plot(sdpx.b, cex = 0.5)
                        points(pw.px[[dt]], cex = 0.5, col = "cyan3", pch = 15)
                     pixel.plot(sdpx.g, cex = 0.5)
                        points(pw.px[[dt]], cex = 0.5, col = "cyan3", pch = 15)
                     pixel.plot(sdpx.w, cex = 0.5)
                        points(pw.px[[dt]], cex = 0.5, col = "cyan3", pch = 15)
                     
                     dev.off()
                 }))

pixel.image(sc[,,"loan"])
image.scale(sc[,,"loan"], range(sc[,,"loan"], na.rm = T), col = sd.colours(), 
             breaks = sd.levels(sc[,,"loan"]))
