
# remake all plots (after changing MAD thresholding approach)

pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7")
l.res <- load.objects("./02_Objects/linear-res/", otype = "l-res")
spot.models <- read.csv("./Other-data/Gaussian-spots.csv", row.names = 1)

fpath <- "./Image-plots/thresholds/"

####################################################################################################

th.px <- function(im, centre = NA) {
    i <- array(findInterval(im, asymmetric.mad(im, fix.centre = centre)) %in% c(0,2), dim = dim(im))
    which(i == 1, arr.ind = T)
}

invisible(lapply(dimnames(pw.m)[[4]],
                 function(dt) {
                     bmp(paste0(fpath, "thresholds-revised-", dt, ".bmp"), height = 960 * 3, width = 960 * 6, pointsize = 28); {
                         par(mfrow = c(3, 6), mar = c(2,2,3,1))
                         im.b <- pw.m[,,"black", dt]; im.g <- pw.m[,,"grey", dt]
                         md.b <- md7[,,"black", dt]; md.g <- md7[,,"grey", dt]
                         linear <- l.res[,,dt]
                         if (is.na(spot.models$fm.rmse[dt])) {
                             gs.res <- im.g - im.b - gaussian.spot.mat(unlist(spot.models[dt, 1:5]))
                         } else {
                             gs.res <- im.g - im.b - gaussian.spot.mat(unlist(setNames(spot.models[dt, 7:11], spot.models[dt, 1:5])))
                         }
                         
                         pixel.image(im.b, title = paste0(dt, " - black image"))
                         pixel.image(im.g, title = paste0(dt, " - grey image"))
                         pixel.image(md.b, title = paste0(dt, " - black residuals"))
                         pixel.image(md.g, title = paste0(dt, " - grey residuals"))
                         pixel.image(linear, title = paste0(dt, " - linear residuals"))
                         pixel.image(gs.res, title = paste0(dt, " - spot residuals"))
                         
                         hist(im.b, breaks = "fd", xlab = "", ylab = "", xlim = c(0,10000),
                              main = paste0("Mode: ", round(modal.density(im.b),0), "; MAD ",
                                            paste(round(abs(modal.density(im.b) - asymmetric.mad(im.b)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(im.b), col = "cyan3", lwd = 2)
                               
                         hist(im.g, breaks = "fd", xlab = "", ylab = "", xlim = c(0,65535),
                              main = paste0("Mode: ", round(modal.density(im.g),0), "; MAD ",
                                            paste(round(abs(modal.density(im.g) - asymmetric.mad(im.g)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(im.g), col = "cyan3", lwd = 2)
                         
                         hist(md.b, breaks = "fd", xlab = "", ylab = "", xlim = c(-1000,1000),
                              main = paste0("MAD ", paste(round(abs(asymmetric.mad(md.b, fix.centre = 0)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(md.b, fix.centre = 0), col = "cyan3", lwd = 2)
                         
                         hist(md.g, breaks = "fd", xlab = "", ylab = "", xlim = c(-1000,1000),
                              main = paste0("MAD ", paste(round(abs(asymmetric.mad(md.g, fix.centre = 0)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(md.g, fix.centre = 0), col = "cyan3", lwd = 2)
                         
                         hist(linear, breaks = "fd", xlab = "", ylab = "", xlim = c(-2000, 2000), 
                              main = paste0("MAD ", paste(round(abs(asymmetric.mad(linear, fix.centre = 0)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(linear, fix.centre = 0), col = "cyan3", lwd = 2)
                         
                         hist(gs.res, breaks = "fd", xlab = "", ylab = "", xlim = c(-2000,2000), 
                              main = paste0("MAD ", paste(round(abs(asymmetric.mad(gs.res, fix.centre = 0)),0), collapse = ", ")))
                         abline(v = asymmetric.mad(gs.res, fix.centre = 0), col = "cyan3", lwd = 2)
                         
                         pixel.plot(th.px(im.b), main = paste0(nrow(th.px(im.b)), " pixels"))
                         pixel.plot(th.px(im.g), main = paste0(nrow(th.px(im.g)), " pixels"))
                         pixel.plot(th.px(md.b, centre = 0), main = paste0(nrow(th.px(md.b, centre = 0)), " pixels"))
                         pixel.plot(th.px(md.g, centre = 0), main = paste0(nrow(th.px(md.g, centre = 0)), " pixels"))
                         pixel.plot(th.px(linear, centre = 0), main = paste0(nrow(th.px(linear, centre = 0)), " pixels"))
                         pixel.plot(th.px(gs.res, centre = 0), main = paste0(nrow(th.px(gs.res, centre = 0)), " pixels"))
                         
                         dev.off()
                     }
                 }))

abs(modal.density(l.res[,,"loan"]) - asymmetric.mad(l.res[,,"loan"], fix.centre = NA, n = 1))
abs(asymmetric.mad(l.res[,,"loan"], fix.centre = 0, n = 1))

