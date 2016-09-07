
# power in vs brightness out

# where power settings need to be changed to obtain same mean/median, what does operator use as target value?
# how does this relate to local defects in the image?

# example: 160705 had higher power settings than 160430
pixel.image(pw.m[,,"white", "160430"])
pixel.image(pw.m[,,"white", "160705"], break.levels = sd.levels(pw.m[,,"white", "160430"]))

hist(pw.m[,,"white", "160430"], breaks = "fd", col = "black")
hist(pw.m[,,"white", "MCT225"], breaks = "fd", add = T, border = adjustcolor("cyan", alpha = 0.4))
hist(pw.m[,,"white", "160705"], breaks = "fd", add = T, border = adjustcolor("red", alpha = 0.4))


points(median(pw.m[,,"white", "160430"], na.rm = T), -600, pch = 20)
points(median(pw.m[,,"white", "160705"], na.rm = T), -600, pch = 20, col = adjustcolor("red", alpha = 0.4))

h.160430 <- hist(pw.m[,,"white", "160430"], breaks = "fd", plot = F)
h.MCT225 <- hist(pw.m[,,"white", "MCT225"], breaks = "fd", plot = F)
h.160705 <- hist(pw.m[,,"white", "160705"], breaks = "fd", plot = F)

plot(h.160430$mids, h.160430$counts, type = "s", ylim = c(0,60000))
points(h.MCT225$mids, h.MCT225$counts, type = "s", col = "blue")
points(h.160705$mids, h.160705$counts, type = "s", col = "red")


####################################################################################################

# fit Gaussian spots to all grey images
# fastest approach: fit constrained model, check for x0 or y0 on boundary, fit free model if so
models <- invisible(lapply(dimnames(pw.m)[[4]],
                       function(dt) {
                           zz <- gaussian.spot.ls(pw.m[,,"grey", dt] - pw.m[,,"black", dt],
                                                  c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                                                  x0.target = c(768, 1280), y0.target = c(768, 1280))
                           rmse <- sqrt(zz$value / sum(!is.na(pw.m[,,"grey", dt] - pw.m[,,"black", dt])))
                           
                           # check for x0, y0 on boundary
                           if (zz$par["x0"] %in% c(768, 1280) | zz$par["y0"] %in% c(768, 1280)) {
                               
                               # if found, fit unconstrained model
                               fm <- gaussian.spot.ls(pw.m[,,"grey", dt] - pw.m[,,"black", dt],
                                                       c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
                               names(fm$par) <- apply(cbind("fm.", names(fm$par)), 1, paste, collapse = "")
                               fm.rmse <- sqrt(fm$value / sum(!is.na(pw.m[,,"grey", dt] - pw.m[,,"black", dt])))
                               df <- data.frame(acq = toString(dt), t(zz$par), rmse = rmse, 
                                                t(fm$par), fm.rmse = fm.rmse, 
                                                stringsAsFactors = F)
                           } else {
                               # otherwise, return parameters & RMSE for constrained model only
                               df <- data.frame(acq = toString(dt), t(zz$par), rmse = rmse, 
                                                stringsAsFactors = F)
                           }
                           df
                       }))

write.csv(rbind.fill(models), "./Other-data/Gaussian-spots.csv", row.names = F)

#================================================================================
mm <- read.csv("./Other-data/Gaussian-spots.csv", row.names = 1)

# temporary function to create pixel image
plot.spot <- function(dt) {
    dt <- toString(dt)
    pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt])
    contour(1:2048, 1:2048, gaussian.spot.mat(unlist(mm[dt, 1:5])), add = T, drawlabels = F, lty = 1, lwd = 2)
    rect(768, 768, 1280, 1280)
    points(mm[dt, "x0"], mm[dt, "y0"], pch = 15, cex = 1.5)
    
    if (!is.na(mm[dt, "fm.rmse"])) {
        contour(1:2048, 1:2048, gaussian.spot.mat(unlist(setNames(mm[dt, 7:11], nm = names(mm)[1:5]))), 
                add = T, drawlabels = F, lty = 2, lwd = 2)
        points(mm[dt, "fm.x0"], mm[dt, "fm.y0"], pch = 17, cex = 1.5)
        title(main = paste0(dt, "; RMSE ", round(mm[dt, "rmse"]), "; free RMSE ", round(mm[dt, "fm.rmse"])))
    } else {
        title(main = paste0(dt, "; RMSE ", round(mm[dt, "rmse"])))
    }
    cat(paste0(dt, " plotting complete \n"))
}

invisible(lapply(rownames(mm),
                 function(nm) {
                     bmp(paste0("./Image-plots/gaussian-spots/spot-", nm, ".bmp"), height = 2048, width = 2048, pointsize = 20)
                     par(mar = c(2,2,3,1))
                     plot.spot(nm)
                     dev.off()
                 }))

# temporary function to plot residuals
dt <- "160430"
plot(which(abs(pw.m[,,"grey", dt] - pw.m[,,"black", dt] - gaussian.spot.mat(unlist(mm[dt, 1:5]))) > 188 * 2, arr.ind = T), pch = 15, cex = 0.4)


# create histogram of images in which spot was fitted without contraint
qq <- invisible(lapply(row.names(mm)[is.na(mm$fm.rmse)],
                       function(dt) {
                          hist(pw.m[,,"grey", dt], breaks = "fd", add = T, border = NA, col = adjustcolor("skyblue", alpha = 0.4)) 
                       }))

ccols <- c("blue", "red", "green3", "orange", "magenta3", "purple", "cyan3", "skyblue")

plot(0, type = "n", xlim = c(0,65535), ylim = range(pretty(unlist(sapply(qq, "[", "counts")))), xlab = "", ylab = "")
invisible(lapply(qq, function(hh) points(hh$mids, hh$counts, type = "s", col = sample(ccols, 1))))

row.names(mm)[!is.na(mm$fm.rmse)]
hist(pw.m[,,"grey", "140129"], breaks = "fd", add = T, border = NA, col = adjustcolor("green3", alpha = 0.4)) 
hist(pw.m[,,"grey", "151015"], breaks = "fd", add = T, border = NA, col = adjustcolor("orange", alpha = 0.4)) 
hist(pw.m[,,"grey", "160705"], breaks = "fd", add = T, border = NA, col = adjustcolor("magenta3", alpha = 0.4)) 
hist(pw.m[,,"grey", "MCT225"], breaks = "fd", add = T, border = NA, col = adjustcolor("red", alpha = 0.4)) 

# rather convoluted approach. Maybe try checking the SD of each image instead?
sds <- apply(pw.m[,,"grey",], 3, sd, na.rm = T)
plot(sds, pch = 15, col = c("black", "red")[(!is.na(mm$fm.rmse)) + 1])
# works for 'main sequence' images, but not for loan & MCT225 panels.

# could also try fitting spot after adjusting vertical offset per subpanel across midline - reduce the offset value as far as possible?

####################################################################################################

# TRY FITTING SPOT TO UPPER & LOWER SUBPANELS SIMULTANEOUSLY - MEASURE DIFFERENCE BETWEEN MODELS?
dt <- "160430"
adj.u <- adj.l <- pw.m[,,"grey", dt] - pw.m[,,"black", dt]
adj.u[,1:1024] <- NA
adj.l[,1025:2048] <- NA

s.160430.u <- gaussian.spot.ls(adj.u,
                 c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                 x0.target = c(768, 1280), y0.target = c(768, 1280))
s.160430.l <- gaussian.spot.ls(adj.l,
                               c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                               x0.target = c(768, 1280), y0.target = c(768, 1280))

pixel.image(pw.m[,,"grey", dt] - pw.m[,,"black", dt])
contour(1:2048, 1:2048, gaussian.spot.mat(s.160430.u$par), add = T, drawlabels = F, lty = 1, lwd = 2)
contour(1:2048, 1:2048, gaussian.spot.mat(s.160430.l$par), add = T, drawlabels = F, lty = 2, lwd = 2)
