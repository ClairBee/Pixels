
# import new data

image.path <- "./Image-data/loan3/"

load.raw.images <- function(image.path) {
    
    img.out <- array(dim = c(2048, 2048, 3), dimnames = list(NULL, NULL, c("black", "grey", "white")))
    
    for (cc in dimnames(img.out)[[3]]) {
        im.data <- xmlToList(xmlParse(list.files(paste0(image.path, cc, "/"), pattern = "\\.xml$", full.names = T)[1]))
        
        x.offset <- as.integer(im.data$CameraProperties$imageOffsetX)
        y.offset <- as.integer(im.data$CameraProperties$imageOffsetY)
        x.size <- as.integer(im.data$CameraProperties$imageSizeX)
        y.size <- as.integer(im.data$CameraProperties$imageSizeY)
        
        im.list <- list.files(paste0(image.path, cc, "/"), pattern = "\\.tif$", full.names = T)
        
        imp <- abind(lapply(im.list, readTIFF, as.is = T), along = 3)
        pwm <- apply(imp, 1:2, mean)
        img.out[x.offset + c(1:x.size), y.offset + c(1:y.size), cc] <- t(pwm[nrow(pwm):1, , drop=FALSE])
    }
    return(img.out)
}

zz <- load.raw.images("./Image-data/loan2/")

saveRDS(img.out, "./02_Objects/images/pwm-loan2.rds")

#==============================================================================

# median-smoothed residuals

md7 <- array(dim = dim(img.out), dimnames = dimnames(img.out))
md7[,,"black"] <- img.out[,,"black"] - r2m(focal(m2r(img.out[,,"black"]), matrix(rep(1, 49), ncol = 7), fun = median))
md7[,,"grey"] <- img.out[,,"grey"] - r2m(focal(m2r(img.out[,,"grey"]), matrix(rep(1, 49), ncol = 7), fun = median))
md7[,,"white"] <- img.out[,,"black"] - r2m(focal(m2r(img.out[,,"white"]), matrix(rep(1, 49), ncol = 7), fun = median))

saveRDS(md7, "./02_Objects/med-diffs/md7-loan2.rds")

#==============================================================================

# linear residuals

lr <- fit.w.lm(img.out)

saveRDS(lr, "./02_Objects/linear-res/l-res-loan2.rds")

#==============================================================================

# fitted Gaussian spot

gs <- gaussian.spot.ls(img.out[,,"grey"] - img.out[,,"black"],
                 c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                 x0.target = c(768, 1280), y0.target = c(768, 1280))

sm <- rbind(spot.models,
            "loan2" = c(gs$par, rmse = sqrt(gs$value / sum(!is.na(img.out))), rep(NA, 6)))

write.csv(sm, "./Other-data/Gaussian-spots.csv")
sm
