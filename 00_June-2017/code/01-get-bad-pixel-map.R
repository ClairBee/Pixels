
library("IO.Pixels"); library("beepr")

img.nm <- "loan3"

############################################################################################################
# IMPORT IMAGE SET & CALCULATE ALL RESIDUALS                                                            ####
#   Import new images                                                                                   ####

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

img <- load.raw.images(paste0("./Image-data/", img.nm, "/"))

saveRDS(img, paste0("./02_Objects/images/pwm-", img.nm, ".rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Median-smoothed residuals                                                                           ####

md7 <- array(dim = dim(img), dimnames = dimnames(img))
md7[,,"black"] <- img[,,"black"] - r2m(focal(m2r(img[,,"black"]), matrix(rep(1, 49), ncol = 7), fun = median))
md7[,,"grey"] <- img[,,"grey"] - r2m(focal(m2r(img[,,"grey"]), matrix(rep(1, 49), ncol = 7), fun = median))
md7[,,"white"] <- img[,,"black"] - r2m(focal(m2r(img[,,"white"]), matrix(rep(1, 49), ncol = 7), fun = median))

saveRDS(md7, paste0("./02_Objects/med-diffs/md7-", img.nm, ".rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Linear residuals                                                                                    ####

lr <- fit.gv.lm(img)

saveRDS(lr, paste0("./02_Objects/linear-res/l-res-", img.nm, ".rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Gaussian spot                                                                                       ####

gs <- gaussian.spot.ls(img[,,"grey"] - img[,,"black"],
                       c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500),
                       x0.target = c(768, 1280), y0.target = c(768, 1280))

sm <- rbind(read.csv("./Other-data/Gaussian-spots.csv", stringsAsFactors = F, row.names = 1),
            img.nm = c(gs$par, rmse = sqrt(gs$value / sum(!is.na(img))), rep(NA, 6)))
rownames(sm)[nrow(sm)] <- img.nm

write.csv(sm, "./Other-data/Gaussian-spots.csv")
sm


############################################################################################################
# BAD PIXEL MAP                                                                                         ####

# load all images
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7", acq.list = img.nm)
pw.m <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = img.nm)
lr <- load.objects("./02_Objects/linear-res/", otype = "l-res", acq.list = img.nm)
gs <- load.objects("./02_Objects/spot-res/", otype = "s-res", acq.list = img.nm)

sc <- shading.corrected(pw.m)

# follow procedure in paper to identify abnormal pixels
{
    hot <- data.frame(which(pw.m[,,"black"] == 65535, arr.ind = T), "type" = "hot")
    dark <- data.frame(which(pw.m[,,"white"] < 15000, arr.ind = T), "type" = "dark")
    
    # globally extreme pixels
    gex.g <- classify.px(pw.m[,,"grey"])
    gex.b <- classify.px(pw.m[,,"black"])
    
    # locally non-uniform pixels
    lnu.b <- classify.px(md7[,,"black"], levels = c("l.dim", NA, "l.bright"))
    lnu.g <- classify.px(md7[,,"grey"], levels = c("l.dim", NA, "l.bright"))
    
    # nonlinear response
    nl <- classify.px(lr, levels = c("nl.dim", NA, "nl.bright"))
    
    # nonstandard response to Gaussian spot
    ng <- classify.px(gs, levels = c("spot.dim", NA, "spot.bright"))
    
    # spots on screen
    ss <- screen.spots(pw.m)
    
    bad.pixel.map <- rbind(hot, dark, gex.g, gex.b, lnu.b, lnu.g, nl, ng)
    bad.pixel.map$type <- ordered(bad.pixel.map$type, 
                                  levels = c("hot", "dark", "v.bright", "bright", "s.bright", "l.bright", "dim", 
                                             "l.dim", "nl.bright", "nl.dim", "spot.bright", "spot.dim"))
    bad.pixel.map <- bad.pixel.map[order(bad.pixel.map$type),]
    
    bpm <- bad.pixel.map[!duplicated(bad.pixel.map[,1:2]),]
    
    # check against archived pixel map to confirm
    # table(bpm$type)
    # table(readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", img.nm, ".rds"))$type)
    
    sp <- readRDS(paste0("./02_Objects/pixel-maps/screen-spots-", "loan", ".rds"))
}

saveRDS(bpm, paste0("./02_Objects/pixel-maps/pixel-map-", img.nm, ".rds"))
saveRDS(ss, paste0("./02_Objects/pixel-maps/screen-spots-", img.nm, ".rds"))

############################################################################################################
# IDENTIFY FEATURES                                                                                     ####

f.cols <- c("cl.body" = "gold", "cl.root" = "red", "dense.region" = "blue", "line.c" = "black", "line.r" = "black",
            "singleton" = "skyblue", "s.spot" = "grey")

fl <- list.files("./02_Objects/pixel-maps", pattern = "pixel-map")
fl <- gsub(".rds", "", gsub("pixel-map-", "", fl))[-(1:2)]

invisible(sapply(fl, function(img.nm) {
    
    bpm <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", img.nm, ".rds"))
    org.ftype <- bpm[,"f.type"]
    bpm <- bpm[,1:3]
    
    # add feature classifications
    sp <- readRDS(paste0("./02_Objects/pixel-maps/screen-spots-", img.nm, ".rds"))
    
    bpm <- label.screen.spots(bpm, sp) 
    
    bpm <- dense.regions(bpm)
    
    bpm <- find.columns(bpm)
    bpm <- find.rows(bpm)
    
    bpm <- find.clusters(bpm)
    bpm$f.type[is.na(bpm$f.type)] <- "singleton"
    
    cat(img.nm, "\n", "Original categories", "\n")
    print(table(org.ftype, useNA = "ifany"))
    
    cat("\n", "New categories")
    print(table(bpm$f.type, useNA = "ifany"))
    
    pixel.plot(bpm[,1:2], col = f.cols[bpm$f.type], main = img.nm)
    legend("top", pch = 15, col = f.cols[unique(bpm$f.type)], unique(bpm$f.type), ncol = 3, bty = "n")
    
    beep()
    yn <- readline("Save this map?")
    if(yn != "n") saveRDS(bpm, paste0("./02_Objects/pixel-maps/pixel-map-", img.nm, ".rds"))
}))
