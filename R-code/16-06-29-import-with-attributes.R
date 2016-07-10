
# change image import so that each acquisition is stored as 3-layer array with dimensions

import.acq <- function(acq.folder, subfolders = c("black", "grey", "white"), panel.dim = c(2048, 2048)) {
    
    acq <- array(dim = c(panel.dim, length(subfolders)),
                 dimnames = list(NULL, NULL, subfolders))
    
    for (sf in subfolders) {
        
        # load xml data into data frame (get image list & offsets)
        im.list <- list.files(paste0(acq.folder, "/", sf), pattern = "\\.xml$", full.names = T)
        df <- rbind.fill(lapply(im.list, 
                                function(im.xml) as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(im.xml)), list(recursive = TRUE))), 
                                                                      FUN = unlist), stringsAsFactors = F)))
        # gives ordered list of images - may be needed if investigating as time series
        df <- data.frame(acq.time = df$Time,
                         filenm = df$.attrs.FileName,
                         offset.x = df$CameraProperties.imageOffsetX,
                         offset.y = df$CameraProperties.imageOffsetY,
                         size.x = df$CameraProperties.imageSizeX,
                         size.y = df$CameraProperties.imageSizeY,
                         kV = df$kV,
                         uA = df$uA,
                         exp.time = df$ExposureTime,
                         stringsAsFactors = F)[order(acq.time = df$Time),]
        
        # quick error check - are all images on same area of detector (should always be T)
        if (!identical(sapply(df[3:6], min), sapply(df[3:6], max))) {
            cat("Error - TIF images do not all have same offset & dimension")
            return(df)
        }
        
        # if all images cover same area of detector, continue
        offset.x <- as.integer(df$offset.x[1])
        offset.y <- as.integer(df$offset.y[1])
        size.x <- as.integer(df$size.x[1])
        size.y <- as.integer(df$size.y[1])
        
        ims <- apply(cbind(acq.folder, "/", sf, "/", df$filenm), 1, paste, collapse = "")
        
        # get pixelwise mean of image
        tmp <- apply(abind(lapply(ims, readTIFF, as.is = T), along = 3), 1:2, mean)
        tmp <- t(tmp[nrow(tmp):1, , drop = FALSE])

        # position correctly in array
        acq[(offset.x + 1) : (offset.x+size.x),
            (panel.dim[2]-size.y-offset.y + 1) : (panel.dim[2]-offset.y),
            sf] <- tmp
    }
    return(acq)
}

im.160430 <- import.acq("/home/clair/Documents/Pixels/Image-data/160430")


saveRDS(im.160430, "./02_Objects/images/acq-160430.rds")

im.MCT225 <-  import.acq("/home/clair/Documents/Pixels/Image-data/MCT225")
saveRDS(im.MCT225, "./02_Objects/images/acq-MCT225.rds")

im.160705 <-  import.acq("/home/clair/Documents/Pixels/Image-data/160705")
saveRDS(im.160705, "./02_Objects/images/acq-160705.rds")

####################################################################################################

# CONVERT EXISTING IMAGES TO SAME FORMAT                                                        ####

# time to load as single image: 10.576
system.time(load.pixel.means())

dt <- "160314"
acq <- array(dim = c(2048, 2048, 3), dimnames = list(NULL, NULL, c("black", "grey", "white")))
acq[3: 1998, 33 : 2028,] <- pw.m[,,,dt]
saveRDS(acq, paste0("./02_Objects/images/pwm-", dt, ".rds"))

zz <- lapply(list.files("./02_Objects/images", pattern = "pwm-[0-9]+\\.rds$", full.names = T), readRDS)
zz <- abind(zz, along = 4)

load.pixel.means.2 <- function(fpath = "./02_Objects/images") {
    ll <- list.files("./02_Objects/images", pattern = "pwm-[0-9]+\\.rds$", full.names = T)

    pw.m <<- abind(lapply(ll, readRDS), along = 4)
    dimnames(pw.m)[[4]] <<- unlist(lapply(ll, substring, 25, 30))
}

system.time(load.pixel.means.2())

####################################################################################################

# ALSO WITH MEDIAN-DIFFERENCED IMAGES                                                           ####

md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

md <- array(dim = c(2048, 2048, 2), dimnames = list(NULL, NULL, c("black", "grey")))
md[3: 1998, 33 : 2028,"black"] <- md.b$"160430"
md[3: 1998, 33 : 2028,"grey"] <- md.g$"160430"

saveRDS(md, "./02_Objects/med-diffs/md-160430.rds")

lapply(names(md.b), function(dt) {
    md[3: 1998, 33 : 2028,"black"] <- md.b[[dt]]
    md[3: 1998, 33 : 2028,"grey"] <- md.g[[dt]]
    saveRDS(md, paste0("./02_Objects/med-diffs/md-", dt, ".rds"))
})