
library("IO.Pixels")

dt <- "loan"
px <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))

####################################################################################################

# IMPORT NEW IMAGES                                                                             ####

# load 20 black, 20 grey, 20 white images
# find median residuals, linear residuals, Gaussian spot residuals; save as objects

####################################################################################################

# IDENTIFY DEFECTS                                                                              ####

# identify abnormal pixels
bad.pixel.map <- function(dt) {
    
    # import prepared data for selected image
    pw.m <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
    md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7", acq.list = dt)
    linear.res <- load.objects("./02_Objects/linear-res/", otype = "l-res", acq.list = dt)
    spot.res <- load.objects("./02_Objects/spot-res/", otype = "s-res", acq.list = dt)
    
    bad.px <- list()
    
    # identify all bad pixels in all images
    if(sum(pw.m[,,"black"] == 65535, na.rm = T) > 0) {
        bad.px$hot <- data.frame(which(pw.m[,,"black"] == 65535, arr.ind = T), type = "hot")
    }
    
    if(sum(pw.m[,,"white"] <= 15000, na.rm = T) > 0) {
        bad.px$dark <- data.frame(which(pw.m[,,"white"] <= 15000, arr.ind = T), type = "dark")
    }
    
        bad.px$b <- classify.px(pw.m[,,"black"])
        bad.px$g <- classify.px(pw.m[,,"grey"])
    
        bad.px$b.res <- classify.px(md7[,,"black"], levels = c("l.dim", NA, "l.bright"))
        bad.px$g.res <- classify.px(md7[,,"grey"], levels = c("l.dim", NA, "l.bright"))
    
        bad.px$linear <- classify.px(linear.res, levels = c("nl.dim", NA, "nl.bright"))
        bad.px$spot <- classify.px(spot.res, levels = c("spot.dim", NA, "spot.bright"))
    
    # combine all bad pixel lists into single list, retaining only most severe classification
        px <- rbind.fill(bad.px)
        px$type <- ordered(px$type, 
                           levels = c("hot", "dark", "v.bright", "bright", "s.bright", "l.bright", "dim",
                                      "l.dim", "nl.bright", "nl.dim", "spot.bright", "spot.dim"))
        
        px <- px[order(px$type),]
        px <- px[!duplicated(px[,1:2]),]
        
        return(px)
}

px <- bad.pixel.map("loan")

pixel.plot(px, col = px.cols()[px$type])
table(px$type)

####################################################################################################

# IDENTIFY FEATURES                                                                             ####

# screen spots
    pw.m <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
    sp <- screen.spots(pw.m)
    px <- label.screen.spots(px, sp) 

# lines
    px <- find.rows(px)
    px <- find.columns(px)

# dense regions
    px <- dense.regions(px)

# clusters
    px <- find.clusters(px)

# remaining pixels must be singletons
    px$f.type[is.na(px$f.type)] <- "singleton"
    
    px$f.type <- factor(px$f.type)
    
table(px$f.type, px$type, useNA = "ifany")
pixel.plot(px, col = c("green3", "green", "grey", "black", "blue", "red", "skyblue")[px$f.type])

saveRDS(px, paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))

####################################################################################################

#                                   ####    WORKING AREA    ####                                ####

dt <- "130613"
pw.m <- load.objects("./02_Objects/images/", otype = "pwm")

process.img <- function(dt) {
    
  
    px <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))
    
    # screen spots
    sp <- screen.spots(pw.m)
    px <- label.screen.spots(px, sp) 
    
    # dense regions
    px <- dense.regions(px) 
    dpx <- px                       # temporary storage
    
    # lines
    cpx <- find.columns(px)
    cpx <- find.rows(cpx)

    # clusters
    cpx <- find.clusters(cpx)
    
    # remaining pixels must be singletons
    px$f.type[is.na(px$f.type)] <- "singleton"
    
    px$f.type <- factor(px$f.type)
    
    saveRDS(px, paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))
}

process.img("130613")
