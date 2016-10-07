
library("IO.Pixels")

dt <- "130613"
px <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))

# rework column detection:
#   - check for double/triple columns
#   - pad image with 0 at edges: otherwise edge px are lost

# also rework dense column detection:
#   - pad with repeated rows/columns to allow classification at edges (eg. 141009 column 3)

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

px <- bad.pixel.map(dt)

pixel.plot(px, col = px.cols()[px$type])
table(px$type)

####################################################################################################

# IDENTIFY FEATURES                                                                             ####

# screen spots
pw.m <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
sp <- screen.spots(pw.m)
px <- label.screen.spots(px, sp) 

# dense regions
px <- dense.regions(px)

# lines
px <- find.rows(px)
px <- find.columns(px)

# clusters
px <- find.clusters(px)

# remaining pixels must be singletons
px$f.type[is.na(px$f.type)] <- "singleton"

px$f.type <- factor(px$f.type)

table(px$f.type, px$type, useNA = "ifany")
pixel.plot(px, col = c("green3", "green", "grey", "black", "blue", "red", "skyblue")[px$f.type])

saveRDS(px, paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))

####################################################################################################

# FEATURES BY SIZE                                                                              ####

dt <- "130613"
px <- readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))

count.columns(px, midline = 1024.5)


# clusters by size & dominant direction
count.clusters(px)
count.clusters(px, check.dir = T)


####################################################################################################

# CLUSTERING BEHAVIOUR                                                                          ####

# convert pixels to point process object

px.pp <- as.ppp

####################################################################################################

# SPATIAL DISTRIBUTION OF FEATURES                                                              ####

# convert pixel map to ppp and crop to avoid dense regions

# quadrat tests of subpanels

# quadrat tests of vertical / horizontal spread



####################################################################################################

#                                   ####    WORKING AREA    ####                                ####

# how much to crop images by, to avoid influence of dense regions?
imgs <- gsub("\\.rds", "", gsub("pwm-", "", list.files("./02_Objects/images")))

all.px <- rbind.fill(lapply(imgs, function(dt) readRDS(paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))))

pixel.plot(all.px[all.px$f.type == "dense.region",])
rect(256, 256, 1792, 1792, border = "red")
rect(100, 100, 1948, 1948, border = "purple")

# crop by 256


blah <- px2ppp(all.px[!(all.px$f.type %in% c("s.spot", "dense.region", "line.c", "line.r")),], crop = list(c(256, 1792), c(256, 1792)))

plot(blah, pch = 15, cex = 0.4, asp = T)
draw.panels(col = "skyblue")
