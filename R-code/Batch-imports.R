
library("IO.Pixels")

dates <- list.dirs("./Image-data/", full.names = F, recursive = F)
dates <- dates[dates != "150702"]
n <- length(dates)

############################################################################################
# IMPORT & ADD A SINGLE IMAGE BATCH TO ALL SUMMARIES                                    ####
library(abind)

new.date <- "160430"
im <- load.daily(new.date)

# pixelwise mean values
{
    # black
    {
        b <- readRDS("./Other-data/Pixelwise-means-black.rds")
        b <- abind(b, pixelwise.mean(im[,,,1]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(b)[[3]], new.date)))
        saveRDS(b, "./Other-data/Pixelwise-means-black.rds")
    }
    
    # grey
    {
        g <- readRDS("./Other-data/Pixelwise-means-grey.rds")
        g <- abind(g, pixelwise.mean(im[,,,2]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(g)[[3]], new.date)))
        saveRDS(g, "./Other-data/Pixelwise-means-grey.rds")
    }
    
    # white
    {
        w <- readRDS("./Other-data/Pixelwise-means-white.rds")
        w <- abind(w, pixelwise.mean(im[,,,3]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(w)[[3]], new.date)))
        saveRDS(w, "./Other-data/Pixelwise-means-white.rds")
    }
}

# pixelwise SD
{
    # black
    {
        b <- readRDS("./Other-data/Pixelwise-sds-black.rds")
        b <- abind(b, pixelwise.sd(im[,,,1]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(b)[[3]], new.date)))
        saveRDS(b, "./Other-data/Pixelwise-sds-black.rds")
    }
    
    # grey
    {
        g <- readRDS("./Other-data/Pixelwise-sds-grey.rds")
        g <- abind(g, pixelwise.sd(im[,,,2]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(g)[[3]], new.date)))
        saveRDS(g, "./Other-data/Pixelwise-sds-grey.rds")
    }
    
    # white
    {
        w <- readRDS("./Other-data/Pixelwise-sds-white.rds")
        w <- abind(w, pixelwise.sd(im[,,,3]), along = 3,
                   new.names = list(NULL, NULL, c(dimnames(w)[[3]], new.date)))
        saveRDS(w, "./Other-data/Pixelwise-sds-white.rds")
    }
}
############################################################################################

# ADD AN IMAGE TO MEDIAN-DIFFERENCE LIST                                                ####

md.black <- readRDS("./Other-data/Med-diffs-black.rds")

md.black$"160430" <- pw.m[,,"black", "160430"] - r2m(focal(m2r(pw.m[,,"black", "160430"]), matrix(rep(1, 9), ncol = 3), fun = median))

saveRDS(md.black, "./Other-data/Med-diffs-black.rds")


############################################################################################

# import all black images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)    ####

    pw.m <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    pw.sd <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
        load.images(dates[i], "black")
        pw.m[,,i] <- pixelwise.mean(eval(parse(text = paste0("b.",dates[i]))))
        pw.sd[,,i] <- pixelwise.sd(eval(parse(text = paste0("b.",dates[i]))))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    remove(i, pb)
    
    saveRDS(pw.m, file = "./Other-data/Pixelwise-means-black.rds")
    saveRDS(pw.sd, file = "./Other-data/Pixelwise-sds-black.rds")



############################################################################################

# import all grey images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)     ####

    pw.m <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    pw.sd <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
        load.images(dates[i], "grey")
        pw.m[,,i] <- pixelwise.mean(eval(parse(text = paste0("g.",dates[i]))))
        pw.sd[,,i] <- pixelwise.sd(eval(parse(text = paste0("g.",dates[i]))))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    remove(i, pb)
    
    saveRDS(pw.m, file = "./Other-data/Pixelwise-means-grey.rds")
    saveRDS(pw.sd, file = "./Other-data/Pixelwise-sds-grey.rds")


############################################################################################

# import all white images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)    ####

    pw.m <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    pw.sd <- array(dim = c(1996, 1996, n), dimnames = list(NULL, NULL, dates))
    
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
        load.images(dates[i], "white")
        pw.m[,,i] <- pixelwise.mean(eval(parse(text = paste0("w.",dates[i]))))
        pw.sd[,,i] <- pixelwise.sd(eval(parse(text = paste0("w.",dates[i]))))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    remove(i, pb)
    
    saveRDS(pw.m, file = "./Other-data/Pixelwise-means-white.rds")
    saveRDS(pw.sd, file = "./Other-data/Pixelwise-sds-white.rds")


############################################################################################
    
# get pixelwise MAD for all images                                                      ####
    
    cols <- c("black", "grey", "white")
    pw.mad <- array(dim = c(1996, 1996, 3, n), dimnames = list(NULL, NULL, cols, dates))

    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (colour in cols) {
        l <- substring(colour, 1, 1)
        for (i in 1:n) {
            load.images(dates[i], colour)
            pw.mad[,, colour, i] <- apply(eval(parse(text = paste0(l, ".", dates[i]))), c(1,2), mad)
            setTxtProgressBar(pb, i)
        }
    }

    close(pb)
    remove(i, pb, colour, cols)
    
    saveRDS(pw.mad, file = "./Other-data/Pixelwise-mads.rds")
    

# numerical summaries of all images                                                     ####
pw.m.b <- readRDS("./Other-data/Pixelwise-means-black.rds")
pw.m.g <- readRDS("./Other-data/Pixelwise-means-grey.rds")
pw.m.w <- readRDS("./Other-data/Pixelwise-means-white.rds")

pw.sd.b <- readRDS("./Other-data/Pixelwise-sds-black.rds")
pw.sd.g <- readRDS("./Other-data/Pixelwise-sds-grey.rds")
pw.sd.w <- readRDS("./Other-data/Pixelwise-sds-white.rds")


range.summary <- function (date, batch) {
    data <- load.images(date, batch)
    c(min = min(data), lq = quantile(data, 0.25), median = median(data), uq = quantile(data, 0.75), 
                                     max = max(data), iqr = IQR(data))
}


df <- data.frame(date = character(), batch = character(),
                 mean = numeric(), sd = numeric(), mad = numeric, pw.sd = numeric(),
                 min = double(), lq = double(), median = double(), uq = double(), max = double(), iqr = double(),
                 stringsAsFactors = F)

for (i in 1:dim(pw.m.b)[[3]]) {
    r <- nrow(df) + 1
    df[r, 1:2] <- c(dimnames(pw.m.b)[[3]][i], "black")
    df[r, 3:6] <- c(mean(pw.m.b[,,i]), sd(pw.m.b[,,i]), mad(pw.m.b[,,i]), mean(pw.sd.b[,,i]))
    df[r, 7:12] <- range.summary(dimnames(pw.m.b)[[3]][i], "black")
    
    r <- nrow(df) + 1
    df[r, 1:2] <- c(dimnames(pw.m.g)[[3]][i], "grey")
    df[r, 3:6] <- c(mean(pw.m.g[,,i]), sd(pw.m.g[,,i]), mad(pw.m.g[,,i]), mean(pw.sd.g[,,i]))
    df[r, 7:12] <- range.summary(dimnames(pw.m.b)[[3]][i], "grey")
    
    r <- nrow(df) + 1
    df[r, 1:2] <- c(dimnames(pw.m.w)[[3]][i], "white")
    df[r, 3:6] <- c(mean(pw.m.w[,,i]), sd(pw.m.w[,,i]), mad(pw.m.w[,,i]), mean(pw.sd.w[,,i]))
    df[r, 7:12] <- range.summary(dimnames(pw.m.b)[[3]][i], "white")
}
rownames(df) <- apply(df[,c(2,1)], 1, paste, collapse = " ")
df <- df[order(rownames(df)),]

 write.csv(df, "./Other-data/Image-summaries.csv", row.names = F)
 df <- read.csv("./Other-data/Image-summaries.csv", as.is = T)

############################################################################################

# GET FF-CORRECTED IMAGES FROM PIXELWISE MEANS                                          ####
 ff.corr <- array(dim = c(1996, 1996, dim(pw.m)[[4]]), dimnames = list(NULL, NULL, dimnames(pw.m)[[4]]))
 
 for (d in dimnames(pw.m)[[4]]) {
     ff.corr[,,d] <- flat.field.corrected(d)
 }
 saveRDS(ff.corr, file = "./Other-data/Flat-field-corrected.rds")
 
 # need to write function to add single image batch to summary files                     ####
 
############################################################################################
 
 # IMPORT NON-STANDARD DATA                                                             ####
 
############################################################################################
 
# MEDIAN DIFFERENCES                                                                   ####
 ff <- gsub(".rds", "", gsub("md7-", "", 
                             list.files("./02_Objects/med-diffs",
                                        pattern = "md7-[a-z, A-Z, 0-9]+\\.rds$", full.names = F)))
 to.add <- dimnames(pw.m)[[4]][!(dimnames(pw.m)[[4]] %in% ff)]
 
 md <- apply(pw.m[,,,to.add[1]], 3, 
             function (im) r2m(focal(m2r(im), matrix(rep(1, 49), ncol = 7), fun = median)))
 md <- array(md, dim = c(2048, 2048, 3))
 saveRDS(pw.m[,,,to.add[1]] - md, paste0("./02_Objects/med-diffs/md7-", to.add[1], ".rds"))
 
 ############################################################################################
 
 # FEATURE IDENTIFICATION                                                                ####
 
 # label all feature types
 {
     process.img <- function(dt) {
         
         px <- bad.pixel.map(dt)
         
         cat(dt, "starting px:", nrow(px), "-")
         
         # screen spots
         sp <- readRDS(paste0("./02_Objects/pixel-maps/screen-spots-", dt, ".rds"))
         px <- label.screen.spots(px, sp) 
         
         # dense regions
         px <- dense.regions(px) 
         
         # lines
         px <- find.columns(px)
         px <- find.rows(px)
         
         # clusters
         px <- find.clusters(px)
         
         # remaining pixels must be singletons
         px$f.type[is.na(px$f.type)] <- "singleton"
         
         px$f.type <- factor(px$f.type)
         
         saveRDS(px, paste0("./02_Objects/pixel-maps/pixel-map-", dt, ".rds"))
         cat(" final px:", nrow(px), "\n")
     }
     
     imgs <- gsub("\\.rds", "", gsub("pwm-", "", list.files("./02_Objects/images")))
     
     invisible(lapply(imgs[20:22], process.img))
     
     process.img("141009")
     
 }