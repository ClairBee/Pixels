
library("IO.Pixels")

dates <- list.dirs("./Image-data/", full.names = F, recursive = F)
dates <- dates[dates != "150702"]
n <- length(dates)

#########################################################################################

# import all black images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
{
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
}

#########################################################################################

# import all grey images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
{
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
}

#########################################################################################

# import all white images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
{
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
}

#########################################################################################

# numerical summaries of all images
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
    df[r, 3:6] <- c(mean(pw.m.b[,,i]), sd(pw.m.b[,,i]), mad(pw.m.b[,,i]), mean(pw.sd.b[,,i])
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

#########################################################################################

# need to write function to add single image batch to summary files