
library("IO.Pixels")

dates <- list.dirs("./Image-data/", full.names = F, recursive = F)
dates <- dates[dates != "150702"]
n <- length(dates)

#########################################################################################

# import all black images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
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

#########################################################################################

# import all grey images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
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

#########################################################################################

# import all white images, get pixelwise means & SDs (elapsed: 1110.51 for 11 dates)
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

#########################################################################################