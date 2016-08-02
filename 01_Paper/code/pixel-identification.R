
library("IO.Pixels"); library("CB.Misc")
fpath <- "./01_Paper/fig/exploratory/"

pw.m <- load.pixel.means()
pw.sd <- load.pixel.sds()

sc <- array(apply(pw.m, 4, shading.corrected), dim = dim(pw.m)[-3], dimnames = dimnames(pw.m)[-3])

# need to get remaining pixelwise SDs (for old data)
# check SDs of dead pixels

####################################################################################################

# METHOD FROM MANUAL                                                                            ####

# get offset-corrected images
os <- pw.m[,,c("grey", "white"), ] - pw.m[,,c("black", "black"), ] # offset-corrected images

# get median-differenced shading corrected images
{
    zz <- apply(sc, 3, function(im) r2m(focal(m2r(im), matrix(rep(1, 81), nrow = 9), fun = median)))
    saveRDS(scmd.99, "./02_Objects/sc-med-diffs/sc-md-9x9.rds")
}
scmd.99 <- readRDS("./02_Objects/sc-med-diffs/sc-md-9x9.rds")


# signal sensitivity
{
    ss <- list(up.b.g = which(os[,,"grey", "160430"] > 1.5 * median(os[,,"grey", "160430"], na.rm = T), arr.ind = T),
               up.b.w = which(os[,,"white", "160430"] > 1.5 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T),
               n.g = which(os[,,"grey", "160430"] < 10000 | os[,,"white", "160430"] < 15000, arr.ind = T),
               up.d.g = which(os[,,"grey", "160430"] < 0.45 * median(os[,,"grey", "160430"], na.rm = T), arr.ind = T),
               up.d.w = which(os[,,"white", "160430"] < 0.45 * median(os[,,"white", "160430"], na.rm = T), arr.ind = T))
               
    spx <- rbind.fill(lapply(names(ss[sapply(ss, nrow) > 0]), function(tp) data.frame(ss[[tp]], type = tp)))
}

# noise
{
    np <- list(bn.g = which(pw.sd[,,"grey", "160430"] > 6 * median(pw.sd[,,"grey", "160430"], na.rm = T), arr.ind = T),
               bn.w = which(pw.sd[,,"white", "160430"] > 6 * median(pw.sd[,,"white", "160430"], na.rm = T), arr.ind = T),
               dn.b = which(pw.sd[,,"black", "160430"] > 6 * median(pw.sd[,,"black", "160430"], na.rm = T), arr.ind = T))
    npx <- rbind.fill(lapply(names(np[sapply(np, nrow) > 0]), function(tp) data.frame(np[[tp]], type = tp)))
}

# uniformity
{
    uu <- list(gl.u = which(sc[,,"160430"] > 1.02 * median(sc[,,"160430"], na.rm = T) | sc[,,"160430"] < 0.98 * median(sc[,,"160430"], na.rm = T), arr.ind = T),
               loc.u = which(sc[,,"160430"] > 1.01 * scmd.99[,,"160430"] | sc[,,"160430"] < 0.99 * scmd.99[,,"160430"], arr.ind = T))
    upx <- rbind.fill(lapply(names(uu[sapply(uu, nrow) > 0]), function(tp) data.frame(uu[[tp]], type = tp)))
}

# combine into single map
{
    # merge & remove duplicates
    px <- rbind(spx, npx, upx)
    px$type <- ordered(px$type, levels = c("n.g", "up.b.g", "up.d.g", "up.d.w", "dn.b", "gl.u", "loc.u"))
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    
    pixel.plot(px, cex = 0.5, col = c("black", "red", "blue", "blue", "green3", "orange", NA)[px$type])
}


# and now, to apply over all images
{
    bpx <- sapply(dimnames(pw.m)[[4]],
                 function (dt) {
                     
                     # only run noise check if pw sd exists 
                     if (dt %in% dimnames(pw.sd)[[4]]) {
                         nn <- list(bright.noise = which(pw.sd[,,"grey", dt] > 6 * median(pw.sd[,,"grey", dt], na.rm = T) | 
                                                             pw.sd[,,"white", dt] > 6 * median(pw.sd[,,"white", dt], na.rm = T), arr.ind = T),
                                    dark.noise = which(pw.sd[,,"black", dt] > 6 * median(pw.sd[,,"black", dt], na.rm = T), arr.ind = T))
                         nn <- rbind.fill(lapply(names(nn[sapply(nn, nrow) > 0]), function(tp) data.frame(nn[[tp]], type = tp)))
                     } else {
                         nn <- NULL
                     }
                     
                     ss <- list(no.gain = which(os[,,"grey", dt] < 10000 | os[,,"white", dt] < 15000, arr.ind = T),
                                up.bright = which(os[,,"grey", dt] > 1.5 * median(os[,,"grey", dt], na.rm = T) |
                                                   os[,,"white", dt] > 1.5 * median(os[,,"white", dt], na.rm = T), arr.ind = T),
                                up.dim = which(os[,,"grey", dt] < 0.45 * median(os[,,"grey", dt], na.rm = T) |
                                                   os[,,"white", dt] < 0.45 * median(os[,,"white", dt], na.rm = T), arr.ind = T),
                                gl.unif = which(sc[,,dt] > 1.02 * median(sc[,,dt], na.rm = T) | sc[,,dt] < 0.98 * median(sc[,,dt], na.rm = T), arr.ind = T),
                                local.unif = which(sc[,,dt] > 1.01 * scmd.99[,,dt] | sc[,,dt] < 0.99 * scmd.99[,,dt], arr.ind = T))
                    
                     ss <- rbind.fill(lapply(names(ss[sapply(ss, nrow) > 0]), function(tp) data.frame(ss[[tp]], type = tp)))
                     px <- rbind(nn, ss)
                     px$type <- ordered(px$type, levels = c("no.gain", "up.bright", "up.dim", "dark.noise", "bright.noise", "gl.unif", "local.unif"))
                     px <- px[order(px$type),]
                     px <- px[!duplicated(px[,1:2]),]
                     }, simplify = F)
}

saveRDS(bpx, paste0(fpath, "official-bad-px.rds"))

px.summ <- rbind.fill.matrix(sapply(bpx, function(px) t(as.matrix(table(px$type))), simplify = F))
row.names(px.summ) <- dimnames(pw.m)[[4]]

bpx <- readRDS(paste0(fpath, "official-bad-px.rds"))

table(bpx$"141009"$type)

ssp <- screen.spots(pw.m[,,"white", "141009"], enlarge = T)

####################################################################################################

# MEDIAN-DIFFERENCING OVER INDIVIDUAL IMAGES                                                    ####

####################################################################################################

# MEDIAN-DIFFERENCING OVER SHADING CORRECTION                                                   ####

####################################################################################################

# NONLINEARITY                                                                                  ####

####################################################################################################

# COMPARE RESULTS                                                                               ####

# check for constant offset in any bad pixels
{
    
}

# check correlation between pixelwise SD and median difference (per image/power)
{
    
}

