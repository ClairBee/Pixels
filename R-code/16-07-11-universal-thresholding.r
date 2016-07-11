

library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)
md <- abind(sapply(c("131122", "MCT225", "160430"), 
                   function(nm) readRDS(paste0("./02_Objects/med-diffs/md-", nm, ".rds")), 
                   simplify = F),
            along = 4)

sc <- array(apply(pw.m, 4, shading.corrected), dim = dim(pw.m[,,1,]), dimnames = dimnames(pw.m[,,1,]))

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

fpath <- "./Notes/Universal-thresholding/fig/"

# haven't picked up bright lines separately b/c convolution currently erratic

# at each step, consider white vs grey/black relationship
# shading correction histogram

####################################################################################################

# SETUP                                                                                         ####

df <- sapply(c("131122", "MCT225", "160430"), 
             function (nm) setNames(data.frame(melt(pw.m[,,"black", nm]),
                                 melt(pw.m[,,"grey", nm]),
                                 melt(pw.m[,,"white", nm]),
                                 melt(md[,,"black", nm]),
                                 melt(md[,,"grey", nm]),
                                 melt(md[,,"white", nm]))[, c("X1", "X2", "value", "value.1", "value.2", "value.3", "value.4", "value.5")], 
                      nm = c("x", "y", "b", "g", "w", "md.b", "md.g", "md.w")),
             simplify = F)

assign.category <- function(dat, local.th = 1000, dark.th = 5000) {
    
    # support functions
    th.u <- function(vals) {
        med <- median(vals, na.rm = T)
        med + (65535 - med)/4
    }
    th.l <- function(vals) {
        med <- median(vals, na.rm = T)
        med * 0.75
    }
    
    # clear any previous categories
    dat$type <- NA
    
    # reassign categories
    dat$type[(dat$md.b < -local.th) | (dat$md.g < -local.th) | (dat$md.w < -local.th)] <- "local.l"
    dat$type[(dat$md.b > local.th) | (dat$md.g > local.th) | (dat$md.w > local.th)] <- "local.u"
    dat$type[(dat$b < th.l(dat$b) | dat$g < th.l(dat$g) | dat$w < th.l(dat$w))] <- "global.l"
    dat$type[(dat$b > th.u(dat$b) | dat$g > th.u(dat$g) | dat$w > th.u(dat$w))] <- "global.u"
    dat$type[dat$g - dat$b < dark.th] <- "dark"
    
    dat$type <- ordered(dat$type, levels = c("dark", "global.u", "global.l", "local.u", "local.l"))
    dat
}

df <- lapply(df, assign.category, local.th = 1500)  

table(df$"131122"$type)
table(df$"MCT225"$type)
table(df$"160430"$type)

c.cols <- c("black", "red", "blue", "gold", "cyan3")
                        
lapply(df, 
       function(px) {
           plot(px[,1:2][!is.na(px$type),], pch = ".", col = c.cols[px$type[!is.na(px$type)]],
                xlab = "", ylab = "")
       })

####################################################################################################

# REGRESSION OF WHITE VALUES                                                                    ####

wlm <- function(dat, fit.region = c(40.5, 2008.5), trunc = T) {
    
    w.lm <- lm(w ~ b * g, 
               data = dat[findInterval(dat$x, fit.region) == 1 & 
                              findInterval(dat$y, fit.region) == 1, ])
    dat$fv <- predict(w.lm, dat[, c("b", "g")])
    
    # optionally, truncate fitted values at 65535
    if (trunc) {
        dat$fv[dat$fv > 65535] <- 65535
    }
    dat$res <- dat$w - dat$fv
    attr(dat, "fit") = list(r2 = round(summary(w.lm)$adj.r.squared, 3), 
                            rmse = round(summary(w.lm)$sigma, 2))
    dat
}

df <- lapply(df, wlm)

lapply(names(df), 
       function(nm) {
           px <- df[[nm]]
           pdf(paste0(fpath, "wlm-all-", nm, ".pdf"))
           par(mar = c(4,4,1,1))
           .smoothScatter(px$fv, px$w, ylim = c(0,65535), xlim = c(0,65535), 
                          xlab = "Fitted", ylab = "Observed")
           abline(0,1, lty = 2, col = adjustcolor("darkred", alpha = 0.4))
           dev.off()
       })

####################################################################################################

# REMOVE DARK PIXELS                                                                            ####

# exclude dark pixels
lapply(names(df), 
       function(nm) {
           px <- df[[nm]]
           px <- px[px$type != "dark",]
           pdf(paste0(fpath, "wlm-ex-dark-", nm, ".pdf"))
           par(mar = c(4,4,1,1))
           .smoothScatter(px$fv, px$w, xlim = c(0, 65535), ylim = c(0,65535), 
                          xlab = "Fitted", ylab = "Observed")
           abline(0,1, lty = 2, col = adjustcolor("darkred", alpha = 0.4))
#           points(df[[nm]][df[[nm]]$type == "dark", c("fv", "w")], pch = 20, 
#                  col = adjustcolor("cyan3", alpha = 0.1))
           dev.off()
       })

# exclude dark  & globally extreme pixels
lapply(names(df), 
       function(nm) {
           px <- df[[nm]]
           px <- px[!px$type %in% c("dark", "global.u", "global.l"),]
           pdf(paste0(fpath, "wlm-ex-global-", nm, ".pdf"))
           par(mar = c(4,4,1,1))
           .smoothScatter(px$fv, px$w, xlim = c(0, 65535), ylim = c(0,65535), 
                          xlab = "Fitted", ylab = "Observed")
           abline(0,1, lty = 2, col = adjustcolor("darkred", alpha = 0.4))
           #           points(df[[nm]][df[[nm]]$type == "dark", c("fv", "w")], pch = 20, 
           #                  col = adjustcolor("cyan3", alpha = 0.1))
           dev.off()
       })

# only unclassified pixels
lapply(names(df), 
       function(nm) {
           px <- df[[nm]]
           px <- px[is.na(px$type),]
           pdf(paste0(fpath, "wlm-unclassified-", nm, ".pdf"))
           par(mar = c(4,4,1,1))
           .smoothScatter(px$fv, px$w, xlim = c(0, 65535), ylim = c(0,65535), 
                          xlab = "Fitted", ylab = "Observed")
           abline(0,1, lty = 2, col = adjustcolor("darkred", alpha = 0.4))
           #           points(df[[nm]][df[[nm]]$type == "dark", c("fv", "w")], pch = 20, 
           #                  col = adjustcolor("cyan3", alpha = 0.1))
           dev.off()
       })

# shading-corrected images excluding identified defects
lapply(names(df),
       function(nm) {
           sc.im <- sc[,,nm]
           sc.im[as.matrix(df[[nm]][!is.na(df[[nm]]$type), c("x", "y")])] <- NA
           jpeg(paste0(fpath, "sc-image-excl-defects-", nm, ".jpg"))
           par(mar = c(2,2,1,1))
           pixel.image(sc.im)
           dev.off()
       })

# histogram of shading-corrected values excluding identified defects
lapply(names(df),
       function(nm) {
           sc.im <- sc[,,nm]
           sc.im[as.matrix(df[[nm]][!is.na(df[[nm]]$type), c("x", "y")])] <- NA
           pdf(paste0(fpath, "sc-hist-excl-defects-", nm, ".pdf"), height = 4)
           par(mar = c(2,2,1,1))
           hist(sc.im, breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = "")
           dev.off()
       })

# check shading-corrected median differences as well?

####################################################################################################

# PLOT REMAINING HIGH-RESIDUAL PIXELS                                                           ####

