
# is 140128 properly aligned? Bright lines seem to cross midpoint.
library("IO.Pixels"); library("CB.Misc")

pw.m <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                     simplify = F),
              along = 4)
md <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                   function(nm) readRDS(paste0("./02_Objects/med-diffs/md-", nm, ".rds")), 
                   simplify = F),
            along = 4)
pw.sd <- abind(sapply(c("MCT225", "160430"), 
                     function(nm) readRDS(paste0("./02_Objects/stdevs/pwsd-", nm, ".rds")), 
                     simplify = F),
              along = 4)
pw.sd <- abind(pw.sd, 
               "140128" = array(NA, dim = c(2048, 2048, 3), dimnames = dimnames(pw.sd[,,,1])),
               "131122" = array(NA, dim = c(2048, 2048, 3), dimnames = dimnames(pw.sd[,,,1])), along = 4)

sc <- array(apply(pw.m, 4, shading.corrected), dim = dim(pw.m[,,1,]), dimnames = dimnames(pw.m[,,1,]))

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")
.sd <- hijack(sd, na.rm = T)

fpath <- "./Notes/Universal-thresholding/fig/"

# haven't picked up bright lines separately b/c convolution currently erratic

# at each step, consider white vs grey/black relationship
# shading correction histogram

####################################################################################################

# SETUP                                                                                         ####

df <- sapply(dimnames(pw.m)[[4]], 
             function (nm) setNames(data.frame(melt(pw.m[,,"black", nm]),
                                 melt(pw.m[,,"grey", nm]),
                                 melt(pw.m[,,"white", nm]),
                                 melt(md[,,"black", nm]),
                                 melt(md[,,"grey", nm]),
                                 melt(md[,,"white", nm]),
                                 melt(sc[,, nm]),
                                 melt(pw.sd[,,"black", nm]),
                                 melt(pw.sd[,,"white", nm]),
                                 melt(sc.md[,, nm])
                                 )[, c("X1", "X2", "value", "value.1", "value.2", "value.3", "value.4", "value.5", "value.6", "value.7", "value.8", "value.9")], 
                      nm = c("x", "y", "b", "g", "w", "md.b", "md.g", "md.w", "sc", "sd.b", "sd.w", "sc.md")),
             simplify = F)

saveRDS(df, )

assign.category <- function(dat, dark.th = 5000) {
    
    # support functions
    th.u <- function(vals) {
        med <- median(vals, na.rm = T)
        med + (65535 - med)/4
    }
    th.l <- function(vals) {
        med <- median(vals, na.rm = T)
        med * 0.75
    }
    local.th <- 2 * c("b" = .sd(dat$b), "g" = .sd(dat$g), "w" = .sd(dat$w))
    
    # clear any previous categories
    dat$type <- NA
    
    # reassign categories
    dat$type[(dat$md.b < -local.th["b"]) | (dat$md.g < -local.th["g"]) | (dat$md.w < -local.th["w"])] <- "local.l"
    dat$type[(dat$md.b > local.th["b"]) | (dat$md.g > local.th["g"]) | (dat$md.w > local.th["w"])] <- "local.u"
    dat$type[(dat$b < th.l(dat$b) | dat$g < th.l(dat$g) | dat$w < th.l(dat$w))] <- "global.l"
    dat$type[dat$g - dat$b < dark.th] <- "dark"
    dat$type[(dat$b > th.u(dat$b) | dat$g > th.u(dat$g) | dat$w > th.u(dat$w))] <- "global.u"
    
    dat$type <- ordered(dat$type, levels = c("dark", "line", "global.u", "global.l", "local.u", "local.l"))
    dat
}

df <- lapply(df, assign.category)  

lapply(df, function(dd) table(dd$type))


c.cols <- c("black", "green3", "red", "blue", "gold", "cyan3")
                 
# plot bad pixels       
lapply(names(df), 
       function(nm) {
           px <- df[[nm]]
           pdf(paste0(fpath, "bad-px-plot-", nm, ".pdf"))
           par(mar = c(2,2,1,1))
           plot(px[,1:2][!is.na(px$type),], pch = ".", cex = 2, col = c.cols[px$type[!is.na(px$type)]],
                xlab = "", ylab = "")
           dev.off()
       })

####################################################################################################

# HISTOGRAMS                                                                                    ####

th.hist <- function(im, crop = NA) {
    par(mar = c(2,2,1,1))
    if (is.na(crop)) {
        hh <- hist(im, breaks = "fd", main = "", xlab = "", ylab = "", xlim = c(0,65535))
    } else {
        hh <- hist(im, breaks = "fd", main = "", xlab = "", ylab = "", ylim = c(0, crop), xlim = c(0,65535))
    }
    
    ym <- max(pretty(hh$counts)) * 1.5
    
    med <- median(im, na.rm = T)
    rect(med + (65535 - med)/2, 0, 65535, ym, col = adjustcolor("red", alpha = 0.3), border = NA)
    rect(med + (65535 - med)/4, 0, med + (65535 - med)/2, ym, col = adjustcolor("orange", alpha = 0.3), border = NA)
    rect(0, 0, med /2, ym, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
    rect(med /2, 0, med * 0.75, ym, col = adjustcolor("cyan3", alpha = 0.3), border = NA)
}

th.hist(pw.m[,,"black", "140128"], crop = 30)
th.hist(pw.m[,,"grey", "140128"], crop = 30)
th.hist(pw.m[,,"white", "140128"], crop = NA)
th.hist(pw.m[,,"white", "MCT225"], crop = NA)


# current 'globally extreme' setting picks up a lot of px from 14-01-28. Part of same 'doughnut' shape
apply(pw.m[,,"white",], 3, sd, na.rm = T)

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

####################################################################################################

# MEDIAN-SWITCH DEFECTS                                                                         ####

med.replace <- function(im, px, w = 5) {
    
    get.px <- cbind(px[1] + c(-floor(w/2):floor(w/2)), px[2])
    get.px <- get.px[get.px[,1] %in% c(1:2048),]
    median(im[get.px], na.rm = T)
}

mr <- abind(sapply(names(df), 
             function(nm) {
                 px <- df[[nm]]
                 px <- as.matrix(px[!is.na(px$type),1:2])
                 ms <- pw.m[,,,nm]
                 ms[,,"black"][px] <- apply(px, 1, med.replace, im = ms[,,"black"])
                 ms[,,"grey"][px] <- apply(px, 1, med.replace, im = ms[,,"grey"])
                 ms[,,"white"][px] <- apply(px, 1, med.replace, im = ms[,,"white"])
                 return(ms)
                 }, simplify = F), along = 4)

hist(mr[,,"black", "160430"], breaks = "fd")

plot(which(mr[,,"white", "160430"] < 35000, arr.ind = T), pch = 15, xlim = c(0,2048), ylim = c(0,2048))
points(which(mr[,,"grey", "160430"] < 15000, arr.ind = T), pch = 20, cex = 0.8, col = "cyan3")
points(which(mr[,,"grey", "160430"] > 20000, arr.ind = T), pch = 20, cex = 0.8, col = "red")
points(which(mr[,,"black", "160430"] > 10000, arr.ind = T), pch = 0, col = "gold")



####################################################################################################

# MANUALLY ASSIGN BRIGHT/DIM LINES (FOR NOW)                                                    ####

# 131122
{
    
}

# 140128
{
    ll <- abind("bright" = find.lines(pw.m[,,"black", "140128"]),
                "dim" = find.lines(pw.m[,,"black", "140128"], dim.lines = T), along = 3)
    ll.px <- which(ll[,,"bright"] > 0 | ll[,,"dim"] > 0, arr.ind = T)
           
    # check defective columns
    ddply(data.frame(ll.px), .("x" = row), summarise,
          ymin = min(col), ymax = max(col), length = length(row))
    
    o.plot(pw.m[745,,"black", "140128"] - pw.m[744,,"black", "140128"], xlim = c(0,1024))
    abline(v = c(536.5, 1024.5), col = "red", lty = 2)
    o.plot(pw.m[747,,"black", "140128"] - pw.m[746,,"black", "140128"], xlim = c(0,1024))
    # 745 bright, 746 is an artefact
    
    o.plot(pw.m[1043,,"black", "140128"] - pw.m[1042,,"black", "140128"], xlim = c(1025, 2048))
    
    o.plot(pw.m[1043,,"black", "140128"] - pw.m[1042,,"black", "140128"], xlim = c(1025, 2048))
    abline(v = c(1133.5, 1024.5), col = "red", lty = 2)
    # 1043 bright
    
          
    image(1:2048, 1:2048, ll, col = c(NA, "red","blue", "green3", "purple"))
}

# 160430
{
    df$"160430"$type[find.lines(pw.m[,,"black", "160430"]) > 0] <- "line"
    plot(df$"160430"[!is.na(df$"160430"$type),1:2], pch = ".", cex = 2, 
         col = c.cols[df$"160430"$type[!is.na(df$"160430"$type)]])
}

# MCT225
{
    
}

apply(which(!is.na(pw.m[,,"black", "140128"]), arr.ind = T), 2, range)

####################################################################################################

# SHADING-CORRECTED IMAGES                                                                      ####

# shading corrected median differences
sc.md <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                      function(nm) readRDS(paste0("./02_Objects/sc-med-diffs/sc-md-", nm, ".rds")), 
                      simplify = F),
               along = 3)

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

# histogram of median differences of shading-corrected values, excluding identified defects
lapply(names(df),
       function(nm) {
           sc.im <- sc.md[,,nm]
           sc.im[as.matrix(df[[nm]][!is.na(df[[nm]]$type), c("x", "y")])] <- NA
           pdf(paste0(fpath, "sc-mid-hist-excl-defects-", nm, ".pdf"), height = 4)
           par(mar = c(2,2,1,1))
           hist(sc.im, breaks = "fd", ylim = c(0,30), xlab = "", ylab = "", main = "")
           dev.off()
       })

####################################################################################################

# PLOT REMAINING HIGH-RESIDUAL PIXELS                                                           ####

lapply(df, function(px) {
    plot(px[is.na(px$type) & abs(px$res) > 700, 1:2], pch = 15, xlim = c(0,2048), ylim = c(0,2048))
})

####################################################################################################

# 'OFFICIAL' THRESHOLDS                                                                         ####

# sensitivity
tmp <- df$"160430"
tmp$t2[tmp$w - tmp$b > 1.5 * median(tmp$w - tmp$b, na.rm = T)] <- "up.bright"
tmp$t2[tmp$w - tmp$b < .45 * median(tmp$w - tmp$b, na.rm = T)] <- "up.dark"
tmp$t2[tmp$w - tmp$b < 5000 & tmp$w < 10000] <- "no.gain"

table(tmp$t2, tmp$type, useNA = "ifany")

# noise
tmp$sd.w <- c(pw.sd[,,"white", "160430"])
tmp$sd.b <- c(pw.sd[,,"black", "160430"])
tmp$t2[tmp$sd.w > 6 * median(tmp$sd.w, na.rm = T)] <- "bright.noise"
tmp$t2[tmp$sd.b > 6 * median(tmp$sd.b, na.rm = T)] <- "dark.noise"

table(tmp$t2, tmp$type, useNA = "ifany")

# uniformity
tmp$sc <- c(sc[,,"160430"])
tmp$sc.md <- c(sc.md[,,"160430"])

tmp$t2[tmp$sc > 1.02 * median(tmp$sc, na.rm = T)] <- "g.nonuniform"
tmp$t2[tmp$sc < 0.98 * median(tmp$sc, na.rm = T)] <- "g.nonuniform"

tmp$t2[tmp$sc > 1.01 * tmp$sc.md] <- "l.nonuniform"
tmp$t2[tmp$sc < 0.98 * tmp$sc.md] <- "l.nonuniform"

table(tmp$t2, tmp$type, useNA = "ifany")
# local nonuniformity is particularly bonkers...

assign.official <- function(dat) {

        # clear any previous categories
        dat$t.off <- NA
        
        # nonuniformity
        dat$t.off[dat$sc > 1.02 * median(dat$sc, na.rm = T)] <- "g.nonuniform"
        dat$t.off[dat$sc < 0.98 * median(dat$sc, na.rm = T)] <- "g.nonuniform"
        
        # noise
        dat$t.off[dat$sd.w > 6 * median(dat$sd.w, na.rm = T)] <- "bright.noise"
        dat$t.off[dat$sd.b > 6 * median(dat$sd.b, na.rm = T)] <- "dark.noise"
        
        # sensitivity
        dat$t.off[dat$w - tmp$b > 1.5 * median(dat$w - dat$b, na.rm = T)] <- "up.bright"
        dat$t.off[dat$w - tmp$b < .45 * median(dat$w - dat$b, na.rm = T)] <- "up.dark"
        dat$t.off[dat$w - tmp$b < 5000 & dat$w < 10000] <- "no.gain"
        
        dat
}

df <- lapply(df, assign.official)  

assign.sc.nonuniform

####################################################################################################

# NONLINEARITY                                                                                  ####

assign.wlm <- function(dat, fit.region = c(40.5, 2008.5), trunc = T, res.th = 500) {
    
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
    dat$t.wlm[abs(dat$res) > res.th] <- "nonlinear"
    dat
}
assign.glm <- function(dat, fit.region = c(40.5, 2008.5), trunc = T, res.th = 200) {
    
    w.lm <- lm(g ~ b * w, 
               data = dat[findInterval(dat$x, fit.region) == 1 & 
                              findInterval(dat$y, fit.region) == 1, ])
    dat$fv.g <- predict(w.lm, dat[, c("b", "w")])
    
    # optionally, truncate fitted values at 65535
    if (trunc) {
        dat$fv.g[dat$fv.g > 65535] <- 65535
    }
    dat$res.g <- dat$g - dat$fv.g
    attr(dat, "fit.g") = list(r2 = round(summary(w.lm)$adj.r.squared, 3), 
                            rmse = round(summary(w.lm)$sigma, 2))
    dat$t.glm[abs(dat$res.g) > res.th] <- "nonlinear"
    dat
}

df <- lapply(df, assign.wlm)
df <- lapply(df, assign.glm)

lapply(df, function(dd) table(dd$t.wlm, dd$type, useNA = "ifany"))
lapply(df, function(dd) table(dd$t.glm, dd$type, useNA = "ifany"))

lapply(names(df), function(nm) {
    dd <- df[[nm]]
    .smoothScatter(dd[,c("res.g", "res")], xlab = "Grey residual", ylab = "white residual", main = nm)
    abline(line(dd[,c("res.g", "res")]), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
}) # grey vs white residuals (all PX)

lapply(names(df), function(nm) {
    dd <- df[[nm]]
    dd <- dd[is.na(dd$type),]
    .smoothScatter(dd[,c("res.g", "res")], xlab = "Grey residual", ylab = "white residual", 
                   main = paste0(nm, " - CB px removed"))
    abline(line(dd[,c("res.g", "res")]), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
}) # grey vs white residuals (CB healthy px)

lapply(names(df), function(nm) {
    dd <- df[[nm]]
    dd <- dd[is.na(dd$t.off),]
    .smoothScatter(dd[,c("res.g", "res")], xlab = "Grey residual", ylab = "white residual", 
                   main = paste0(nm, " - official bad px removed"))
    abline(line(dd[,c("res.g", "res")]), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
}) # grey vs white residuals (official healthy px)

####################################################################################################

# SHADING-CORRECTED MEDIAN DIFFERENCES                                                          ####

assign.sc.nonuniform <- function(dat, sc.md.th = 300) {
    
    # clear any previous categories
    dat$t.sc.md <- NA
    
    # nonuniformity
    dat$t.sc.md[dat$sc.md > sc.md.th] <- "sc.l.bright"
    dat$t.sc.md[dat$sc.md < -sc.md.th] <- "sc.l.dim"
    
    dat
}

df <- lapply(df, assign.sc.nonuniform)

lapply(df, function(dd) table(dd$t.glm, dd$type, useNA = "ifany"))

lapply(df, function(dd) s.hist(dd$sc))


####################################################################################################

# ASSESS RESULTS                                                                                ####
 
# histogram of all shading-corrected values
lapply(names(df),
       function(nm) {
           dd <- df[[nm]]
           par(mar = c(2,2,1,1))
           hh <- hist(dd$sc, breaks = "fd", ylim = c(0,30), col = "black",
                      xlab = "", ylab = "", main = paste0(nm, " shading corrections"))
       })

# histogram of shading-corrected values excluding identified defects
lapply(names(df),
       function(nm) {
           dd <- df[[nm]]
           par(mar = c(2,2,1,1))
           hh <- hist(dd$sc[is.na(dd$type)], breaks = "fd", ylim = c(0,30), col = "black",
                      xlab = "", ylab = "", main = paste0(nm, " excl. CB bpx"))
       })
{
    tmp <- df$"160430"
    tmp <- tmp[which(is.na(tmp$type) & (tmp$sc > 18250 | tmp$sc < 17250)),]
    
}

# histogram of shading-corrected values excluding official defects
lapply(names(df),
       function(nm) {
           dd <- df[[nm]]
           par(mar = c(2,2,1,1))
           hh <- hist(dd$sc[is.na(dd$t.off)], breaks = "fd", ylim = c(0,30), col = "black",
                      xlab = "", ylab = "", main = paste0(nm, " excl. official bpx"))
       })

# histogram of shading-corrected values excluding official defects
lapply(names(df),
       function(nm) {
           dd <- df[[nm]]
           par(mar = c(2,2,1,1))
           hh <- hist(dd$sc[is.na(dd$t.sc.md)], breaks = "fd", ylim = c(0,30), col = "black",
                      xlab = "", ylab = "", main = paste0(nm, " excl. SC-md nonuniform bpx"))
       })

lapply(df, function(dd) table(dd$t.sc.md, dd$type, useNA = "ifany"))


# non-uniformity doesn't exclude all px with high or low SD
# histogram of shading-corrected values excluding official defects
lapply(names(df),
       function(nm) {
           dd <- df[[nm]]
           par(mar = c(2,2,1,1))
           hh <- hist(dd$sc[is.na(dd$t.off)], breaks = "fd", ylim = c(0,30), col = "black",
                      xlab = "", ylab = "", main = paste0(nm, " excl. official bpx"))
       })

.smoothScatter(df$"160430"[,c("res.g", "sc.md")], xlab = "Grey residual", ylab = "SC median diff")
.smoothScatter(df$"160430"[,c("res.g", "sc.md")], nrpoints = 100,
               col = px.cols[df$"160430"$type], pch = 20,
               xlim = c(-5000, 5000), ylim = c(-10000, 10000), xlab = "Grey residual", ylab = "SC median diff")
px.cols <- c(NA, "darkblue", "green3", "magenta3", "skyblue", "pink", "cyan3")

points(df$"160430"[!is.na(df$"160430"$type), c("res.g", "sc.md")], pch = ".", cex = 2,
       col = adjustcolor(px.cols[df$"160430"$type[!is.na(df$"160430"$type)]], alpha = 0.3))

dd <- df$"160430"
hist(dd$sc.md, breaks = "fd", ylim = c(0,30), col = "darkblue", border = "darkblue")
hist(dd$sc.md[dd$type != "dark"], breaks = "fd", add = T, col = "magenta3", border = "magenta3")
hist(dd$sc.md[!dd$type %in% c("dark", "global.u")], breaks = "fd", add = T, col = "green3", border = "green3")
hist(dd$sc.md[!dd$type %in% c("dark", "global.u", "line")], breaks = "fd", add = T, col = "skyblue", border = "skyblue")

points(df$"160430"[df$"160430"$type == "dark", c("res.g", "sc.md")], 
       col = adjustcolor("darkblue", alpha = 0.4), pch = 20)
points(df$"160430"[df$"160430"$type == "global.u", c("res.g", "sc.md")], 
       col = adjustcolor("magenta3", alpha = 0.4), pch = 20)
points(df$"160430"[df$"160430"$type == "global.l", c("res.g", "sc.md")], 
       col = adjustcolor("skyblue", alpha = 0.4), pch = 20)
points(df$"160430"[df$"160430"$type == "local.u", c("res.g", "sc.md")], 
       col = adjustcolor("pink", alpha = 0.4), pch = 20)
points(df$"160430"[df$"160430"$type == "local.l", c("res.g", "sc.md")], 
       col = adjustcolor("pink", alpha = 0.4), pch = 20)


####################################################################################################

# LOCAL VARIATION FOR THRESHOLDING                                                              ####

local.sd <- r2m(focal(m2r(pw.m[,,"black", "160430"]), matrix(c(rep(1, 84), 0, rep(1, 84)), ncol = 13), fun = sd))
mean(local.sd, na.rm = T)

local.mad <- r2m(focal(m2r(pw.m[,,"black", "160430"]), matrix(rep(1, 169), ncol = 13), fun = mad))
mean(local.mad, na.rm = T)
