
fpath <- "./Notes/Universal-thresholding/fig/"

df <- readRDS(paste0(fpath, "all-px.rds"))

####################################################################################################

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")
.sd <- hijack(sd, na.rm = T)

####################################################################################################

# CLASSIFICATION FUNCTIONS                                                                      ####

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

assign.official <- function(dat) {
    
    # clear any previous categories
    dat$t.off <- NA
    
    # nonuniformity
    dat$t.off[dat$sc > 1.02 * median(dat$sc, na.rm = T)] <- "g.nonuniform"
    dat$t.off[dat$sc < 0.98 * median(dat$sc, na.rm = T)] <- "g.nonuniform"
    
    # noise
    dat$t.off[dat$sd.w > 6 * median(dat$sd.w, na.rm = T)] <- "bright.noise"
    dat$t.off[dat$sd.b > 6 * median(dat$sd.b, na.rm = T)] <- "dark.noise"
    
    # sensitivity (based on gain image)
    dat$t.off[dat$w - dat$b > 1.5 * median(dat$w - dat$b, na.rm = T)] <- "up.bright"
    dat$t.off[dat$w - dat$b < .45 * median(dat$w - dat$b, na.rm = T)] <- "up.dark"
    dat$t.off[dat$w - dat$b < 5000 & dat$w < 10000] <- "no.gain"
    
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

med.replace <- function(im, px, w = 5) {
    
    get.px <- cbind(px[1] + rep(c(-floor(w/2):floor(w/2)), w),
                    px[2] + sort(rep(c(-floor(w/2):floor(w/2)), w)))
    get.px <- get.px[get.px[,1] %in% c(1:2048) & get.px[,2] %in% c(1:2048),]
    median(im[get.px], na.rm = T)
}

pixel.plots <- function(cn, ...) {
    om <- par()$mar
    par(mfrow = c(2,2), mar = c(0,0,1,1))
        lapply(names(df), function(nm) {
            px <- df[[nm]]
            px <- px[!is.na(px[, cn]),]
            pixel.plot(px, main = paste0(nm, " - ", cn), ...)
        })
    par(mfrow = c(1,1), mar = om)
}

####################################################################################################

# CLASSIFY BAD PIXELS                                                                           ####

# official classification
df <- lapply(df, assign.official)
{
    pixel.plots("t.off", cex = 0.5)
    lapply(df, function(dd) table(dd$t.off))
}

# non-linear pixels (in grey image)
df <- lapply(df, assign.glm)
{
    pixel.plots("t.glm", cex = 0.5, col = "grey")
    lapply(df, function(dd) table(dd$t.glm, dd$t.off, useNA = "ifany"))
}

# new thresholds
df <- lapply(df, assign.category)
{
    pixel.plots("type", cex = 0.5, col = "darkblue")
    lapply(df, function(dd) table(dd$t.glm, dd$type, useNA = "ifany"))
    lapply(df, function(dd) table(dd$t.off, dd$type, useNA = "ifany"))
}

# median-switch values where defects identified
{
    pw.m <- abind(sapply(c("131122", "140128", "MCT225", "160430"), 
                         function(nm) readRDS(paste0("./02_Objects/images/pwm-", nm, ".rds")), 
                         simplify = F),
                  along = 4)
    
    fixed <- abind(sapply(names(df), 
                          function(nm) {
                              px <- df[[nm]]
                              px <- as.matrix(px[!is.na(px$type),1:2])
                              ms <- pw.m[,,,nm]
                              ms[,,"black"][px] <- apply(px, 1, med.replace, im = ms[,,"black"])
                              ms[,,"grey"][px] <- apply(px, 1, med.replace, im = ms[,,"grey"])
                              ms[,,"white"][px] <- apply(px, 1, med.replace, im = ms[,,"white"])
                              return(ms)
                          }, simplify = F), along = 4)

    # histograms of original vs median-switched values
    {
        hist(pw.m[,,"black", "MCT225"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"black", "MCT225"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"grey", "MCT225"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"grey", "MCT225"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"white", "MCT225"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"white", "MCT225"], breaks = "fd", add = T, col = "black")
        
        #----------------------------------------------------------------------------------
        
        hist(pw.m[,,"black", "160430"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"black", "160430"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"grey", "160430"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"grey", "160430"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"white", "160430"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"white", "160430"], breaks = "fd", add = T, col = "black")
        
        #----------------------------------------------------------------------------------
        
        hist(pw.m[,,"black", "140128"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"black", "140128"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"grey", "140128"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"grey", "140128"], breaks = "fd", add = T, col = "black")
        
        hist(pw.m[,,"white", "140128"], breaks = "fd", ylim = c(0,30), border = "skyblue")
        hist(fixed[,,"white", "140128"], breaks = "fd", add = T, col = "black")
    }
}

# LINE CLASSIFICATION NEEDS TO BE SORTED OUT
{
    plot.line.diffs <- function(ll, im) {
        dd <- ddply(data.frame(which(ll > 0, arr.ind = T)), .(x = row), summarise,
                    ymin = min(col), ymax = max(col), length = length(row))
        print(dd)
        apply(dd, 1, 
              function(px) {
                  plot(im[px["x"],] - im[px["x"]-1,],
                       type = "l", main = paste0("Column ", px["x"]), xlab = "", ylab = "",
                       xlim = ((px["ymin"] > 1024) + c(0,1)) * 1024)
                  points(im[px["x"],] - im[px["x"]-1,], pch = ".", cex = 2,
                         col = c(NA, "red")[c(1:2048) %in% c(px["ymin"]:px["ymax"]) + 1])
                  abline(h = 0, col = "darkred")
                  abline(h = c(300) * c(-1,1), col = "cyan3", lty = 2)
                  abline(v = 1024.5, col = "blue")
              })
    }

    # classify lines in median-switched data

    # 131122
    {
        ll.131122 <- find.lines(fixed[,,"grey", "131122"], dim.lines = F, threshold = 4000) + 
                     find.lines(fixed[,,"grey", "131122"], dim.lines = T, threshold = 4000)
        plot.line.diffs(ll.131122, pw.m[,,"grey", "131122"])
        
        ll.131122 <- find.lines(md[,,"black", "131122"], dim.lines = F, threshold = 4000) + 
                     find.lines(md[,,"black", "131122"], dim.lines = T, threshold = 4000)
        plot.line.diffs(ll.131122, pw.m[,,"black", "131122"])
        
        plot(pw.m[1396,,"black", "131122"], type = "l", xlim = c(1024, 2048))
        lines(pw.m[1397,,"black", "131122"], col = "blue")
        lines(pw.m[1398,,"black", "131122"], col = "red")
        
        plot(md[1396, , "black", "131122"], type = "l", xlim = c(1024, 2048))
        lines(md[1397,,"black", "131122"], col = "blue")
        lines(md[1398,,"black", "131122"], col = "red")
        
    }
    # larger smoothing kernel may help?
    md5 <- pw.m[,,"black", "131122"] - r2m(focal(m2r(pw.m[,,"black", "131122"]), matrix(rep(1, 25), ncol = 5), fun = median))
    md7 <- pw.m[,,"black", "131122"] - r2m(focal(m2r(pw.m[,,"black", "131122"]), matrix(rep(1, 49), ncol = 7), fun = median))

    ll7 <- find.lines(md7, threshold = 4000) + find.lines(md7, dim.lines = T, threshold = 4000)
    
    plot.line.diffs(ll7, pw.m[,,"black", "131122"])
    
    plot(pw.m[522,,"black", "131122"], xlim = c(0,1024), type = "l")
    points(pw.m[522,,"black", "131122"], pch = ".", col = c(NA, "red")[c(1:1024) %in% c(735:912) + 1])
}

####################################################################################################

# EXTRACT & COMPARE BAD PIXEL MAPS                                                              ####

# extract each type of bad pixel map
bpx.off <- lapply(df, function(px) px[!is.na(px$t.off),])       # 8431, 686, 3281, 6280
bpx.cb <- lapply(df, function(px) px[!is.na(px$type),])         # 9664, 97941, 3351, 7852 
bpx.nl <- lapply(df, function(px) px[!is.na(px$t.glm),])        # 9004, 686, 2586, 3105

# high shading correction median diff (300px brighter than neighbours)  # 9872, 680, 2296, 6681 
bpx.scmd <- lapply(df, function(px) px[which(abs(px$sc.md) > 200),])

pixel.plot(bpx.cb[[i]], cex = 0.4, col = "magenta3")
points(bpx.off[[i]][,1:2], pch = 15, cex = 0.4, col = "blue")
points(bpx.nl[[i]][,1:2], pch = 15, cex = 0.4, col = "green3")
points(bpx.scmd[[i]][,1:2], pch = 15, cex = 0.4, col = "cyan3")

# 131122: CB has more general scatter, NL pixks up some scatter & more points in corners
# 140128: CB has lots of bright pixels (because of wide spread of values), others v similar
# MCT225: broadly similar situation to 131122
# 160430: all v. similar, general noisiness in detector forming circular pattern


    
####################################################################################################

# QUADRAT TESTS                                                                                 ####

apply(which(!is.na(fixed[,,"black", "140128"]), arr.ind = T), 2, range)

im.params <- list("131122" = list(x.dim = c(25, 2024), y.dim = c(225, 1824),
                                  x.panels = panel.edges()$x, y.panels = panel.edges()$y),
                  "140128" = list(x.dim = c(1, 2000), y.dim =  c(29, 2028),
                                  x.panels = panel.edges()$x, y.panels = panel.edges()$y),
                  "MCT225" = list(x.dim = c(25, 2024), y.dim =  c(25, 2024),
                                  x.panels = panel.edges()$x, y.panels = range(panel.edges()$y)),
                  "160430" = list(x.dim = c(3,1998), y.dim = c(33,2028),
                                  x.panels = panel.edges()$x, y.panels = panel.edges()$y))

# official pixels
sapply(names(bpx.off), function(nm) {
    px <- bpx.off[[nm]]
    im.dims <- im.params[[nm]]
    quadrat.test(ppp(px$x, px$y, im.dims$x.dim, im.dims$y.dim),
                 xbreaks = im.dims$x.panels - 0.5, ybreaks = im.dims$y.panels - 0.5)
}, simplify = F)

# CB pixels
sapply(names(bpx.off), function(nm) {
    px <- bpx.cb[[nm]]
    im.dims <- im.params[[nm]]
    quadrat.test(ppp(px$x, px$y, im.dims$x.dim, im.dims$y.dim),
                 xbreaks = im.dims$x.panels - 0.5, ybreaks = im.dims$y.panels - 0.5)
}, simplify = F)


####################################################################################################