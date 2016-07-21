
# update pixel.plot function to handle ppp objects
# next step: check quadratic trend model (K-function etc) for each pixel map
# look at deviations from quadratic trend
# correlation between nonlinearity & brightness is how strong?

library("IO.Pixels"); library("CB.Misc")
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
        med + (65535 - med)/2
    }
    th.l <- function(vals) {
        med <- median(vals, na.rm = T)
        med * 0.5
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

# MANUALLY ADD EXTRA LINES OF BAD PIXELS                                                        ####

# 131122
{
    pixel.plot(bpx.scmd$"131122", xlim = c(1800, 1900), ylim = c(1000, 1100))
    click()
    
    # find lines maually (all of these have had neighbour columns tested as well)
    # 17 dark columns to identify - all single columns
    o.plot(pw.m[88,,"grey", "131122"] - pw.m[87,,"grey", "131122"])     # 88 (upper)
    o.plot(pw.m[159,,"grey", "131122"] - pw.m[158,,"grey", "131122"]) # 159 (lower)
    o.plot(pw.m[522,,"grey", "131122"] - pw.m[521,,"grey", "131122"]) # 522 (lower)
    o.plot(pw.m[610,,"grey", "131122"] - pw.m[609,,"grey", "131122"]) # 610 (upper)
    o.plot(pw.m[736,,"grey", "131122"] - pw.m[735,,"grey", "131122"]) # 736 (lower)
    o.plot(pw.m[959,,"grey", "131122"] - pw.m[958,,"grey", "131122"]) # 959 (lower)
    o.plot(pw.m[1080,,"grey", "131122"] - pw.m[1079,,"grey", "131122"]) # 1080 (lower)
    o.plot(pw.m[1302,,"grey", "131122"] - pw.m[1301,,"grey", "131122"]) # 1302 (upper)
    o.plot(pw.m[1314,,"grey", "131122"] - pw.m[1313,,"grey", "131122"]) # 1312 (lower)
    o.plot(pw.m[1327,,"grey", "131122"] - pw.m[1326,,"grey", "131122"]) # 1327 (upper)
    o.plot(pw.m[1398,,"grey", "131122"] - pw.m[1397,,"grey", "131122"]) # 1398 (upper)
    o.plot(pw.m[1396,,"grey", "131122"] - pw.m[1395,,"grey", "131122"]) # 1396 (upper)
    o.plot(pw.m[1526,,"grey", "131122"] - pw.m[1525,,"grey", "131122"]) # 1526 (upper)
    o.plot(pw.m[1704,,"grey", "131122"] - pw.m[1703,,"grey", "131122"]) # 1704 (lower)
    o.plot(pw.m[1766,,"grey", "131122"] - pw.m[1767,,"grey", "131122"]) # 1766 (lower)
    o.plot(pw.m[1893,,"grey", "131122"] - pw.m[1892,,"grey", "131122"]) # 1893 (upper)
    o.plot(pw.m[1998,,"grey", "131122"] - pw.m[1997,,"grey", "131122"]) # 1998 (upper)
    
    dc <- list("u" = c(88, 1302, 1327, 1398, 1396, 1526, 1998, 1893),
               "l" = c(159, 522, 610, 736, 959, 1080, 1312, 1704, 1766))
    
    sp.131122 <- array(readRDS("./02_Objects/med-diffs/md7-131122.rds"),
                       dim = c(2048, 1024, 2, 3), 
                       dimnames = list(NULL, NULL, c("l", "u"), c("black", "grey", "white")))
    pixel.image(sp.131122[,, "u", "grey"])
    
    pdf("./dark-line-736-131122.pdf")
    plot(pw.m[159,,"white", "131122"], type = "l", col = "gold", xlim = c(0,1024), ylim = c(4800,6000))
    lines(pw.m[159,,"grey", "131122"], col = "green3")
    lines(pw.m[159,,"black", "131122"])
    legend("topleft", col = c("gold", "green3", "black"), lty = 1, legend = c("white", "grey", "black"), bty = "n")
    dev.off()
    
    # now use change-point detection (ecp.divisive/bcp) to identify line ends
    
    
    # 
    ddply(bpx.scmd$"131122", .(x, l.id), summarise, ymin = min(y), ymax = max(y),
          range = max(y) - min(y) + 1, length = length(x), coverage = length / range)
    bpx.scmd$"131122"[bpx.scmd$"131122"$x == 1997 & is.na(bpx.scmd$"131122"$l.id),]
    
    pixel.image(pw.m[,,"grey", "131122"], xlim = c(80, 120), ylim = c(1130, 1140))
    points(bpx.scmd$"131122"[bpx.scmd$"131122"$l.id == "line.64",], pch = 0)
}

# MCT225
{
    
}


####################################################################################################

# EXTRACT & COMPARE BAD PIXEL MAPS                                                              ####

# extract each type of bad pixel map
bpx.off <- lapply(df, function(px) px[!is.na(px$t.off),])       # 8431, 686, 3281, 6280
bpx.cb <- lapply(df, function(px) px[!is.na(px$type),])         # 9664, 97941, 3351, 7852 
bpx.nl <- lapply(df, function(px) px[!is.na(px$t.glm),])        # 9004, 686, 2586, 3105

# high shading-corrected median diff (300px brighter than neighbours)  # 9872, 680, 2296, 6681 
bpx.scmd <- lapply(df, function(px) px[which(abs(px$sc.md) > 200),])

pixel.plot(bpx.cb[[i]], cex = 0.4, col = "magenta3")
points(bpx.off[[i]][,1:2], pch = 15, cex = 0.4, col = "blue")
points(bpx.nl[[i]][,1:2], pch = 15, cex = 0.4, col = "green3")
points(bpx.scmd[[i]][,1:2], pch = 15, cex = 0.4, col = "cyan3")

# 131122: CB has more general scatter, NL pixks up some scatter & more points in corners
# 140128: CB has lots of bright pixels (because of wide spread of values), others v similar
# MCT225: broadly similar situation to 131122
# 160430: all v. similar, general noisiness in detector forming circular pattern

# identify column features so that they can be easily removed
tag.lines <- function(px) {
    
    # identify & group all adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px[,1:2], type = 1), im.dim = c(2048, 2048))), dir = 4)
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))), 
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    # filter out & retain only long line segments (> 5px)
    ll <- ddply(xy, .(col = x, id), summarise, 
                ymin = min(y), ymax = max(y), length = max(y)-min(y) + 1)
    ll <- ll[ll$length > 5,]
    
    # check if any lines actually occur
    if (nrow(ll) > 0) {
        # turn line summary into list of coordinates
        zz <- do.call("rbind", apply(ll, 1, 
                                     function(l) cbind(x = l["col"], 
                                                       y = c(l["ymin"]:l["ymax"]), 
                                                       l.id = paste0("line.", l["id"]))))
        
        # tag identified pixels in original list of coordinates & return
        px <- merge(px, data.frame(zz), by = c(1:2), all.x = T)
    } else {
        px$"l.id" <- NA
    }

    return(px)
}

bpx.off <- lapply(bpx.off, tag.lines)
bpx.cb <- lapply(bpx.cb, tag.lines)
bpx.nl <- lapply(bpx.nl, tag.lines)
bpx.scmd <- lapply(bpx.scmd, tag.lines)

# identify clusters to that they can easily be removed
tag.clusters <- function(px) {
    
    # identify & group all adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px[,1:2], type = 1), im.dim = c(2048, 2048))), dir = 4)
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))), 
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    # filter out & retain only clusters containing multiple clusters
    # could use more sophisticated way to tag centres
    ll <- ddply(xy, .(id), summarise, 
                xmean = floor(mean(x)), ymean = floor(mean(y)), size = length(x))
    
    # check if any clusters actually occur, tag in pixel list 
    if (nrow(ll) > 0) {
        xy <- merge(xy, ll[ll$size > 1,], by = "id", all.x = T)
        
        xy$cl.type[!is.na(xy$size)] <- "cl.body"
        xy$cl.type[which(xy$x == xy$xmean & xy$y == xy$ymean)] <- "cl.root"
 
    } else {
        xt$cl.type <- NA
    }
    
    px <- merge(px, xy[,c("x", "y", "cl.type", "cl.id" = "id")], by = c(1:2), all = T)
    
    px$cl.type[!is.na(px$l.id)] <- "line.body"
    px$cl.type[is.na(px$cl.type)] <- "singleton"
    return(px)
}

bpx.off <- lapply(bpx.off, tag.clusters)
bpx.cb <- lapply(bpx.cb, tag.clusters)
bpx.nl <- lapply(bpx.nl, tag.clusters)
bpx.scmd <- lapply(bpx.scmd, tag.clusters)

# convert to list of ppp objects for easier reference
im.dims <- lapply(df, function(px) list(x = range(px$x[!is.na(px$b)]), y = range(px$y[!is.na(px$b)])))

bpx2ppp <- function(bpx, excl.feat = c("cl.body", "line.body")) {
    sapply(names(bpx), function(nm) {
        px <- bpx[[nm]]
        
        # filter lines if necessary
        px <- px[!(px$cl.type %in% excl.feat),]

        pp <- ppp(px$x, px$y, im.dims[[nm]]$x, im.dims[[nm]]$y)
    }, simplify = F)
}

ppp.off <- bpx2ppp(bpx.off, excl.feat = c("cl.body", "line.body"))
ppp.cb <- bpx2ppp(bpx.cb, excl.feat = c("cl.body", "line.body"))
ppp.nl <- bpx2ppp(bpx.nl, excl.feat = c("cl.body", "line.body"))
ppp.scmd <- bpx2ppp(bpx.scmd, excl.feat = c("cl.body", "line.body"))

                                    #           all pixels                      roots only
                                    # 131122 140128 MCT225 160430       131122 140128 MCT225 160430
sapply(ppp.off, "[[", "n")          #     97     81     48   5892           75     57     32   5694
sapply(ppp.cb, "[[", "n")           #    529    244    118   6862          443    175     75   6557
sapply(ppp.nl, "[[", "n")           #   1701     81    374   2654         1327     56    305   2545
sapply(ppp.scmd, "[[", "n")         #    153     77    546   6681          103     55    276   6574


####################################################################################################

# DESCRIBE CLUSTER DISTRIBUTION (SIZE & BRIGHTNESS)                                             ####
px <- bpx.cb$"160430"
px <- px[!(px$cl.type %in% c("line.body")), ]
hh <- ddply(px, .(id), summarise,
            xm = floor(mean(x)), ym = floor(mean(y)), size = length(x),
            bmean = mean(b), gmean = mean(g), wmean = mean(w),
            bsum = sum(b), gsum = sum(g), wsum = sum(w))
table(hh$size)

plot(hh[,c("xm", "ym")], 
     pch = c(20, 16, 16)[findInterval(hh$size, c(0, 1.5, 4.5, 999))],
     col = c("gold", "blue", "red")[findInterval(hh$size, c(0, 1.5, 4.5, 999))])

####################################################################################################

# QUADRAT TESTS                                                                                 ####

# need to rework this code - MCT225 panels are wrong

# official pixels
sapply(names(bpx.off), function(nm) {
    quadrat.test(ppp.off[[nm]],
                 xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)
}, simplify = F)

# CB pixels
sapply(names(bpx.off), function(nm) {
    quadrat.test(ppp.cb[[nm]],
                 xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)
}, simplify = F)

# nonlinear pixels
sapply(names(bpx.off), function(nm) {
    quadrat.test(ppp.nl[[nm]],
                 xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)
}, simplify = F)

# shading-corrected median-differenced pixels
sapply(names(bpx.off), function(nm) {
    quadrat.test(ppp.scmd[[nm]],
                 xbreaks = panel.edges()$x - 0.5, ybreaks = panel.edges()$y - 0.5)
}, simplify = F)

####################################################################################################

# ENVELOPE FUNCTIONS                                                                            ####

plot(envelope(ppp.cb$"160430", Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T))
plot(envelope(ppp.cb$"160430", Fest, nsim = 99, nrank = 2, verbose = F, fix.n = T))
plot(envelope(ppp.cb$"160430", Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T, transform = expression(. - pi * r ** 2)))

plot(envelope(ppp.cb$"131122", Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T, transform = expression(. - pi * r ** 2)))
plot(envelope(ppp.cb$"140128", Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T, transform = expression(. - pi * r ** 2)))
plot(envelope(ppp.cb$"MCT225", Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T, transform = expression(. - pi * r ** 2)))

pixel.plot(bpx.cb$"MCT225")
pixel.image(pw.m[,,"grey", "MCT225"], xlim = c(1500, 1600), ylim = c(1150, 1200))

pixel.plot(bpx.cb$"MCT225"[bpx.cb$"MCT225"$cl.type == "cl.root", 1:2], xlim = c(1500,1600), ylim = c(1200,1400))
points(bpx.cb$"MCT225"[bpx.cb$"MCT225"$cl.type == "singleton", 1:2], pch = 15, col = "red")
points(ppp.cb$"MCT225"$x, ppp.cb$"MCT225"$y, col = "red", pch = 20)
pixel.plot(ppp.cb$"MCT225"$x, ppp.cb$"MCT225"$y))

####################################################################################################

# NONPARAMETRIC DENSITY                                                                         ####

# function to rescale ppm values for easier interpretation of plots
scale.ppm <- function(px.ppm, scale.by = 128 * 1024) {
    
    fv <- predict(px.ppm)
    fv$v <- fv$v * scale.by
    fv
}


nonpara <- function(bpx.ppp, scale.by = 128*1024, sig = bw.diggle(bpx.ppp), return.model = F) {
    
    # bandwidth est using MSE approach
    np <- density(bpx.ppp, sigma = sig)
    
    # rescale intensity for easier interpretation
    np$v <- np$v * 128*1024
    
    plot(np, main = "")
    
    if (return.model) return(np)
}

lapply(ppp.off, bw.diggle)
lapply(ppp.off, nonpara, sig = 80)

lapply(ppp.cb, bw.diggle)
lapply(ppp.cb, nonpara, sig = 80)

lapply(ppp.nl, bw.diggle)
lapply(ppp.nl, nonpara, sig = 80)

lapply(ppp.scmd, bw.diggle)
lapply(ppp.scmd, nonpara, sig = 80)

####################################################################################################

# PARAMETRIC MODELLING                                                                          ####

# envelope plotting function
env.plot <- function(px.ppp, px.ppm, dist.fun, normalise = F, ...) {
    
    if (normalise) {
        trans <- expression(. - pi * r ** 2)
    } else {
        trans <- NULL
    }
    
    plot(envelope(px.ppm, dist.fun, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         col = "blue", legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), ...)
    plot(envelope(px.ppp, dist.fun, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, legend = F, shadecol = adjustcolor("gold", alpha = 0.2))
    
    legend("topleft", bty = "n",
           pt.bg = adjustcolor(c("gold", "cyan3", NA, NA, NA), alpha = 0.2), 
           col = c(NA, NA, "red", "blue", "black"),
           pch = c(22, 22, NA, NA, NA), 
           lty = c(NA, NA, 2, 2, 1),
           legend = c("Envelope for data", "Envelope for model", "Theoretical (data)", "Theoretical (model)", "Observed"))
}

# quadratic trend modelling over CB error points
{
    zz <- lapply(ppp.cb, ppm, ~ x + y + I(x^2) + I(y^2) + I(x *y))
    lapply(zz, plot, se = F, how = "contour")
    
    env.plot(ppp.cb$"131122", zz$"131122", Kest, normalise = T)
    env.plot(ppp.cb$"140128", zz$"140128", Kest, normalise = T)
    env.plot(ppp.cb$"160430", zz$"160430", Kest, normalise = T)
    env.plot(ppp.cb$"MCT225", zz$"MCT225", Kest, normalise = T)
    
    env.plot(ppp.cb$"131122", zz$"131122", Gest, normalise = F)
    env.plot(ppp.cb$"140128", zz$"140128", Gest, normalise = F)
    env.plot(ppp.cb$"160430", zz$"160430", Gest, normalise = F)
    env.plot(ppp.cb$"MCT225", zz$"MCT225", Gest, normalise = F)
    
    env.plot(ppp.cb$"131122", zz$"131122", Fest, normalise = F)
    env.plot(ppp.cb$"140128", zz$"140128", Fest, normalise = F)
    env.plot(ppp.cb$"160430", zz$"160430", Fest, normalise = F)
    env.plot(ppp.cb$"MCT225", zz$"MCT225", Fest, normalise = F)
}

# quadratic trend modelling over SC-nonuniform points
{
    zz.scmd <- lapply(ppp.scmd, ppm, ~ x + y + I(x^2) + I(y^2) + I(x *y))
    lapply(zz.scmd, plot, se = F, how = "contour")
    
    env.plot(ppp.cb$"131122", zz.scmd$"131122", Kest, normalise = T)
    env.plot(ppp.cb$"140128", zz.scmd$"140128", Kest, normalise = T)
    env.plot(ppp.cb$"160430", zz.scmd$"160430", Kest, normalise = T)
    env.plot(ppp.cb$"MCT225", zz.scmd$"MCT225", Kest, normalise = T)
    
    env.plot(ppp.cb$"131122", zz$"131122", Gest, normalise = F)
    env.plot(ppp.cb$"140128", zz$"140128", Gest, normalise = F)
    env.plot(ppp.cb$"160430", zz$"160430", Gest, normalise = F)
    env.plot(ppp.cb$"MCT225", zz$"MCT225", Gest, normalise = F)
    
    env.plot(ppp.cb$"131122", zz$"131122", Fest, normalise = F)
    env.plot(ppp.cb$"140128", zz$"140128", Fest, normalise = F)
    env.plot(ppp.cb$"160430", zz$"160430", Fest, normalise = F)
    env.plot(ppp.cb$"MCT225", zz$"MCT225", Fest, normalise = F)
}

lapply(lapply(ppp.nl, ppm, ~ x + y + I(x^2) + I(y^2) + I(x *y)), plot, se = F)

