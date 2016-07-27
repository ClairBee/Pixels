
# IN WHICH I SHALL FIT AN INHOMOGENEOUS POISSON PROCESS TO THE BAD PIXEL MAP

library("IO.Pixels"); library("CB.Misc")

# convert latest bad pixel map to point process
bp <- readRDS("./Notes/Final-classifications/fig/bad-px-by-feature-incl-local.rds")
bpx <- bp$"160430"
bpx <- bpx[bpx$f.type %in% c("cl.root", "singleton") & !bpx$type %in% c("l.bright", "l.dim"),]
bp.ppp <- ppp(bpx$row, bpx$col, c(1,1996), c(1,1996))

fpath <- "./Notes/IPP-fitting/fig/"
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "yellow", NA, "gold", "grey", NA, "blue", "blue", "green3")

####################################################################################################

# FUNCTIONS                                                                                     ####

scale.ppm <- function(px.ppm, scale.by = 128 * 1024) {
    
    fv <- predict(px.ppm)
    fv$v <- fv$v * scale.by
    fv
}

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


####################################################################################################

# QUADRAT TEST - INDEPENDENCE WITHIN SUBPANELS                                                  ####

# use 'step-up' approach: if any quadrat rejects CSR (after p-value adjustment), reject for whole panel
# could also use distance between kernel density & model as measure of GOF?

# test CSR across all quadrats
quadrat.test(bp.ppp, xbreaks = c(0:16) * 128 + 0.5, ybreaks = c(0:2) * 1024 + 0.5)
quadratcount(bp.ppp, xbreaks = c(0:16) * 128 + 0.5, ybreaks = c(0:2) * 1024 + 0.5)

quadrat.test(bp.ppp[owin(c(0,128) + 1 * 128, c(0,1024))], nx = 1, ny = 8)
hh <- quadrat.test(bp.ppp[owin(c(0,128) + 1 * 128, c(1024,2048))], nx = 1, ny = 8)

res <- array(dim = c(16, 2, 2), dimnames = list(NULL, c("l", "u"), c("statistic", "p.val")))

# expected values < 5, so use Monte Carlo test
for (ul in 0:1) {
    for (p in 0:15) {
        tmp <- quadrat.test(bp.ppp[owin(c(1,128) + p * 128, c(1,1024) + ul * 1024)], 
                            nx = 1, ny = 8, method = "MonteCarlo")
        res[p+1,ul+1,"statistic"] = tmp$statistic
        res[p+1,ul+1,"p.val"] = tmp$p.value
    }
}

plot(bp.ppp[owin(c(1,128) + 15 * 128, c(1,1024) + 0 * 1024)], pch = 20)
quadrat.test(bp.ppp[owin(c(1,128) + 1 * 128, c(1,1024) + 1 * 1024)], nx = 1, ny = 8)

plot(bp.ppp, pch = 20, main = "Raw p-values"); {
    rect(c(0:15) * 128 + 0.5, rep(1024.5, 16), c(1:16) * 128 + 0.5, rep(2048.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(res[,"u","p.val"] < 0.05) + 1])
    rect(c(0:15) * 128 + 0.5, rep(0.5, 16), c(1:16) * 128 + 0.5, rep(1024.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(res[,"l","p.val"] < 0.05) + 1])
}

plot(bp.ppp, pch = 20, main = "Holm correction"); {
    rect(c(0:15) * 128 + 0.5, rep(1024.5, 16), c(1:16) * 128 + 0.5, rep(2048.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,,"p.val"])[17:32] < 0.05) + 1])
    rect(c(0:15) * 128 + 0.5, rep(0.5, 16), c(1:16) * 128 + 0.5, rep(1024.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,"l","p.val"])[1:16] < 0.05) + 1])
}

plot(bp.ppp, pch = 20, main = "FDR correction"); {
    rect(c(0:15) * 128 + 0.5, rep(1024.5, 16), c(1:16) * 128 + 0.5, rep(2048.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,,"p.val"], method = "fdr")[17:32] < 0.05) + 1])
    rect(c(0:15) * 128 + 0.5, rep(0.5, 16), c(1:16) * 128 + 0.5, rep(1024.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,"l","p.val"], method = "fdr")[1:16] < 0.05) + 1])
}

plot(bp.ppp, pch = 20, main = "Hochberg correction"); {
    rect(c(0:15) * 128 + 0.5, rep(1024.5, 16), c(1:16) * 128 + 0.5, rep(2048.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,,"p.val"], method = "hochberg")[17:32] < 0.05) + 1])
    rect(c(0:15) * 128 + 0.5, rep(0.5, 16), c(1:16) * 128 + 0.5, rep(1024.5, 16), border = NA, 
         col = c(NA, adjustcolor("red", alpha = 0.4))[(p.adjust(res[,"l","p.val"], method = "hochberg")[1:16] < 0.05) + 1])
}

####################################################################################################

# NONPARAMETRIC - QUARTIC SMOOTHING                                                             ####

nonpara <- density(bp.ppp, sigma = bw.diggle(bp.ppp))    # bandwidth est using MSE approach
nonpara$v <- nonpara$v * 128*1024       # rescale intensity for easier interpretation
plot(nonpara)

pdf(paste0(fpath, "nonparametric-bw120.pdf")); {
    par(mar = c(0,0,0,0))
    contour(nonpara, main = "")
    draw.panels(col = "grey", lty = 2)
    points(bpx, pch = 20, cex = 0.8, col = Cat.cols[bpx$type])
    dev.off()
}

pdf(paste0(fpath, "nonparametric-bw120-image.pdf")); {
    par(mar = c(0,0,0,1))
    image(nonpara, main = "")
    draw.panels(lty = 2)
    points(bpx, pch = 20, cex = 0.8)
    dev.off()
    crop.pdf(paste0(fpath, "nonparametric-bw120-image.pdf"))
}

# plot(envelope(bp.ppp, Kinhom, sigma = bw.diggle, simulate = expression(rpoispp(nonpara))))

####################################################################################################

# QUADRATIC TREND                                                                               ####

qt.ppm <- ppm(bp.ppp ~ x + y + I(x^2) + I(y^2) + I(x *y))

pdf(paste0(fpath, "quadratic-trend.pdf")); {
    par(mar = c(0,0,0,0))
    contour(scale.ppm(qt.ppm), main = "")
    draw.panels(col = "grey", lty = 2)
    points(bpx, pch = 20, cex = 0.8, col = Cat.cols[bpx$type])
    dev.off()
    crop.pdf(paste0(fpath, "quadratic-trend.pdf"))
}

pdf(paste0(fpath, "quadratic-image.pdf")); {
    par(mar = c(0,0,0,1))
    image(scale.ppm(qt.ppm), main = "")
    draw.panels(lty = 2)
    points(bpx, pch = 20, cex = 0.8)
    dev.off()
    crop.pdf(paste0(fpath, "quadratic-image.pdf"))
}

pdf(paste0(fpath, "quadratic-trend-K.pdf"), width = 7, height = 4); {
    par(mar = c(4, 4, 1, 1))
    env.plot(bp.ppp, qt.ppm, Kest, normalise = T, main = "")
    dev.off()
}


env.plot(bp.ppp, qt.ppm, Fest, normalise = F, main = "")
env.plot(bp.ppp, qt.ppm, Gest, normalise = F, main = "")
env.plot(bp.ppp, qt.ppm, Jest, normalise = F, main = "")
env.plot(bp.ppp, qt.ppm, Kest, normalise = T, main = "")

# difference between model & kernel density
vv <- nonpara
vv$v <- (predict(qt.ppm)$v * 128 * 1024 - vv$v) 
plot(vv)

####################################################################################################

# QUADRAT TEST OF FITTED MODELS                                                                 ####

# CSR (raw data)
quadrat.test(bp.ppp, xbreaks = c(0:16) * 128 + 0.5, ybreaks = c(0:2) * 1024 + 0.5)
quadrat.test(qt.ppm, xbreaks = c(0:16) * 128 + 0.5, ybreaks = c(0:2) * 1024 + 0.5)
# does not reject fitted model. Hooray!

####################################################################################################

# SUBPANELS                                                                                     ####

# plot using only subpanels as covariates
sp.ppm <- ppm(bp.ppp ~ tt, covariates = list(tt = tess(xgrid = panel.edges()$x-0.5, ygrid = panel.edges()$y-0.5)))

pdf(paste0(fpath, "subpanel-flat-trend.pdf")); {
    par(mar = c(0,0,0,0))
    contour(scale.ppm(sp.ppm), main = "")
    draw.panels(col = "grey", lty = 2)
    points(bpx, pch = 20, cex = 0.8, col = Cat.cols[bpx$type])
    dev.off()
    crop.pdf(paste0(fpath, "subpanel-flat-trend.pdf"))
}

pdf(paste0(fpath, "subpanel-flat-image.pdf")); {
    par(mar = c(0,0,0,1))
    image(scale.ppm(sp.ppm), main = "")
    points(bpx, pch = 20, cex = 0.8)
    dev.off()
    crop.pdf(paste0(fpath, "subpanel-flat-image.pdf"))
}

pdf(paste0(fpath, "subpanel-flat-trend-K.pdf"), width = 7, height = 4); {
    par(mar = c(4, 4, 1, 1))
    env.plot(bp.ppp, sp.ppm, Kest, normalise = T, main = "")
    dev.off()
}

# difference between model & kernel density
vv <- nonpara
vv$v <- (predict(sp.ppm)$v * 128 * 1024 - vv$v) 
plot(vv)

####################################################################################################

# DISTANCE BETWEEN FITTED MODEL AND NONPARAMETRIC DENSITY                                       ####

nonpara <- density(bp.ppp, sigma = bw.diggle(bp.ppp))    # bandwidth est using MSE approach
qt.ppm <- ppm(bp.ppp ~ x + y + I(x^2) + I(y^2) + I(x *y))
qt.pred <- predict(qt.ppm)
diff <- qt.pred; diff$v <- qt.pred$v - nonpara$v
plot(diff)

log(sqrt(mean(diff$v^2)))

sp.pred <- predict(sp.ppm)
diff.sp <- sp.pred; diff.sp$v <- sp.pred$v - nonpara$v

log(sqrt(mean(diff.sp$v^2)))

fv.diff <- qt.pred; fv.diff$v <- qt.pred$v - sp.pred$v
plot(fv.diff)
plot(diff.sp)
plot(diff)
log(sqrt(mean(fv.diff$v^2)))


####################################################################################################

# SUBPANELS WITH GRADIENT                                                                       ####

spg.ppm <- ppm(bp.ppp ~ (x + y) * tt, 
               covariates = list(tt = tess(xgrid = panel.edges()$x-0.5, ygrid = panel.edges()$y-0.5)))

pdf(paste0(fpath, "subpanel-gradient-trend.pdf")); {
    par(mar = c(0,0,0,0))
    contour(scale.ppm(spg.ppm), main = "")
    draw.panels(col = "grey", lty = 2)
    points(bpx, pch = 20, cex = 0.8, col = Cat.cols[bpx$type])
    dev.off()
    crop.pdf(paste0(fpath, "subpanel-gradient-trend.pdf"))
}

pdf(paste0(fpath, "subpanel-gradient-image.pdf")); {
    par(mar = c(0,0,0,1))
    image(scale.ppm(spg.ppm), main = "")
    points(bpx, pch = 20, cex = 0.8)
    dev.off()
    crop.pdf(paste0(fpath, "subpanel-gradient-image.pdf"))
}

# difference between model & kernel density
vv <- nonpara
vv$v <- (predict(spg.ppm)$v * 128 * 1024 - vv$v) 
plot(vv)

####################################################################################################

# CSR PER SUBPANEL?                                                                             ####

# plot K-function per subpanel
pdf(paste0(fpath, "K-per-subpanel.pdf")); {
    par(mar = c(4, 4, 3, 1), mfrow = c(4, 4))
    for (ul in 1:2) {
        for (p in 1:16) {
            p.ref <- paste0(c("L","U")[ul], formatC(p, width = 2, flag = "0"))
            sp <- bpx[bpx$sp == p.ref,]
            plot(envelope(ppp(sp$row, sp$col, panel.edges()$x[c(p, p+1)]-0.5, panel.edges()$y[c(ul, ul+1)]-0.5), 
                          nsim = 99, nrank = 2, transform = expression(. - pi * r ** 2), fix.n = T),
                 main = p.ref)
        }
    }
    dev.off()
}

# mostly CSR: panels L8, U13, U16 not

plot(bpx[,1:2], pch = 20)
draw.panels(col = "grey")
rect(895, 1, 1022, 992, col = adjustcolor("cyan3", alpha = 0.1), border = NA)
rect(1919, 993, 1997, 1996, col = adjustcolor("cyan3", alpha = 0.1), border = NA)
rect(1535, 993, 1663, 1996, col = adjustcolor("cyan3", alpha = 0.1), border = NA)
                 
####################################################################################################

# VARIOUS REJECTED MODELS                                                                       ####

# linear x + y
{
    lin.ppm <- ppm(bp.ppp ~ x + y)
    
    contour(scale.ppm(lin.ppm), main = "")
    image(scale.ppm(lin.ppm), main = "")
    env.plot(bp.ppp, lin.ppm, Kest, normalise = T, main = "")
}

# point source at 1024, 1024: linear/quadratic/cubic gradient
{
    centre.spot <- function(x, y) {(1024 - x)^2 + (1024 - y)^2}
    circ.ppm.1 <- ppm(bp.ppp ~ centre.spot)
    circ.ppm.2 <- ppm(bp.ppp ~ poly(centre.spot, 2))
    circ.ppm.3 <- ppm(bp.ppp ~ poly(centre.spot, 3))
    
    par(mfrow = c(2, 3))
    contour(scale.ppm(circ.ppm.1), main = "")
    contour(scale.ppm(circ.ppm.2), main = "")
    contour(scale.ppm(circ.ppm.3), main = "")
    
    image(scale.ppm(circ.ppm.1), main = "")
    image(scale.ppm(circ.ppm.2), main = "")
    image(scale.ppm(circ.ppm.3), main = "")
    
    par(mfrow = c(1,1))
    plot(envelope(circ.ppm.1, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         main = "", col = "magenta3", shadecol = adjustcolor("pink", alpha = 0.2), legend = F)
    plot(envelope(circ.ppm.2, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, col = "darkgreen", shadecol = adjustcolor("chartreuse3", alpha = 0.2))
    plot(envelope(circ.ppm.3, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, shadecol = adjustcolor("gold", alpha = 0.2))
    plot(envelope(bp.ppp, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, shadecol = adjustcolor("cyan3", alpha = 0.2), main = "")
}


####################################################################################################


####################################################################################################

# COMPARE ALL K-FUNCTIONS                                                                       ####

# transects (mainly for reference)
{
    cc <- 120
    
    o.plot(nonpara$v[cc,], col = "blue")
    o.plot(scale.ppm(qt.ppm)$v[cc,], add = T, col = "darkgreen")
    o.plot(scale.ppm(sp.ppm)$v[cc,], add = T, col = "red")
    o.plot(scale.ppm(tmp.ppm)$v[cc,], add = T, col = "gold")
    
}

# normalised K-function
pdf(paste0(fpath, "normed-K.pdf")); {
    trans <- expression((. - pi * r ** 2))
    
    par(mar = c(4, 4, 1, 1))
    
    plot(envelope(bp.ppp, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         col = "blue", legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), main = "")
    plot(envelope(qt.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, col = "darkgreen", shadecol = adjustcolor("chartreuse3", alpha = 0.2))
    plot(envelope(sp.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, shadecol = adjustcolor("gold", alpha = 0.2))
    plot(envelope(spg.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, col = "magenta3", shadecol = adjustcolor("pink", alpha = 0.2))
    
    legend("topleft", bty = "n",
           col = c("blue", "darkgreen", "red", "magenta3", "black"),
           lty = c(2, 2, 2, 2, 1), 
           fill = adjustcolor(c("cyan3", "chartreuse3", "gold", "pink", NA), alpha = 0.4), border = NA,
           legend = c("CSR", "Quadratic trend", "Subpanels (flat)", "Subpanels (linear)", "Observed"))
 
    dev.off()
}

# focus on small scales (0-100)
pdf(paste0(fpath, "normed-K-small-scale.pdf")); {
    trans <- expression((. - pi * r ** 2))
    
    par(mar = c(4, 4, 1, 1))
    
    plot(envelope(bp.ppp, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         xlim = c(0,100), col = "blue", legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), main = "")
    plot(envelope(qt.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, col = "darkgreen", shadecol = adjustcolor("chartreuse3", alpha = 0.2))
    plot(envelope(sp.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, shadecol = adjustcolor("gold", alpha = 0.2))
    plot(envelope(spg.ppm, Kest, nsim = 99, nrank = 2, transform = trans, verbose = F, fix.n = T), 
         add = T, col = "magenta3", shadecol = adjustcolor("pink", alpha = 0.2))
    
    legend("topleft", bty = "n",
           col = c("blue", "darkgreen", "red", "magenta3", "black"),
           lty = c(2, 2, 2, 2, 1), 
           fill = adjustcolor(c("cyan3", "chartreuse3", "gold", "pink", NA), alpha = 0.4), border = NA,
           legend = c("CSR", "Quadratic trend", "Subpanels (flat)", "Subpanels (linear)", "Observed"))
    
    
    dev.off()
}

# G-function
pdf(paste0(fpath, "G.pdf")); {
    trans <- expression((. - pi * r ** 2))
    
    par(mar = c(4, 4, 1, 1))
    
    plot(envelope(bp.ppp, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), xlim = c(0,20), 
         col = "blue", legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), main = "")
    plot(envelope(qt.ppm, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         add = T, col = "darkgreen", shadecol = adjustcolor("chartreuse3", alpha = 0.2))
    plot(envelope(sp.ppm, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         add = T, shadecol = adjustcolor("gold", alpha = 0.2))
    plot(envelope(spg.ppm, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         add = T, col = "magenta3", shadecol = adjustcolor("pink", alpha = 0.2))
    
    legend("topleft", bty = "n",
           col = c("blue", "darkgreen", "red", "magenta3", "black"),
           lty = c(2, 2, 2, 2, 1), 
           fill = adjustcolor(c("cyan3", "chartreuse3", "gold", "pink", NA), alpha = 0.4), border = NA,
           legend = c("CSR", "Quadratic trend", "Subpanels (flat)", "Subpanels (linear)", "Observed"))
    
    dev.off()
}

####################################################################################################

####################################################################################################

# CLUSTER MODEL - MATERN                                                                        ####

bpx <- bp$"160430"
bpx <- bpx[!bpx$type %in% c("l.bright", "l.dim", "line.b", "line.d"),]
# rescale ppp for easier interpretation of parameters
bp.ppp <- ppp(bpx$row / 1996, bpx$col / 1996, c(0,1), c(0,1))

kppm.mat <- kppm(bp.ppp ~ x + y + I(x^2) + I(y^2) + I(x *y),
                 clusters = c("MatClust"))
plot(kppm.mat)
summary(kppm.mat)

# envelope plots
{
    plot(envelope(kppm.mat, Fest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    
    plot(envelope(kppm.mat, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    
    plot(envelope(kppm.mat, Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    plot(envelope(kppm.mat, Kest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), xlim = c(0,100))
}

ipp <- function(x, y, fitted.ppm = qt.ppm) {
    eval(parse(text = paste(coef(fitted.ppm)[1],
          apply(cbind("+", coef(fitted.ppm)[-1], "*", 
                names(coef(fitted.ppm))[-1]),
          1, paste, collapse = ""), collapse = "")))
}

plot(envelope(bp.ppp, Fest, nrank = 2, verbose = F, 
              simulate = expression(rMatClust(kappa = ipp, scale = 1/1996, mu = 99)), nsim = 99))

ppm.mat <- as.ppm(kppm.mat)
clusterkernel(kppm.mat)     # sqrt(x^2 + y^2)
clusterradius(kppm.mat)     # 4.567697

parameters(kppm.mat.2)

####################################################################################################

# CLUSTER MODEL - THOMAS                                                                        ####

kppm.thom <- kppm(bp.ppp ~ x + y + I(x^2) + I(y^2) + I(x *y),
                 clusters = "Thomas")
plot(kppm.thom)

# envelope plots
{
    plot(envelope(kppm.thom, Fest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    
    plot(envelope(kppm.thom, Gest, nsim = 99, nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    
    plot(envelope(kppm.thom, Kest, nsim = 99, transform = expression(. - pi * r ** 2), 
                  nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2))
    plot(envelope(kppm.thom, Kest, nsim = 99, transform = expression(. - pi * r ** 2), 
                  nrank = 2, verbose = F, fix.n = T), 
         legend = F, shadecol = adjustcolor("cyan3", alpha = 0.2), xlim = c(0,100))
}

####################################################################################################

# PAIRWISE CLUSTERING                                                                           ####

####################################################################################################

# OLD DATA                                                                                      ####

# convert latest bad pixel map to point process
bpx <- readRDS("./Other-data/Old-data/bad-px-features-131122.rds")

bpx <- bpx[bpx$f.type %in% c("cl.root", "singleton") & !bpx$type %in% c("edge", "l.bright", "l.dim", "s.dim", "s.bright"),]
bp.ppp <- ppp(bpx$row, bpx$col, c(1,2000), c(200,1800))

# nonparametric
nonpara <- density(bp.ppp, sigma = bw.diggle(bp.ppp))    # bandwidth est using MSE approach
nonpara$v <- nonpara$v * 128*1024       # rescale intensity for easier interpretation
plot(nonpara)

####################################################################################################

# MODEL FITTING WITH PADDED ARRAY                                                               ####

qq <- apply(which(!is.na(pwm[,,"black", "160430"]), arr.ind = T), 2, range)
attr(bp$"160430", "crop.region") <- list(x = qq[,1], y = qq[,2])
attr(bp$"160430", "array.dim") <- c(2048, 2048)

bpx <- bp$"160430"
bpx <- bpx[!bpx$type %in% c("line.b", "line.d", "l.bright", "l.dim"),]
bpx <- bpx[!bpx$f.type %in% c("line.body", "cl.body"),]

crop.ppp <- ppp(bpx$row + 2, bpx$col + 32, 
                attr(bpx, "crop.region")$x, attr(bpx, "crop.region")$y)

nonpara <- density(crop.ppp, sigma = bw.diggle(crop.ppp))    # bandwidth est using MSE approach
nonpara$v <- nonpara$v * 128*1024       # rescale intensity for easier interpretation
plot(nonpara)
rect(0.5, 0.5, 2047.5, 2047.5)          # add border showing true panel edge
draw.panels(p = panel.edges(left.crop = 0, upper.crop = 0, x.dim = 2048, y.dim = 2048))

# can use this approach to assign attribute to each bad pixel map & track original image dimensions

####################################################################################################
