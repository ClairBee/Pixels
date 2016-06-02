
library("IO.Pixels"); library("CB.Misc")

d <- load.daily("160430")

load.pixel.means()
load.pixel.sds()

bp <- readRDS("./Notes/Standard-deviations/bad-px-maps-with-noise.rds")
fpath <- "./Notes/Standard-deviations/fig2/"

####################################################################################################

# TOTAL VARIATION                                                                               ####

tv <- apply(d[,,,"black"], 1:2, function(x) sum(abs(c(x, NA) - c(NA, x))[2:20]))

hist(tv, breaks = "fd", xlim = c(0,2000))

hist(pw.sd[,,"black", "160430"], breaks = "fd", xlim = c(0,200))

s <- sample(1:1996^2, 100000, replace = F)

plot(pw.sd[,,"black", "160430"][s], tv[s], pch = "*")

cor(pw.sd[,,"black", "160430"][s], tv[s])   # 0.927: strongly correlated with SD

####################################################################################################

# EXPLORATORY PLOTS OF SD                                                                       ####

# plot highest 20 SDs from each channel, and lowest 20

noisy <- unique(rbind(which(pw.sd[,,"black", "160430"] > sort(pw.sd[,,"black", "160430"], decreasing = T)[21], arr.ind = TRUE),
               which(pw.sd[,,"grey", "160430"] > sort(pw.sd[,,"grey", "160430"], decreasing = T)[21], arr.ind = TRUE),
               which(pw.sd[,,"white", "160430"] > sort(pw.sd[,,"white", "160430"], decreasing = T)[21], arr.ind = TRUE)))

quiet <- unique(rbind(which(pw.sd[,,"black", "160430"] > 0 & pw.sd[,,"black", "160430"] < sort(pw.sd[,,"black", "160430"][pw.sd[,,"black", "160430"] > 0], decreasing = F)[21], arr.ind = TRUE),
                      which(pw.sd[,,"grey", "160430"] > 0 & pw.sd[,,"grey", "160430"] < sort(pw.sd[,,"grey", "160430"][pw.sd[,,"grey", "160430"] > 0], decreasing = F)[21], arr.ind = TRUE),
                      which(pw.sd[,,"white", "160430"] > 0 & pw.sd[,,"white", "160430"] < sort(pw.sd[,,"white", "160430"][pw.sd[,,"white", "160430"] > 0], decreasing = F)[21], arr.ind = TRUE)))


zz.n <- round(rbind.fill(apply(noisy, 1, 
                       function(p) data.frame(b.mean = mean(d[p[1],p[2],,"black"]),
                                              b.median = median(d[p[1],p[2],,"black"]),
                                              n.median = median(apply(d[,,,"black"], 3, function(x) x[matrix(c(p[1] + rep(-1:1, 3), p[2] + sort(rep(-1:1, 3))), ncol = 2)])),
                                              b.sd = sd(d[p[1],p[2],,"black"]),
                                              b.mad = mad(d[p[1],p[2],,"black"]),
                                              b.min = min(d[p[1],p[2],,"black"]),
                                              b.max = max(d[p[1],p[2],,"black"]),
                                              b.tv = sum(abs(c(d[p[1],p[2],,"black"], NA) - c(NA, d[p[1],p[2],,"black"]))[2:20])))),0)
      
zz.q <- round(rbind.fill(apply(quiet, 1, 
                               function(p) data.frame(b.mean = mean(d[p[1],p[2],,"black"]),
                                                      b.median = median(d[p[1],p[2],,"black"]),
                                                      n.median = median(apply(d[,,,"black"], 3, function(x) x[matrix(c(p[1] + rep(-1:1, 3), p[2] + sort(rep(-1:1, 3))), ncol = 2)])),
                                                      b.sd = sd(d[p[1],p[2],,"black"]),
                                                      b.mad = mad(d[p[1],p[2],,"black"]),
                                                      b.min = min(d[p[1],p[2],,"black"]),
                                                      b.max = max(d[p[1],p[2],,"black"]),
                                                      b.tv = sum(abs(c(d[p[1],p[2],,"black"], NA) - c(NA, d[p[1],p[2],,"black"]))[2:20])))),0)

             
pdf(paste0(fpath, "ts-plots-noisy-black.pdf")); {
    par(mfrow = c(5, 4), mar = c(2, 2, 1, 1))
    for (i in 1:nrow(noisy)) {
        m <- median(d[noisy[i,1],noisy[i,2],,"black"])
        
        o.plot(d[noisy[i,1],noisy[i,2],,"black"], ylim = m + c(-5000, 5000))
#        rect(-1, m - s, 21, m + s, col = adjustcolor("gold", alpha = 0.2), border = NA)
        abline(h = m, col = "red")
        mm <- median(apply(d[,,,"black"], 3, function(x) x[matrix(c(noisy[i,1] + rep(-1:1, 3), noisy[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
#        rect(-1, mm - s, 21, mm + s, col = adjustcolor("cyan3", alpha = 0.2), border = NA, density = 20)
        abline(h = mm, col =  "cyan3")
        
   }
    dev.off()
}

s <- sd(d[,,,"black"])

pdf(paste0(fpath, "ts-plots-quiet-black.pdf")); {
    par(mfrow = c(5, 3), mar = c(2, 2, 1, 1))
    for (i in 1:nrow(quiet)) {
        {
            m <- median(d[quiet[i,1],quiet[i,2],,"black"])
            
            o.plot(d[quiet[i,1],quiet[i,2],,"black"], ylim = m + c(-5000, 5000))
            #        rect(-1, m - s, 21, m + s, col = adjustcolor("gold", alpha = 0.2), border = NA)
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"black"], 3, function(x) x[matrix(c(quiet[i,1] + rep(-1:1, 3), quiet[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            #        rect(-1, mm - s, 21, mm + s, col = adjustcolor("cyan3", alpha = 0.2), border = NA)
            abline(h = mm, col =  "cyan3")
            text(0.5, m - 4500, paste0("SD: ", round(sd(d[quiet[i,1],quiet[i,2],,"black"]), 1)), pos = 4)
        }
        {
            m <- median(d[quiet[i,1],quiet[i,2],,"grey"])
            o.plot(d[quiet[i,1],quiet[i,2],,"grey"], ylim = m + c(-5000, 5000))
            abline(h = m, col = "red")
            abline(h = median(apply(d[,,,"grey"], 3, function(x) x[matrix(c(quiet[i,1] + rep(-1:1, 3), quiet[i,2] + sort(rep(-1:1, 3))), ncol = 2)])), col =  "cyan3")
            text(0.5, m - 4500, paste0("SD: ", round(sd(d[quiet[i,1],quiet[i,2],,"grey"]), 1)), pos = 4)
        }
        {
            m <- median(d[quiet[i,1],quiet[i,2],,"white"])
            o.plot(d[quiet[i,1],quiet[i,2],,"white"], ylim = m + c(-5000, 5000))
            abline(h = m, col = "red")
            abline(h = median(apply(d[,,,"white"], 3, function(x) x[matrix(c(quiet[i,1] + rep(-1:1, 3), quiet[i,2] + sort(rep(-1:1, 3))), ncol = 2)])), col =  "cyan3")
            text(0.5, m - 4500, paste0("SD: ", round(sd(d[quiet[i,1],quiet[i,2],,"white"]), 1)), pos = 4)
        }
    }
    dev.off()
}

acf(d[197:199,84:86,,"black"])

####################################################################################################

# POINTS NOT OTHERWISE IDENTIFIED                                                               ####



bb <- cbind(bp$"160430"[bp$"160430"$type == "noisy.b",1:2], 
            sd.b = pw.sd[,,"black", "160430"][as.matrix(bp$"160430"[bp$"160430"$type == "noisy.b",1:2])])

bb <- bb[rev(order(bb$sd.b)),]
bb[1:length(which(bp$"160430"$type == "noisy.gw")),]

px <- rbind(bp$"160430"[bp$"160430"$type == "noisy.gw",1:2],
            bb[1:length(which(bp$"160430"$type == "noisy.gw")),1:2])

gp <- sample.healthy(bp$"160430", n = nrow(px))

{
    zz.px <- round(rbind.fill(apply(px[,1:2], 1, 
                                    function(p) data.frame(b.mean = mean(d[p[1],p[2],,"black"]),
                                                           b.median = median(d[p[1],p[2],,"black"]),
                                                           n.median = median(apply(d[,,,"black"], 3, function(x) x[matrix(c(p[1] + rep(-1:1, 3), p[2] + sort(rep(-1:1, 3))), ncol = 2)])),
                                                           b.sd = sd(d[p[1],p[2],,"black"]),
                                                           b.mad = mad(d[p[1],p[2],,"black"]),
                                                           b.min = min(d[p[1],p[2],,"black"]),
                                                           b.max = max(d[p[1],p[2],,"black"]),
                                                           b.tv = sum(abs(c(d[p[1],p[2],,"black"], NA) - c(NA, d[p[1],p[2],,"black"]))[2:20])))),0)
    
    zz.gp <- round(rbind.fill(apply(gp[,1:2], 1, 
                                    function(p) data.frame(b.mean = mean(d[p[1],p[2],,"black"]),
                                                           b.median = median(d[p[1],p[2],,"black"]),
                                                           n.median = median(apply(d[,,,"black"], 3, function(x) x[matrix(c(p[1] + rep(-1:1, 3), p[2] + sort(rep(-1:1, 3))), ncol = 2)])),
                                                           b.sd = sd(d[p[1],p[2],,"black"]),
                                                           b.mad = mad(d[p[1],p[2],,"black"]),
                                                           b.min = min(d[p[1],p[2],,"black"]),
                                                           b.max = max(d[p[1],p[2],,"black"]),
                                                           b.tv = sum(abs(c(d[p[1],p[2],,"black"], NA) - c(NA, d[p[1],p[2],,"black"]))[2:20])))),0)
    
}

yr <- 4000
lim.b <- median(pw.sd[,,"black","160430"]) + 6 * sd(pw.sd[,,"black","160430"])
lim.g <- median(pw.sd[,,"grey","160430"]) + 6 * sd(pw.sd[,,"grey","160430"])
lim.w <- median(pw.sd[,,"white","160430"]) + 6 * sd(pw.sd[,,"white","160430"])

pdf(paste0(fpath, "ts-plots-noisy-px.pdf")); {
    par(mfrow = c(5, 3), mar = c(2, 2, 1, 1))
    for (i in 1:nrow(px)) {
        {
            m <- median(d[px[i,1],px[i,2],,"black"])
            o.plot(d[px[i,1],px[i,2],,"black"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"black"], 3, function(x) x[matrix(c(px[i,1] + rep(-1:1, 3), px[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[px[i,1],px[i,2],,"black"])
            if (s > lim.b) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[px[i,1],px[i,2],,"black"], NA) - c(NA, d[px[i,1],px[i,2],,"black"]))[2:20]), 1)), pos = 2)
        }
        {
            m <- median(d[px[i,1],px[i,2],,"grey"])
            o.plot(d[px[i,1],px[i,2],,"grey"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"grey"], 3, function(x) x[matrix(c(px[i,1] + rep(-1:1, 3), px[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[px[i,1],px[i,2],,"grey"])
            if (s > lim.g) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[px[i,1],px[i,2],,"grey"], NA) - c(NA, d[px[i,1],px[i,2],,"grey"]))[2:20]), 1)), pos = 2)
        }
        {
            m <- median(d[px[i,1],px[i,2],,"white"])
            o.plot(d[px[i,1],px[i,2],,"white"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"white"], 3, function(x) x[matrix(c(px[i,1] + rep(-1:1, 3), px[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[px[i,1],px[i,2],,"white"])
            if (s > lim.w) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[px[i,1],px[i,2],,"white"], NA) - c(NA, d[px[i,1],px[i,2],,"white"]))[2:20]), 1)), pos = 2)
        }
    }
    dev.off()
}



pdf(paste0(fpath, "ts-plots-normal-px.pdf")); {
    par(mfrow = c(5, 3), mar = c(2, 2, 1, 1))
    for (i in 1:nrow(gp)) {
        {
            m <- median(d[gp[i,1],gp[i,2],,"black"])
            o.plot(d[gp[i,1],gp[i,2],,"black"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"black"], 3, function(x) x[matrix(c(gp[i,1] + rep(-1:1, 3), gp[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[gp[i,1],gp[i,2],,"black"])
            if (s > lim.b) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[gp[i,1],gp[i,2],,"black"], NA) - c(NA, d[gp[i,1],gp[i,2],,"black"]))[2:20]), 1)), pos = 2)
        }
        {
            m <- median(d[gp[i,1],gp[i,2],,"grey"])
            o.plot(d[gp[i,1],gp[i,2],,"grey"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"grey"], 3, function(x) x[matrix(c(gp[i,1] + rep(-1:1, 3), gp[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[gp[i,1],gp[i,2],,"grey"])
            if (s > lim.g) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[gp[i,1],gp[i,2],,"grey"], NA) - c(NA, d[gp[i,1],gp[i,2],,"grey"]))[2:20]), 1)), pos = 2)
        }
        {
            m <- median(d[gp[i,1],gp[i,2],,"white"])
            o.plot(d[gp[i,1],gp[i,2],,"white"], ylim = m + c(-yr, yr))
            abline(h = m, col = "red")
            mm <- median(apply(d[,,,"white"], 3, function(x) x[matrix(c(gp[i,1] + rep(-1:1, 3), gp[i,2] + sort(rep(-1:1, 3))), ncol = 2)]))
            abline(h = mm, col =  "cyan3")
            s <- sd(d[gp[i,1],gp[i,2],,"white"])
            if (s > lim.w) {cc <- "red"} else {cc <- "black"}
            text(0.5, m - 3500, paste0("SD: ", round(s, 1)), pos = 4, col = cc)
            text(20.5, m - 3500, paste0("TV: ", round(sum(abs(c(d[gp[i,1],gp[i,2],,"white"], NA) - c(NA, d[gp[i,1],gp[i,2],,"white"]))[2:20]), 1)), pos = 2)
        }
    }
    dev.off()
}

# cutoff for 'noisy' by 6-sigma method: black 92, grey 341, white 633
apply(pw.sd[,,,"160430"], 3, function(x) median(x) + 6 * sd(x))

# AUTOCORRELATION                                                                               ####

a <- apply(d[,,,"black"], 1:2, function(x) acf(x, plot = F)$acf)

apply(a[2:14,,], 1, mean, na.rm = T)
# -0.03071157 -0.03708759 -0.03954781 -0.03989383 -0.03960967 -0.03846892 -0.03679143 -0.03502878 -0.03258498
# -0.03034245 -0.02767410 -0.02482435 -0.02184372
# significance level ~ 2/sqrt(21) ~ 0.436

a <- apply(d[,,,"grey"], 1:2, function(x) acf(x, plot = F)$acf)
apply(a[2:14,,], 1, mean, na.rm = T)
# -0.04665190 -0.04489819 -0.04316645 -0.04119162 -0.03899877 -0.03702583 -0.03351993 -0.03135371 -0.02974974
# -0.02698320 -0.02435307 -0.02255749 -0.01923282

pixel.image(a[2,,])
# no sign of anything in particular

a <- apply(d[,,,"black"], 1:2, function(x) pacf(x, plot = F)$acf)
apply(a[2:14,,], 1, mean, na.rm = T)

plot(apply(a[2:13,,], 1, mean, na.rm = T), ylim = c(-1,1))

p <- apply(d[1:100,1:100,,"black"], 1:2, function(x) pacf(x, plot = F)$acf)

####################################################################################################

# XTERISTICS OF PERSISTENTLY NOISY PIXELS                                                       ####

# get all noisy px
{
    # identify all 6-sigma noisy pixels
    npx <- apply(pw.sd, 4, apply, 3,
                 function (x) data.frame(which(x > median(x) + 6 * sd(x), arr.ind = T),
                                         type = "noisy"))
    
    # convert to array for easier manipulation
    npx.im <- abind(lapply(lapply(npx, lapply, bpx2im), abind, along = 3), along = 4)
    
    # get persistently noisy pixels from black images: 147 pixels
    p.b <- apply(npx.im[,,"black",], 1:2, sum)
    npx.b <- which(p.b == 12, arr.ind = T)
    
    apply(npx.im[,,"black",], 3, sum)
    
    # 141009 141118 141217 150108 150113 150126 150529 150730 150828 151015 160314 160430 
    #   1037    899   1168    567   1724   1770   1378   2167   2472   2770   3206   4579 
    
    # no strong clustering: usual spatial characteristics
    plot(npx.b, pch = 15, cex = 0.4, asp = T, xlim = c(0,1996), ylim = c(0,1996))
}

# load all daily images (only keeping noisy px, to save memory)
{
    dd <- array(dim = c(147, 20, 3, 12))
    for (i in 1:12) {
        dd[,,,i] <- array(apply(load.daily(dimnames(pw.sd)[[4]][i]), 4, apply, 3, "[", npx.b), dim = c(147, 20, 3))
    }
    dimnames(dd) <- list(NULL, NULL, dimnames(pw.sd)[[3]], dimnames(pw.sd)[[4]])
    saveRDS(dd, paste0(fpath, "persisting-noisy-px.rds"))
}

dd <- readRDS(paste0(fpath, "persisting-noisy-px.rds"))

# plots of transects
{
    mp.cols <- c(rep("green3", 2), rep("cyan3", 2), rep("blue", 2), 
                 rep("magenta3", 2), rep("orange", 2), rep("gold", 2))
    
    pdf(paste0(fpath, "persistent-noisy-px-black.pdf")); {
        par(mar = c(2, 4, 1, 1), mfrow = c(5, 4))
        for (i in 1:nrow(npx.b))     {
            tmp <- dd[i,,"black", ]
            if (max(pretty(tmp)) - min(pretty(tmp)) < 6000) {
                r <- c(min(pretty(tmp)), min(pretty(tmp)) + 6000)
            } else {
                r <- range(pretty(tmp))
            }
            matplot(tmp, type = "l", col = mp.cols, lty = c(4, 1), ylim = r,
                    ylab = paste(c("(", npx.b[i,1], ", ", npx.b[i,2], ")"), collapse = ""))
        }
        dev.off()
    }

    # check classification of persistently noisy pixels
    {
        lapply(bp, function(x) x[x[,1] == 188 & x[,2] == 1542,])    # noisy px at 188, 1542
        lapply(bp, function(x) x[x[,1] == 479 & x[,2] == 84,])      # noisy px at 479, 84
        lapply(bp, function(x) x[x[,1] == 1860 & x[,2] == 42,])     # noisy px at 1860, 42
        lapply(bp, function(x) x[x[,1] == 411 & x[,2] == 161,])     # noisy px at 1860, 42
        lapply(bp, function(x) x[x[,1] == 583 & x[,2] == 217,])     # noisy / l. bright px at 583, 217
        lapply(bp, function(x) x[x[,1] == 183 & x[,2] == 237,])     # noisy / l. bright px at 183, 237
        lapply(bp, function(x) x[x[,1] == 280 & x[,2] == 285,])     # l. bright
        lapply(bp, function(x) x[x[,1] == 57 & x[,2] == 405,])      # l. bright
        lapply(bp, function(x) x[x[,1] == 319 & x[,2] == 508,])     # l. bright > bright > v.bright
        lapply(bp, function(x) x[x[,1] == 27 & x[,2] == 546,])      # l. bright
        lapply(bp, function(x) x[x[,1] == 4 & x[,2] == 631,])       # l. bright > bright > v.bright
    }
}


# get key xteristics of all pixels identified as noisy
{
    # coords & counts of pixels ever identified as noisy
    p <- count(unique(which(npx.im == 1, arr.ind = T)[,c(1:2, 4)])[,1:2])
    pp <- unique(p[,1:2])
    
    dd.all <- array(dim = c(9831, 20, 3, 12))
    for (i in 1:12) {
        dd.all[,,,i] <- array(apply(load.daily(dimnames(pw.sd)[[4]][i]), 4, apply, 3, "[", as.matrix(pp)), dim = c(9831, 20, 3))
    }
    dimnames(dd.all) <- list(NULL, NULL, dimnames(pw.sd)[[3]], dimnames(pw.sd)[[4]])
    saveRDS(dd.all, paste0(fpath, "all-noisy-px.rds"))
}
# BUT where do these persistently noisy pixels sit in terms of non-persistent pixels?

# histograms of key xteristics, split by 'all px' / 7+ images / 10+ images. Compare distribution
tv <- function(x) sum(abs(c(x, NA) - c(NA, x))[2:20])
linear.coef <- function(x) line(x)$coef[2]
qfit.l <- function(z) lm(x ~ y + y2, data.frame(x = z, y = 1:20, y2 = 1:20^2))$coef["y2"]
qfit.q <- function(z) lm(x ~ y + y2, data.frame(x = z, y = 1:20, y2 = 1:20^2))$coef["y"]

qq <- data.frame(n = p$freq,
                 sd.b = apply(dd.all[,,"black","160430"], 1, sd),
                 sd.g = apply(dd.all[,,"grey","160430"], 1, sd),
                 sd.w = apply(dd.all[,,"white","160430"], 1, sd),
                 tv.b = apply(dd.all[,,"black", "160430"], 1, tv),
                 tv.g = apply(dd.all[,,"grey", "160430"], 1, tv),
                 tv.w = apply(dd.all[,,"white", "160430"], 1, tv),
                 mean.b = apply(dd.all[,,"black", "160430"], 1, mean),
                 mean.g = apply(dd.all[,,"grey", "160430"], 1, mean),
                 mean.w = apply(dd.all[,,"white", "160430"], 1, mean),
                 median.b = apply(dd.all[,,"black", "160430"], 1, median),
                 median.g = apply(dd.all[,,"grey", "160430"], 1, median),
                 median.w = apply(dd.all[,,"white", "160430"], 1, median),
                 min.b = apply(dd.all[,,"black", "160430"], 1, min),
                 min.g = apply(dd.all[,,"grey", "160430"], 1, min),
                 min.w = apply(dd.all[,,"white", "160430"], 1, min),
                 max.b = apply(dd.all[,,"black", "160430"], 1, max),
                 max.g = apply(dd.all[,,"grey", "160430"], 1, max),
                 max.w = apply(dd.all[,,"white", "160430"], 1, max),
                 lc.b = apply(dd.all[,,"black", "160430"], 1, linear.coef),
                 lc.g = apply(dd.all[,,"grey", "160430"], 1, linear.coef),
                 lc.w = apply(dd.all[,,"white", "160430"], 1, linear.coef),
                 qlc.b = apply(dd.all[,,"black", "160430"], 1, qfit.l),
                 qlc.g = apply(dd.all[,,"grey", "160430"], 1, qfit.l),
                 qlc.w = apply(dd.all[,,"white", "160430"], 1, qfit.l),
                 qqc.b = apply(dd.all[,,"black", "160430"], 1, qfit.q),
                 qqc.g = apply(dd.all[,,"grey", "160430"], 1, qfit.q),
                 qqc.w = apply(dd.all[,,"white", "160430"], 1, qfit.q))

# used for EEG clustering (Siuly et al, 2011):
#   minimum; maximum; mean; median; modus; first quartile; third quartile; inter-quartile range; standard deviation

# used by Barton:
#   mean; minimum; maximum; linear coefficient; quadratic coefficient

# histograms of various characteristics
{
    # custom histogram function to speed this up
    c.hist <- function(qq.col, d = 4, ...) {
        
        b <- c(floor(min(qq.col, 0)/d):ceiling(max(qq.col)/d))*d
        hh.low <- hist(qq.col[qq$n <= 6], breaks = b, plot = F)
        hh.mid <- hist(qq.col[qq$n > 6 & qq$n <= 9], breaks = b, plot = F)
        hh.high <- hist(qq.col[qq$n > 9], breaks = b, plot = F)
        
        hh.low$counts <- hh.low$density
        hh.mid$counts <- -hh.mid$density
        hh.high$counts <- -hh.mid$density
        
        plot(hh.low, ylim = c(min(hh.mid$counts, hh.high$counts), max(hh.low$counts)), col = "black", ...)
        plot(hh.mid, add = T, col = "gold", border = "gold")
        plot(hh.high, add = T, col = adjustcolor("red", alpha = 0.4), border = adjustcolor("red", 0.4))
        abline(h = 0)
        
        legend("bottomright", col = c("black", "gold", "red"), pch = 15, 
               legend = c("n <= 6", "6 < n <= 9", "n > 9"), bty = "n", cex = 0.6)
    }
    
    pdf(paste0(fpath, "xteristic-plots.pdf")); {
        par(mfrow = c(4, 3), mar = c(2, 2, 3, 1))
        # standard deviations
        {
            c.hist(qq$sd.b, xlim = c(0,1000), main = "Black SD", xlab = "", ylab = "")
            c.hist(qq$sd.g, xlim = c(0,1000), main = "Grey SD", xlab = "", ylab = "")
            c.hist(qq$sd.w, xlim = c(0,1000), main = "White SD", xlab = "", ylab = "")
        }
        
        # total variance
        {
            c.hist(qq$tv.b, xlim = c(0,25000), d = 10, main = "Black total variance", xlab = "", ylab = "")
            c.hist(qq$tv.g, xlim = c(0,25000), d = 10, main = "Grey total variance", xlab = "", ylab = "")
            c.hist(qq$tv.w, xlim = c(0,25000), d = 10, main = "White total variance", xlab = "", ylab = "")
        }
        
        # median value
        {
            c.hist(qq$median.b, xlim = c(0,25000), d = 10, main = "Black median", xlab = "", ylab = "")
            c.hist(qq$median.g, xlim = c(10000,40000), d = 10, main = "Grey median", xlab = "", ylab = "")
            c.hist(qq$median.w, xlim = c(35000, 65535), d = 10, main = "White median", xlab = "", ylab = "")
        }
        
        # minimum value
        {
            c.hist(qq$min.b, xlim = c(0,25000), d = 10, main = "Black minimum", xlab = "", ylab = "")
            c.hist(qq$min.g, xlim = c(10000,40000), d = 10, main = "Grey minimum", xlab = "", ylab = "")
            c.hist(qq$min.w, xlim = c(35000, 65535), d = 10, main = "White minimum", xlab = "", ylab = "")
        }
        
        # maximum value
        {
            c.hist(qq$max.b, xlim = c(0,25000), d = 10, main = "Black maximum", xlab = "", ylab = "")
            c.hist(qq$max.g, xlim = c(10000,40000), d = 10, main = "Grey maximum", xlab = "", ylab = "")
            c.hist(qq$max.w, xlim = c(35000, 65535), d = 10, main = "White maximum", xlab = "", ylab = "")
        }
        
        # linear coefficient (gradient)
        {
            c.hist(qq$lc.b, xlim = c(-100, 100), d = 1, main = "Black linear gradient", xlab = "", ylab = "")
            c.hist(qq$lc.g, xlim = c(-100, 100), d = 1, main = "Grey linear gradient", xlab = "", ylab = "")
            c.hist(qq$lc.w, xlim = c(-100, 100), d = 1, main = "White linear gradient", xlab = "", ylab = "")
        }
        
        # quadratic coefficient (linear component)
        {
            c.hist(qq$qqc.b, xlim = c(-100, 100), d = 1, main = "Black quadratic coefficient", xlab = "", ylab = "")
            c.hist(qq$qqc.g, xlim = c(-100, 100), d = 1, main = "Grey quadratic coefficient", xlab = "", ylab = "")
            c.hist(qq$qqc.w, xlim = c(-100, 100), d = 1, main = "White quadratic coefficient", xlab = "", ylab = "")
        }
        
        dev.off()
    }
    
}

# let's try hclust over those 5 proposed xteristics, in the black channel only...
qq.mini <- qq[,c("mean.b", "min.b", "max.b", "lc.b", "qqc.b")]

qq.nmini <- sapply(qq.mini, function(x) (x - mean(x)) / sd(x))

qq.dist <- dist(qq.nmini)
qq.hc <- hclust(qq.dist)

plot(qq.hc)
zz.2 <- cutree(qq.hc, k = 2)
table(zz.2)

dd.mini.2 <- dd.all[zz.2 == 1,,"black", "160430"]

pdf(paste0(fpath, "tmp.pdf")); {
    par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))
    for (i in 1:length(which(zz.2 == 1))) {
        o.plot(dd.mini.2[i,])
    }
    dev.off()
}


# propose some thresholds, use the 147 persistent noisy px as 'true positives', assess performance

####################################################################################################

# CTE TRAILS                                                                                    ####

o.plot(d[1:30,588,1,"black"], ylim = c(0,10000))
lines(d[1:30,589,1,"black"], col = adjustcolor("blue", alpha = 0.3))

o.plot(d[1:30,588,1,"black"] - d[1:30,587,1,"black"], ylim = c(-1000, 1000))
lines(d[1:30,588,1,"black"] - d[1:30,589,1,"black"], col = adjustcolor("blue", alpha = 0.3))

abline(h = 0, col = "red")

# magnitude of CTE trails in Hubble images? What are we looking for? Suspect more than this....
