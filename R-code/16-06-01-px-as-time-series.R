
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
