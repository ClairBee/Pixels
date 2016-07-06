
# ALTERNATIVE THRESHOLDING APPROACHES

library("IO.Pixels"); library("CB.Misc")

fpath <- "./Notes/Regression-thresholding/fig/"

acq.160430 <- readRDS("./02_Objects/images/pwm-160430.rds")
acq.131122 <- readRDS("./02_Objects/old-data/pwm-131122.rds")
acq.140128 <- readRDS("./02_Objects/old-data/pwm-140128.rds")

.smoothScatter <- hijack(smoothScatter, nrpoints = 0,
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))

md.160430 <- readRDS("./02_Objects/med-diffs/md-160430.rds")

####################################################################################################

# FUNCTIONS TO GET VARIOUS BAD PIXEL LISTS                                                      ####

# basic thresholding (incl. locally bright/dim px)
basic.bpx <- function(dt) {
    dt <- toString(dt)
    acq <- eval(parse(text = paste0("acq.", dt)))
    
    # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    # if not already present, get median differences
    if (!exists(eval(parse(text = paste0("md.", dt))))) {
        md <- readRDS(paste0("./02_Objects/med-diffs/md-", dt, ".rds"))
    } else {
        md <- eval(parse(text = paste0("md.", dt)))
    }
    
    px <- rbind(data.frame(edge.px(acq, edge.width = 40), type = "edge"),
                    data.frame(no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]), type = "no.resp"),
                    data.frame(which(acq[,,"black"] == 65535, arr.ind = T), type = "hot"),
                    data.frame(which(acq[,,"grey"] == 0, arr.ind = T), type = "dead"),
                    data.frame(screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40), type = "screen.spot"),
                    data.frame(which(find.lines(acq[, , "black"], midline = 1024.5) > 0, arr.ind = T), type = "line.b"),
                    data.frame(which(acq[, , "black"] > th["v.bright", "black"], arr.ind = T), type = "v.bright"),
                    data.frame(which(acq[, , "black"] > th["bright", "black"], arr.ind = T), type = "bright"),
                    data.frame(which(acq[, , "black"] < th["v.dim", "black"], arr.ind = T), type = "v.dim"),
                    data.frame(which(acq[, , "black"] < th["dim", "black"], arr.ind = T), type = "dim"),
                    data.frame(which(acq[, , "grey"] > th["v.bright", "grey"], arr.ind = T), type = "v.bright"),
                    data.frame(which(acq[, , "grey"] > th["bright", "grey"], arr.ind = T), type = "bright"),
                    data.frame(which(acq[, , "grey"] < th["v.dim", "grey"], arr.ind = T), type = "v.dim"),
                    data.frame(which(acq[, , "grey"] < th["dim", "grey"], arr.ind = T), type = "dim"),
                    data.frame(which(md[, , "black"] > 1200, arr.ind = T), type = "l.bright"),
                    data.frame(which(md[, , "grey"] > 1200, arr.ind = T), type = "l.bright"),
                    data.frame(which(md[, , "black"] < -1200, arr.ind = T), type = "l.dim"),
                    data.frame(which(md[, , "grey"] < -1200, arr.ind = T), type = "l.dim"))
    
    Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "edge", "l.bright", "line.d", "screen.spot", "v.dim", "dim", "l.dim")
    px$type <- ordered(px$type, levels = Cat)
    
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    return(px)
}

# qualitative thresholding
qual.bpx <- function(dt, w.th = 1000) {
    dt <- toString(dt)
    acq <- eval(parse(text = paste0("acq.", dt)))
    
    # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    wlm <- fit.w.lm(acq)$df
    
    px <- rbind(data.frame(which(acq[,,"black"] > th["bright", "black"], arr.ind = T), type = "warm"),
                data.frame(which(acq[,,"grey"] > th["bright", "grey"], arr.ind = T), type = "bright"),
                data.frame(setNames(wlm[abs(wlm$res) > w.th, 1:2], nm = c("row", "col")), type = "nonlinear"),
                data.frame(which(acq[,,"grey"] - acq[,,"black"] < 500 & acq[,,"grey"] > acq[,,"black"], arr.ind = T), type = "stuck"),
                data.frame(which(acq[,,"black"] < th["dim", "black"], arr.ind = T), type = "cool"),
                data.frame(which(acq[,,"grey"] < th["dim", "grey"], arr.ind = T), type = "dim"))
    
    px$type <- ordered(px$type, levels = c("stuck", "warm", "bright", "nonlinear", "cool", "dim"))
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    
    return(px)
}


qq <- qual.bpx(160430)
hh <- basic.bpx(160430)

plot(qq[,1:2], pch = 20)

####################################################################################################

# ILLUSTRATE LINEAR RELATIONSHIP EVEN IN MOST DAMAGED IMAGES                                    ####

# FUNCTIONS

fit.w.lm <- function(im) {
    
    # fit linear model excluding edge pixels
    df <- setNames(data.frame(melt(im[,,"black"]), 
                              melt(im[,,"grey"]),
                              melt(im[,,"white"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    w.lm <- lm(w ~ b * g,
                  data = df[findInterval(df$x, c(40.5, 2008.5)) == 1 &
                                findInterval(df$y, c(40.5, 2008.5)) == 1,])
    
    df$fv <- predict(w.lm, df[,c("b", "g")])
    df$res <- df$w - df$fv
    
    list(df = df, r2 = round(summary(w.lm)$adj.r.squared, 3), rmse = round(summary(w.lm)$sigma, 2))
}

wlm.plot <- function(dt) {
    
    dt <- toString(dt)
    
    wlm <- eval(parse(text = paste0("lm.", dt)))
    
    pdf(paste0(fpath, "white-fit-", dt, ".pdf"))
    par(mar = c(4, 4, 1, 1))
    smoothScatter(wlm$df$w, wlm$df$fv, nrpoints = 0, xlim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed", ylab = "Fitted value")
    abline(line(wlm$df$w, wlm$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
    
    write(paste0("Adj $r^2$ ", round(wlm$r2, 3), "; ",
                 "RMSE ", round(wlm$rmse, 2)),
          paste0(fpath, "fitted-wv-all-", dt, ".txt"))
}

#--------------------------------------------------------------------------

lm.160430 <- fit.w.lm(acq.160430)
lm.131122 <- fit.w.lm(acq.131122)
lm.140128 <- fit.w.lm(acq.140128)

wlm.plot(131122)
wlm.plot(140128)
wlm.plot(160430)

# remove RMSE & R^2, now only need data frame
lm.131122 <- lm.131122$df; lm.140128 <- lm.140128$df; lm.160430 <- lm.160430$df;

####################################################################################################

# REMOVE NON-LINEAR PX & PLOT SHADING CORRECTION                                                ####

lm.131122$type[abs(lm.131122$res) > 1000] <- "nl"
lm.140128$type[abs(lm.140128$res) > 1000] <- "nl"
lm.160430$type[abs(lm.160430$res) > 1000] <- "nl"

sc.131122 <- shading.corrected(acq.131122)
sc.140128 <- shading.corrected(acq.140128)
sc.160430 <- shading.corrected(acq.160430)

pdf(paste0(fpath, "white-vs-sc-131122.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.131122[,,"white"], sc.131122, ylim = range(pretty(sc.131122, na.rm = T)),
                   xlab = "White value", ylab = "Shading-corrected value", xlim = c(0,65535))
    dev.off()
}
pdf(paste0(fpath, "white-vs-sc-140128.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.140128[,,"white"], sc.140128, ylim = c(-10000,max(pretty(sc.140128, na.rm = T))),
                   xlab = "White value", ylab = "Shading-corrected value", xlim = c(0,65535))
    dev.off()
}
pdf(paste0(fpath, "white-vs-sc-160430.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.160430[,,"white"], sc.160430, ylim = range(pretty(sc.160430, na.rm = T)),
                   xlab = "White value", ylab = "Shading-corrected value", xlim = c(0,65535))
    dev.off()
}


pdf(paste0(fpath, "white-vs-sc-without-nonlinear-131122.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.131122[,,"white"][as.matrix(lm.131122[is.na(lm.131122$type), 1:2])],
                   sc.131122[as.matrix(lm.131122[is.na(lm.131122$type), 1:2])],
                   ylim = range(pretty(sc.131122, na.rm = T)), xlim = c(0,65535),
                   xlab = "White value", ylab = "Shading-corrected value")
    dev.off()
}
write(sum(!is.na(lm.131122$type)), paste0(fpath, "nl-px-131122.txt"))

pdf(paste0(fpath, "white-vs-sc-without-nonlinear-140128.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.140128[,,"white"][as.matrix(lm.140128[is.na(lm.140128$type), 1:2])],
                   sc.140128[as.matrix(lm.140128[is.na(lm.140128$type), 1:2])],
                   ylim = range(pretty(sc.140128, na.rm = T)), xlim = c(0,65535),
                   xlab = "White value", ylab = "Shading-corrected value")
    dev.off()
}
write(sum(!is.na(lm.140128$type)), paste0(fpath, "nl-px-140128.txt"))

pdf(paste0(fpath, "white-vs-sc-without-nonlinear-160430.pdf")); {
    par(mar = c(4,4,1,1))
    .smoothScatter(acq.160430[,,"white"][as.matrix(lm.160430[is.na(lm.160430$type), 1:2])],
                   sc.160430[as.matrix(lm.160430[is.na(lm.160430$type), 1:2])],
                   ylim = range(pretty(sc.160430, na.rm = T)), xlim = c(0,65535),
                   xlab = "White value", ylab = "Shading-corrected value")
    dev.off()
}
write(sum(!is.na(lm.160430$type)), paste0(fpath, "nl-px-160430.txt"))

####################################################################################################

# REMOVE STUCK PIXELS                                                                           ####



# NONLINEARITY                                                                                  ####

.smoothScatter(acq.131122[,,"black"][which(acq.131122[,,"white"] < 10000, arr.ind = T)],
               acq.131122[,,"grey"][which(acq.131122[,,"white"] < 10000, arr.ind = T)])
summary(c(acq.131122[,,"grey"][which(acq.131122[,,"white"] < 10000, arr.ind = T)] - 
              acq.131122[,,"black"][which(acq.131122[,,"white"] < 10000, arr.ind = T)]))
abline(0,1)

abline(500,1, lty = 2); abline(-500,1,lty = 2)

# produce data frame of all necessary variables
{
    # data frame of all variables for active region of image
    df <- setNames(data.frame(melt(acq[,,"black"]), 
                              melt(acq[,,"grey"]),
                              melt(acq[,,"white"]),
                              melt(sc))[,c("X1", "X2", "value", "value.1", "value.2", "value.3")],
                   nm = c("x", "y", "b", "g", "w", "sc"))
    df$type[findInterval(df$x, c(40.5, 2008.5)) != 1 | findInterval(df$y, c(40.5, 2008.5)) != 1] <- "edge"
    
    # fit linear model, excluding edge pixels & NA padding
    w.lm <- rlm(w ~ b * g, data = df[!is.na(df$b) & is.na(df$type),])
    
    df$fv <- predict(w.lm, df[,c("b", "g")])
    df$res <- df$w - df$fv
}

df$type[abs(df$res) > 1000 & is.na(df$type)] <- "nonlinear"
df$type[df$b > 3 * (median(acq[,,"black"], na.rm = T) + 21845)/4] <- "bright"

smoothScatter(acq[,,"black"][as.matrix(df[is.na(df$type),1:2])], 
              sc[as.matrix(df[is.na(df$type), 1:2])], nrpoints = 0, 
              colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"))

plot(df[df$type %in% c("bright", "nonlinear"), 1:2], pch = 20)

sc.lin <- sc
sc.lin[as.matrix(df[abs(df$res) > 1000, 1:2])] <- NA
3 * (median(acq[,,"black"], na.rm = T) + 21845)/4

hist(sc.lin[!is.na(sc.lin) & !is.infinite(sc.lin)], breaks = "fd", ylim = c(0,5))

# only 5px out of spec
which(sc.lin > 20000, arr.ind = T)
acq[1751,467,]
acq[1670,992,]

df[which(df$b > 20000),]

df[df$x == 1751 & df$y == 467,]
)

####################################################################################################

# STUCK PIXELS                                                                                  ####