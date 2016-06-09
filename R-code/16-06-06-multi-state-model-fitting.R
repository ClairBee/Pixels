
library("IO.Pixels"); library("CB.Misc")

fpath <- "./Drafts/State-space/fig/"
bpm.fpath <- "./Drafts/State-space/bad-px-maps/"

load.pixel.means()

md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

Cat <- c("no response", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "s.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim", "s.dim")

####################################################################################################

# PLOT THRESHOLD SIGMA VS %                                                                     ####

# sigma = SD of per-array mean value, not pixelwise SD
sig <- apply(pw.m, 3:4, sd)

th.prop <- function(im, n.sigma, midpoint = "median", table = F) {
        if (midpoint == "median") {
            m <- median(im)
        } else {
            m <- mean(im)
        }
        
        s <- sd(im)
        if (length(n.sigma) > 1) {
            p <- unlist(lapply(n.sigma, function(x) 100 * sum(im > m + x * s) / length(im)))
        } else {
            p <- 100 * sum(im > m + n.sigma * s) / length(im)
        }
        
        if (table) {
            return(data.frame(n.sigma, 
                              prop = p,
                              GV = m + n.sigma * s))
        } else {
            return(p)
        }
    }

    ccols <- c("black", "blue", "purple", "magenta3", "orange", "gold", "green", "green3", "cyan3", "skyblue", "slateblue1", "orchid")
    ns <- c(2:10)
    
    pdfheight <- 5; pdfwidth = 7
    
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-black.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"black", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"black",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-black.pdf"))
    }
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-grey.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"grey", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"grey",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-grey.pdf"))
    }
    pdf(paste0(fpath, "bad-px-vs-sigma-threshold-white.pdf"), width = pdfwidth, height = pdfheight); {
        plot(ns, th.prop(pw.m[,,"white", 12], n.sigma = ns), type = "o", pch = 20, 
             xlab = expression(paste("Thresholded at median value + n * ", sigma)), 
             ylab = "% of pixels exceeding threshold")
        for(i in 11:2) {
            points(ns, th.prop(pw.m[,,"white",i], n.sigma = ns), type = "o", pch = 20, col = rev(ccols)[i])
        }
        legend("topright", col = ccols, legend = rev(sapply(dimnames(pw.m)[[4]], fancy.date)), pch = 20, bty = "n")
        dev.off()
        crop.pdf(paste0(fpath, "bad-px-vs-sigma-threshold-white.pdf"))
    }

# get % at each setting by sigma
zz <- sapply(dimnames(pw.m)[[4]], function(i) th.prop(pw.m[,,"black", i], n.sigma = ns))
rownames(zz) <- ns
    
matplot(t(zz), type = "o", col = rev(ccols), lty = 1, pch = 20)

# rescale all rows to have same variance & centre. Very similar shape of increase at all levels.
matplot(scale(t(zz), center = T, scale = T), type = "l", col = rev(ccols), lty = 1,
        xaxt = "n")
axis(1, at = c(1:12), labels = sapply(dimnames(pw.m)[[4]], fancy.date), srt = 90)    

####################################################################################################

# CHECK DEVELOPMENT OF PER-PANEL GRADIENT OVER TIME                                             ####

# fit linear model to each image in turn

panel.lm.black <- lapply(apply(pw.m[,,"black",], 3, panel.lm, robust = T), "[", "models")

x.grad <- do.call("rbind", lapply(lapply(panel.lm.black, "[[", "models"), "[",,"x"))
y.grad <- do.call("rbind", lapply(lapply(panel.lm.black, "[[", "models"), "[",, "y"))

colnames(x.grad) <- colnames(y.grad) <- apply(cbind(c(rep("U", 16), rep("L", 16)), rep(c(1:16), 2)), 1, paste, collapse = "")

par(mfrow = c(2, 1))
matplot(x.grad, type = "l")
matplot(y.grad, type = "l")
par(mfrow = c(1,1))

ccols <- c("black", "blue", "purple", "magenta3", "orange", "gold", "green", "green3", "cyan3", "skyblue", "slateblue1", "orchid")

plot(x.grad[,"U1"], y.grad[,"U1"], type = "o", pch = 21, bg = ccols, xlim = c(0,5), ylim = c(0,1.5),
     xlab = "x-gradient fitted", ylab = "y-gradient fitted")
points(x.grad[,"U2"], y.grad[,"U2"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U3"], y.grad[,"U3"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U4"], y.grad[,"U4"], type = "o", pch = 21, bg = ccols)
points(x.grad[,"U5"], y.grad[,"U5"], type = "o", pch = 21, bg = ccols)

points(x.grad[,"U14"], y.grad[,"U14"], type = "o", pch = 21, bg = ccols)

plot(x.grad[1,1:16], type = "l", ylim = c(-7,5))
lines(x.grad[1,17:32])
abline(h = 0, lty = 2)
for (i in 2:12) {
    lines(x.grad[i,1:16], col = rev(ccols)[i])
    lines(x.grad[i,17:32], col = rev(ccols)[i])
}

plot(y.grad[1,1:16], type = "l", ylim = c(-1.5,1))
lines(y.grad[1,17:32])
abline(h = 0, lty = 2)
for (i in 2:12) {
    lines(y.grad[i,1:16], col = rev(ccols)[i])
    lines(y.grad[i,17:32], col = rev(ccols)[i])
}


####################################################################################################

# THRESHOLD BEST CASE VS WORST CASE                                                             ####

####################################################################################################

# CLASSIFY ALL PIXELS                                                                           ####

# standard approach uses white, grey and black images for thresholding
bp.std <- readRDS(paste0(bpm.fpath, "bpx-initial.rds"))

# standard thresholds, but without using white image for thresholding (only for no response/)
{
    bp <- lapply(dimnames(pw.m)[[4]],
                 function(x) rbind(data.frame(edge.px(pw.m), type = ordered("edge", levels = Cat)),
                                   data.frame(no.response(x), type = "no response"),
                                   data.frame(which(pw.m[, , "black", x] == 65535, arr.ind = T), type = "hot"),
                                   data.frame(which(pw.m[, , "grey", x] == 0, arr.ind = T), type = "dead"),
                                   screen.spots.xy(x),
                                   get.dim.bright.px(pw.m[,,"grey", x]),
                                   get.dim.bright.px(pw.m[,,"black", x]),
                                   data.frame(which(find.lines(pw.m[, , "black", x]) > 0, arr.ind = T), type = "line.b"),
                                   data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * -2) == 0, arr.ind = T), type = "l.dim"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * 2) > 0, arr.ind = T), type = "l.bright"),
                                   data.frame(which(threshold(md.g[[x]], level = mad(pw.m[,,"grey", x]) * -2) == 0, arr.ind = T), type = "l.dim")))
    
    bp <- lapply(lapply(bp, 
                        function(x) x[order(x$type),]),
                 function(x) x[!duplicated(x[,1:2]),])
    names(bp) <- dimnames(pw.m)[[4]]
    
    saveRDS(bp, paste0(bpm.fpath, "bpx-standard-excl-white.rds"))
}
bp.std.bg <- readRDS(paste0(bpm.fpath, "bpx-standard-excl-white.rds"))

####################################################################################################

# FIT MARKOV MODEL IN R                                                                         ####

# need a dataframe with time of observation, observed state, and subject id.
# observations must be grouped by subject, and ordered by time within each subject.

# convert date strings to actual dates
dt <- as.Date(names(bp), format = "%y%m%d")

# convert pixel map of types to data frame of all pixels, with pixel id.
bpx.im <- abind(lapply(bp, bpx2im), along = 3)

df.all <- data.frame(id = sort(rep(1:length(bpx.im[,,1]), 12)),
                     acq = rep(as.Date(names(bp), "%y%m%d"), length(bpx.im[,,1])),
                     state = c(t(apply(bpx.im, 3, c))))


# also some data cleaning needed; use switcher table to set correct values
{
    # edges are set to type 99: censored state (state not observed)
    Cat.conv <- data.frame(org = c(0:15), 
                           new = c(1, 2, 3, 4, 4, 5, 6, 5, 1, 1, 1, 99, 3, 3, 3, 1))
    
    # check conversion
    cbind(c("ok", Cat)[Cat.conv$org + 1], Cat.conv$new)
    
    # make conversion (only takes few seconds)
    df$state2 <- Cat.conv$new[match(df$state, Cat.conv$org)]
}

# create an MSM

st <- statetable.msm(state, id, data = df.all)       # 75s elapsed
colnames(st) <- c("ok", Cat)[as.numeric(colnames(st)) + 1]
rownames(st) <- c("ok", Cat)[as.numeric(rownames(st)) + 1]

write.csv(st, paste0(bpm.fpath, "states-all-standard.csv"), quote = F)

Q <- rbind(c(0, 0, 0.001, 0.01, 0.01, 0.01),
           c(0, 0, 0, 0, 0, 0.01),
           c(0, 0.1, 0, 0.04, 0, 0),
           c(0.15, 0, 0.01, 0, 0.01, 0.01),
           c(0, 0.01, 0.01, 0.01, 0, 0.01),
           c(0.15, 0.01, 0.01, 0.01, 0.01, 0))

system.time(msm.1 <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, censor = 99))

####################################################################################################

# FIT MSM TO SINGLE SUBPANEL                                                                    ####

bp <- readRDS(paste0(bpm.fpath, "bpx-initial.rds"))
bpx.im <- abind(lapply(bp, bpx2im), along = 3)

sp.u4 <- bpx.im[383:510,993:1996,]
df <- data.frame(id = sort(rep(1:length(sp.u4[,,1]), 12)),
                 acq = rep(as.Date(names(bp), "%y%m%d") - as.Date(names(bp)[1], "%y%m%d"), length(sp.u4[,,1])),
                 state = c(t(apply(sp.u4, 3, c))))

# simplify states
Cat.conv <- data.frame(org = c(0:15), 
                       new = c(1,7,6,4,4,3,8,2,1,1,1,99,6,6,5,1))
ns <- c("ok", "l.bright", "bright", "hot", "l.dim", "dim", "no.resp", "line")

# check conversion
cbind(c("ok", Cat)[Cat.conv$org + 1], ns[Cat.conv$new])

df$state2 <- Cat.conv$new[match(df$state, Cat.conv$org)]

st <- statetable.msm(state2, id, data = df)
rownames(st) <- ns[as.numeric(rownames(st))]
colnames(st) <- ns[as.numeric(colnames(st))]
st
# edge pixels are treated as censored

Q <- rbind (c(1,1,0,0,1,0,0,1),      # can move from healthy to locally bright, locally dim, or lines
            c(1,1,1,0,0,0,0,1),      # locally bright can recover or become bright (or part of a line)
            c(0,1,1,1,0,0,0,1),      # bright > locally bright, or hot
            c(0,0,1,1,0,0,0,1),      # 'hot' may be misclassified as very bright?
            c(1,0,0,0,1,1,0,0),      # locally dim can recover or become dim (or part of a line)
            c(0,0,0,0,1,1,1,0),      # dim > locally dim or non-responsive? (or part of a line)
            c(0,0,0,0,0,1,1,0),      # no response > absorbing (in this model)
            c(0,0,0,0,0,0,0,1))      # lines will remain as such (until classification changes)
            
system.time(msm <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, censor = 99))
# obstype: default is 1 (all states observed at arbitrary times)

# still aborts after 75 seconds or so...

# let's try a single column. This may be possible.
{
    cc <- 30
    df <- data.frame(id = sort(rep(1:length(sp.u4[cc,,1]), 12)),
                     acq = rep(as.Date(names(bp), "%y%m%d") - as.Date(names(bp)[1], "%y%m%d"), length(sp.u4[cc,,1])),
                     state = c(t(apply(sp.u4[cc,,], 2, c))))
    
    Cat.conv <- data.frame(org = c(0:15), 
                           new = c(1,1,1,1,1,1,1,2,3,4,5,5,7,8,9,10))
    df$state2 <- Cat.conv$new[match(df$state, Cat.conv$org)]
    
    QQ <- rbind(c(1,0.001, 0.001,0.001,0),
                c(0.001,1,0,0,0),
                c(0.001,0,1,0.001,0),
                c(0.001,0,0,1,0.001),
                c(0.001,0,0,0,1))
    
    QQ <- sweep(st, 1, rowSums(st), "/")
    system.time(msm.col <- msm(state2 ~ acq, subject = id, data = df, qmatrix = QQ))
}
# didn't converge. ARGH

# try splitting subpanel in 2 (and then 2 again & again to obtain squares, if necessary)...
{
    df <- data.frame(id = sort(rep(1:length(sp.u4[,1:512,1]), 12)),
                     acq = rep(as.Date(names(bp), "%y%m%d") - as.Date(names(bp)[1], "%y%m%d"), length(sp.u4[,1:512,1])) / 365,
                     state = c(t(apply(sp.u4[,1:512,], 3, c))))
    
    # simplify states
    Cat.conv <- data.frame(org = c(0:15), 
                           new = c(1,7,6,4,4,3,8,2,1,1,1,99,6,6,5,1))
    ns <- c("ok", "l.bright", "bright", "hot", "l.dim", "dim", "no.resp", "line")
    df$state2 <- Cat.conv$new[match(df$state, Cat.conv$org)]
    table(df$state2)

    st <- statetable.msm(state2, id, data = df)
    rownames(st) <- ns[as.numeric(rownames(st))]
    colnames(st) <- ns[as.numeric(colnames(st))]
    st
    round(sweep(st, 1, rowSums(st), "/"), 3)
    
    Q <- rbind (c(1,0.001,0,0,0.0001,0,0,0.001),      # can move from healthy to locally bright, locally dim, or lines
                c(0.154,1,0.013,0,0,0,0,0.001),       # locally bright can recover or become bright (or part of a line)
                c(0,0.09,1,0.038,0,0,0,.019),         # bright > locally bright, or hot
                c(0,0,.042,1,0,0,0,0.0001),           # 'hot' may be misclassified as very bright?
                c(0.25,0,0,0,1,.0001,0,0),            # locally dim can recover or become dim (or part of a line)
                c(0,0,0,0,0.01,1,0.0001,0),           # dim > locally dim or non-responsive? (or part of a line)
                c(0,0,0,0,0,0.01,1,0.0001),           # no response > absorbing (in this model)
                c(0.0001,0,0,0,0,0,0,1))              # lines will remain as such (until classification changes)
    
    system.time(msm.sp2 <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, 
                               control = list(fnscale = 5000)))

    # numerical overflow in calculating likelihood -> add 'control = list(fnscale = 5000)'
    # BUT failed to converge before reaching iteration limit (after ~4mins)
    
    system.time(msm.sp2 <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, 
                               method = "SANN",                     # simulated annealing in optim instead
                               control = list(fnscale = 5000)))
    
    
    
    system.time(msm.sp2 <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, 
                               opt.method = "fisher",
                               control = list(fnscale = 5000, damp = 1)))
    # back to the overflow problem, despite scaling
    # added damp = 1: system is computationally singular/
    
    system.time(msm.sp2 <- msm(state2 ~ acq, subject = id, data = df, qmatrix = Q, 
                               opt.method = "nlm", fscale = 5000))
    # numerical overflow again
    
}

# view msm obtained
{
    statetable.msm(state2, id, data = df)
    
    qm <- prep.csv(qmatrix.msm(msm.sp2, ci = "none"), dp = 3)
    colnames(qm) <- rownames(qm) <- ns    
    write.csv(qm, paste0(bpm.fpath, "qm-2.csv")) 
    
    pm <- prep.csv(pmatrix.msm(msm.sp2), dp = 3)
    colnames(pm) <- rownames(pm) <- ns
    write.csv(pm, paste0(bpm.fpath, "pm-2.csv")) 
    
}
####################################################################################################

# COMPARE MSM FOR DIFFERENT BAD PIXEL MAPS                                                      ####

####################################################################################################

# BRIGHT VS HOT PIXELS                                                                          ####

# hot pixels: those that produce a high current all the time (even in black images)

hot <- get.dim.bright.px(pw.m[,,"black", "160430"])

# bright pixels: those that produce a high current in response to spot
bright <- get.dim.bright.px(pw.m[,,"grey", "160430"])

hb <- table(black = bpx2im(hot), grey = bpx2im(bright))
colnames(hb) <-  c("ok", "v.bright", "bright", "s.bright", "v.dim", "dim", "s.dim")[as.numeric(colnames(hb)) + 1]
rownames(hb) <-  c("ok", "v.bright", "bright", "s.bright", "v.dim", "dim", "s.dim")[as.numeric(rownames(hb)) + 1]
hb
