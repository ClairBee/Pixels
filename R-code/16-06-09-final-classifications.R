
# COMPARISON OF PROPOSED CLASSIFICATION METHODS & SELECTION OF FINAL APPROACH

# considering observable state spaces only

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Final-classifications/fig/"

load.pixel.means()
load.pixel.sds()

# median-differenced images to find local non-uniformity
md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

####################################################################################################

# INCLUDE A NOTE ON LOOKING FOR DIM & HORIZONTAL LINES, THEN FILTERING MANUALLY
    # dim line segments may appear if there are two v bright points adjacent to it within a short distance of one another
    # use 
# WHAT IS TARGET VALUE OF AIR, AND HOW IS THIS RELATED TO GREY MEAN VALUE? IS GREY REPRESENTATIVE?

####################################################################################################

# COMMON ELEMENTS OF BAD PIXEL MAP                                                              ####

# bad pixel types that don't change with thresholds: 
    # edges, non-responsives, abs. hot/dead, screen spots
{
    bp.const <- apply(pw.m, 4,
                      function(acq) rbind(data.frame(edge.px(acq), type = "edge"),
                                          data.frame(no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]), type = "no.resp"),
                                          data.frame(which(acq[,,"black"] == 65535, arr.ind = T), type = "hot"),
                                          data.frame(which(acq[,,"white"] == 0, arr.ind = T), type = "dead")))
    saveRDS(bp.const, paste0(fpath, "bad-px-constant.rds"))
    
    bp.screen <- apply(pw.m, 4,
                       function(acq) screen.spots(acq[,,"white"], enlarge = T))

    bp.screen <- lapply(bp.screen, 
                        function(df) if (is.null(df)) {
                            NULL
                        } else {
                            data.frame(df, type = "screen.spot")
                        })
    saveRDS(bp.screen, paste0(fpath, "bad-px-screenspots.rds"))
    }

    
# lines (may want to observe these pixels singly, so keep separate)
{
    bp.lines <- apply(pw.m, 4,
                      function(acq) data.frame(which(find.lines(acq[, , "black"]) > 0, arr.ind = T), type = "line.b"))
    saveRDS(bp.lines, paste0(fpath, "bad-px-bright-lines.rds"))
    
    # no useful dim columns, or row artefacts of either kind.
        qq <- apply(pw.m, 4, function(acq) which(find.lines(acq[, , "black"], dim.lines = T) > 0, arr.ind = T))
        qq <- apply(pw.m, 4, function(acq) which(find.lines(acq[, , "black"], horizontal = T) > 0, arr.ind = T))
        qq <- apply(pw.m, 4, function(acq) which(find.lines(acq[, , "black"], dim.lines = T, horizontal = T) > 0, arr.ind = T))
}

# locally dim & bright pixels
{
    bp.local <- lapply(names(md.b),
                       function(dt) rbind(data.frame(which(md.b[[dt]] > 2 * mad(pw.m[,,"black", dt]), arr.ind = T),
                                                     type = "l.bright"),
                                          data.frame(which(md.b[[dt]] < -2 * mad(pw.m[,,"black", dt]), arr.ind = T),
                                                     type = "l.dim"),
                                          data.frame(which(md.g[[dt]] > 2 * mad(pw.m[,,"grey", dt]), arr.ind = T),
                                                     type = "l.bright"),
                                          data.frame(which(md.g[[dt]] < -2 * mad(pw.m[,,"grey", dt]), arr.ind = T),
                                                     type = "l.dim")))
    saveRDS(bp.local, paste0(fpath, "bad-px-local-bright-dim-pixels.rds"))
}

####################################################################################################

# FUNCTIONS                                                                                     ####

# combine & rationalise partial bad pixel maps
merge.bp.lists <- function(bp.list, cat.order) {
    
    bpx <- list()
    for (i in 1:length(bp.list[[1]])) {
        bpx[[i]] <- rbind.fill(lapply(bp.list, "[[", i))
        bpx[[i]]$type = ordered(bpx[[i]]$type, levels = cat.order)
        bpx[[i]] <- bpx[[i]][order(bpx[[i]]$type),]
        bpx[[i]] <- bpx[[i]][!duplicated(bpx[[i]][,1:2]),]
    }
    names(bpx) <- names(bp.list[[which.max(lapply(lapply(bp.list, names), length))]])
    return(bpx)
}

# plot transition proportion at each acquisition
tr.prop.plot <- function(tr.mat, times, from, to, add = F, ...) {
    
    if (add) {
        lines(times[2:length(times)] - times[1], 
              tr.mat[from, to,], type = "o", pch = 20, ...)
    } else {
        plot(times[2:length(times)] - times[1],
             tr.mat[from, to,], type = "o", pch = 20, 
             ylab = "Proportion transitioning" , xlab = "Time after first acquisition", ...)
    }
}

# get gradient of line fitted to transition proportion at each acquisition
tr.prop.slope <- function(tr.mat, times, from, to) {
    
    line(times[2:length(times)] - times[1],
         tr.mat[from, to,])$coef[2]
}

####################################################################################################

# MEAN-VALUED CLASSIFICATION                                                                    ####
    # this is the 'standard' classification approach used throughout the project so far.
    # simplest approach - use as baseline against which to look for improvements.

# define possible bad pixel states 
Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "edge", "l.bright", "line.d", "screen.spot", "v.dim", "dim", "l.dim")
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "yellow", "grey", "violet", NA, "green3", "lightskyblue", "grey")


# classify bad pixels by global thresholding
{
    bp.th <- apply(pw.m, 4, 
                   function(acq) rbind(data.frame(unique(rbind.fill(apply(acq, 3, function(img) data.frame(which(img > median(img) + (65535 - median(img)) / 2, arr.ind = T))))),
                                                  type = "v.bright"),
                                       data.frame(unique(rbind.fill(apply(acq, 3, function(img) data.frame(which(img > median(img) + (65535 - median(img)) / 4, arr.ind = T))))),
                                                  type = "bright"),
                                       data.frame(unique(rbind.fill(apply(acq, 3, function(img) data.frame(which(img < median(img) / 2, arr.ind = T))))),
                                                  type = "v.dim"),
                                       data.frame(unique(rbind.fill(apply(acq, 3, function(img) data.frame(which(img < median(img) * 0.75, arr.ind = T))))),
                                                  type = "dim")))
    saveRDS(bp.th, paste0(fpath, "bad-px-thresholded-mean-values-incl-white.rds"))
}
bp <- merge.bp.lists(list(readRDS(paste0(fpath, "bad-px-thresholded-mean-values-incl-white.rds")), 
                          readRDS(paste0(fpath, "bad-px-constant.rds")),
                          readRDS(paste0(fpath, "bad-px-screenspots.rds")),
                          readRDS(paste0(fpath, "bad-px-bright-lines.rds")),
                          readRDS(paste0(fpath, "bad-px-local-bright-dim-pixels.rds"))), cat.order = Cat)

# any pixels marked as 'screen spot' are actually okay (after duplicates are removed), so filter these out
# classification is needed intitially to discriminate between genuinely dim pixels and shadows from spot
bp <- lapply(bp, function(bpx) bpx[bpx$type != "screen.spot",])

# also not interested in edge pixels, so remove these as well
bp <- lapply(bp, function(bpx) bpx[bpx$type != "edge",])


# get transition matrix
tr <- abind(get.transitions(bp), along = 3)

# filter out categories that don't occur at all
tr <- tr[colSums(apply(tr, 1:2, sum)) * rowSums(apply(tr, 1:2, sum)) > 0,
         colSums(apply(tr, 1:2, sum)) * rowSums(apply(tr, 1:2, sum)) > 0,]

# plot transitions as function of real time - are rates constant?
{
    # get matrix of all transitions, to see which actually occur
    apply(tr, 1:2, sum)
    
    tr.prop <- array(apply(tr, 3, function(x) sweep(x, 1, rowSums(x), "/")), dim = dim(tr), dimnames = dimnames(tr))
    
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "hot", to = "hot", ylim = c(0,1),
                 main = "Transitions from hot state at each acquisition")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "hot", to = "v.bright", add = T, col = "red")
    
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "v.bright", to = "v.bright", ylim = c(0,1),
                 main = "Transitions from very bright state at each acquisition")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "v.bright", to = "hot", add = T, col = "red")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "v.bright", to = "bright", add = T, col = "blue")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "v.bright", to = "l.bright", add = T, col = "green3")
    
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "l.bright", to = "l.bright", ylim = c(0,1),
                 main = "Transitions from locally bright state at each acquisition")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "l.bright", to = "normal", add = T, col = "gold")
    tr.prop.plot(tr.prop, as.Date(names(bp), "%y%m%d"), from = "l.bright", to = "bright", add = T, col = "red")
    

    tr.prop.slope(tr.prop, as.Date(names(bp), "%y%m%d"), from = "hot", to = "v.bright")
}
# no particular sign of rates worsening over time. Proportions remain reasonably constant.

# check behaviour of each class under flat field correction
# possible to make some sort of 3d plot in b/w/g 

# check stability of classes
{
    tr.prop <- array(apply(tr, 3, function(x) 100 * x / rowSums(x)),
                     dim = dim(tr), dimnames = dimnames(tr))
    
    prep.csv(apply(tr, 1:2, mean), dp = 0)          # mean transition rates (in pixels)
    prep.csv(apply(tr, 1:2, sd), dp = 0)            # mean transition SDs (in pixels)
    
    write.csv(prep.csv(apply(tr.prop, 1:2, mean, na.rm = T)),
              paste0(fpath, "mean-valued-tr-props-mean.csv"), quote = F)

    write.csv(prep.csv(apply(tr.prop, 1:2, sd, na.rm = T)),
              paste0(fpath, "mean-valued-tr-props-sd.csv"), quote = F)
    
    tr.mean <- apply(abind(prep.csv(apply(tr.prop, 1:2, mean), 2), 
                           array(" \\textit{(", dim = dim(tr.prop)[1:2]),
                           prep.csv(apply(tr.prop, 1:2, sd),2), 
                           array(")}", dim = dim(tr.prop)[1:2]),
                           along = 3), 1:2, paste, collapse = "")
    tr.mean[tr.mean == "- \\textit{(-)}"] <- "-"
    write.csv(tr.mean, paste0(fpath, "mean-valued-tr-props.csv"), quote = F)
    diag(tr.mean) <- apply(cbind("\\textbf{", diag(tr.mean), "}"), 1, paste, collapse = "")
    }

# cluster pixels according to path taken?
{
    # get list of all pixels identified as unhealthy: 147132 pixels
    qq <- rmerge.df.list(unique(rbind.fill(lapply(bp, "[", 1:2))), bp, by = c(1:2), all = T)
    
    qq[,sapply(qq, is.factor)] <- lapply(qq[,sapply(qq, is.factor)], as.numeric)
    qq[is.na(qq)] <- 0
    qq$path <- apply(qq[,3:14], 1, paste, collapse = "")
    
    # identify static pixels
    zz <- apply(qq[1:5,3:14], 1, function(r) all(sapply(r, function(x) x == r[1])))
    
    # plot trajectories through state space
    qq$path.type[which(apply(qq[,3:14], 1, function(r) all(r == r[1])))] <- "static"
    table(qq$path[which(qq$path.type == "static")])
    
    # check direction of movement through states
    qq$zz <- apply(qq[,3:14], 1, function(vv) max(as.numeric(gsub(0,99, vv[2:12])) - as.numeric(gsub(0,99,vv[1:11]))))

#********************************************************************************************    
    
    # probably better to convert states so that 0 is at centre
    
    state.conv <- data.frame(c.name = c("normal", Cat), org = c(0:13),
                        new = c(0,-4,-5,4,3,2,10,0,1,-10,0,-3,-2,-1))
    state.conv[order(state.conv[,3]),]

    hh <- apply(qq[,3:14], 2, function(r) state.conv$new[match(r, state.conv$org)])
    
    tt <- hh[,2:12] - hh[,1:11]
    
    tt.min <- apply(tt, 1, min)
    
    qq$path.type <- "unclassified"
    qq$path.type[which(apply(tt, 1, min) == 0 & apply(hh, 1, min) == 0)] <- "increasing"
    qq$path.type[which(apply(tt, 1, max) == 0 & apply(hh, 1, max) == 0)] <- "decreasing"
    qq$path.type[which(rowSums(tt) == 0)] <- "static"
    
    count(qq$path.type)
    
    head(qq[qq$path.type == "decreasing",], 20)
    head(qq[qq$path.type == "increasing",], 20)
    
}

# observed vs expected prevalence, chi-squared test?
{}

####################################################################################################

# MEAN-VALUED CLASSIFICATION (GREY/BLACK ONLY)                                                  ####

####################################################################################################

# QUALITATIVE CLASSIFICATION                                                                    ####

####################################################################################################

# FEATURE-BASED CLASSIFICATION                                                                  ####

####################################################################################################

# SUPPORT FOR CHOICE OF IMAGES USED IN CLASSIFICATION                                           ####

# identify bright lines in all images & compare output
{
    
}

####################################################################################################

# MEAN VS VARIANCE PER SUBPANEL                                                                 #### 

panel.edges()

plot(pw.m[1:126, 1:992, "grey", 1], pw.sd[1:126, 1:992, "grey", 1]^2, pch = 20, 
     xlab = "mean", ylab = "variance", xlim = c(10000, 30000), ylim = c(0,250000))
points(pw.m[1:126, 1:992, "grey", 2], pw.sd[1:126, 1:992, "grey", 2]^2, pch = 20,
       col = adjustcolor("blue", alpha = 0.4))
points(pw.m[1:126, 1:992, "grey", 3], pw.sd[1:126, 1:992, "grey", 3]^2, pch = 20,
       col = adjustcolor("green3", alpha = 0.4))
points(pw.m[1:126, 1:992, "grey", 4], pw.sd[1:126, 1:992, "grey", 4]^2, pch = 20,
       col = adjustcolor("cyan3", alpha = 0.4))
points(pw.m[1:126, 1:992, "grey", 5], pw.sd[1:126, 1:992, "grey", 5]^2, pch = 20,
       col = adjustcolor("purple", alpha = 0.4))
points(pw.m[1:126, 1:992, "grey", 6], pw.sd[1:126, 1:992, "grey", 6]^2, pch = 20,
       col = adjustcolor("magenta3", alpha = 0.4))

points(pw.m[1:126, 1:992, "grey", 12], pw.sd[1:126, 1:992, "grey", 12]^2, pch = 20,
       col = "gold")
points(pw.m[1:126, 1:992, "grey", 11], pw.sd[1:126, 1:992, "grey", 11]^2, pch = 20,
       col = "orange")

abline(line(pw.m[1:126, 1:992, "grey", 1], pw.sd[1:126, 1:992, "grey", 1]^2), col = "black")
abline(line(pw.m[1:126, 1:992, "grey", 2], pw.sd[1:126, 1:992, "grey", 2]^2), col = adjustcolor("blue", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 3], pw.sd[1:126, 1:992, "grey", 3]^2), col = adjustcolor("green3", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 4], pw.sd[1:126, 1:992, "grey", 4]^2), col = adjustcolor("cyan3", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 5], pw.sd[1:126, 1:992, "grey", 5]^2), col = adjustcolor("purple", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 6], pw.sd[1:126, 1:992, "grey", 6]^2), col = adjustcolor("magenta3", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 7], pw.sd[1:126, 1:992, "grey", 7]^2), col = adjustcolor("orange", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 8], pw.sd[1:126, 1:992, "grey", 8]^2), col = adjustcolor("gold", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 9], pw.sd[1:126, 1:992, "grey", 9]^2), col = adjustcolor("yellow", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 10], pw.sd[1:126, 1:992, "grey", 10]^2), col = adjustcolor("magenta3", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 11], pw.sd[1:126, 1:992, "grey", 11]^2), col = adjustcolor("magenta3", alpha = 0.4))
abline(line(pw.m[1:126, 1:992, "grey", 12], pw.sd[1:126, 1:992, "grey", 12]^2), col = adjustcolor("magenta3", alpha = 0.4))

zz <- lapply(dimnames(pw.m)[[4]], 
             function (dt) coef(line(pw.m[511:638, 993:1996, "grey", dt], pw.sd[511:638, 993:1996, "grey", dt]^2))[2])
