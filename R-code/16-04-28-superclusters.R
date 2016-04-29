
library("IO.Pixels")

fpath <- "./Notes/Superclusters/fig/"
bp <- readRDS("./Models/Simple-parametric/Bad-px-new-thresholds.rds")

load.pixel.means()
res <- readRDS("./Models/Simple-parametric/Residuals-simple-parametric.rds")
ff <- readRDS("./Other-data/Flat-field-corrected.rds")

####################################################################################################

# extend identification of bright/dim pixels to black image. (May have to ignore s.bright pixels)
# probably label these as different type, not bright/dim, to distinguish from behaviour under spot

# revise thresholding classification in line with sd.levels()
# (or revise plot colours in line with thresholding classification)

# TRY USING MEAN +-2SD AS CUTOFF FOR BRIGHT/DIM PIXELS (ONLY INCLUDE WHEN PART OF SUPERCLUSTER?)

mean(res[,,"grey", "160314"]) + (c(-2,2) * sd(res[,,"grey", "160314"]))

####################################################################################################
# GET BAD PIXELS ACROSS ALL IMAGES                                                              ####
# approx. 8 minutes to run across 11 image sets
    {
        # bp <- lapply(dimnames(pw.m)[[4]], bad.px)
        # saveRDS(bp, "./Models/Simple-parametric/Bad-px-new-thresholds.rds")
    }

####################################################################################################

# CLUSTERS & SUPERCLUSTERS                                                                      ####
plot.bad.px(bp$"160314")

bpm <- superclusters(cluster.by.type(enhance.spots(bp$"160314")))
plot.bad.px(bpm)

{
    # no enlargements to screen spots in 160314 acquisitions, so type should be identical before & after:
    all(bpm[order(bpm$x, bpm$y), c("x", "y", "type")] == bp$"160314"[order(bp$"160314"$row, bp$"160314"$col),])
} # checked & confirmed all pixels & pixel types retained correctly

length(unique(bpm$sc.id))   # 123 superclusters identified in this acquisition, including slightly bright px
head(table(bpm$sc.id, bpm$type))

table(qq$type, useNA = "ifany")

#----------------------------------------------------------------------------------------
# plot all superclusters identified: bright & dim px from grey/white only
{
    pdf(paste0(fpath, "Superclusters-excl-black-thresholds.pdf"))
    par(mfrow = c(4,4), mar = c(2,2,1,1))
    for (i in unique(bpm$sc.id[!is.na(bpm$sc.id)])) {
        # get range of cells to plot
        focus <- get.focus(bpm[which(bpm$sc.id == i),], 7)
        
        for (col in dimnames(pw.m)[[3]]) {
            pixel.image(pw.m[,, col, "160314"], xlim = range(focus[,1]), ylim = range(focus[,2]), cex.axis = 0.6)
            text(focus, cex = 0.5, labels = round(pw.m[,, col, "160314"][focus]/1000,0))
            points(bpm[which(bpm$sc.id == i),1:2], pch = 0, cex = 2)
        }
        pixel.image(ff[,,"160314"], xlim = range(focus[,1]), ylim = range(focus[,2]), cex.axis = 0.6)
        text(focus, cex = 0.5, labels = round(ff[,, "160314"][focus]/1000,0))
        points(bpm[which(bpm$sc.id == i),1:2], pch = 0, cex = 2)
    }
    dev.off()
}

{
    pdf(paste0(fpath, "Supercluster-residuals-excl-black-thresholds.pdf"))
    par(mfrow = c(4,4), mar = c(2,2,1,1))
    for (i in unique(bpm$sc.id[!is.na(bpm$sc.id)])) {
        # get range of cells to plot
        focus <- get.focus(bpm[which(bpm$sc.id == i),], 7)
        
        for (col in dimnames(res)[[3]]) {
            pixel.image(res[,, col, "160314"], xlim = range(focus[,1]), ylim = range(focus[,2]), cex.axis = 0.6)
            text(focus, cex = 0.5, labels = round(res[,, col, "160314"][focus]/1000,0))
            points(bpm[which(bpm$sc.id == i),1:2], pch = 0, cex = 2)
        }
        pixel.image(ff[,,"160314"], xlim = range(focus[,1]), ylim = range(focus[,2]), cex.axis = 0.6)
        text(focus, cex = 0.5, labels = round(ff[,, "160314"][focus]/1000,0))
        points(bpm[which(bpm$sc.id == i),1:2], pch = 0, cex = 2)
    }
    dev.off()
}
#----------------------------------------------------------------------------------------
# supercluster composition: what does each supercluster contain? Patterns?
bpm <- merge(bpm, 
             data.frame(sc.id = unique(bpm$sc.id[!is.na(bpm$sc.id)]),
                        sc.comp = sapply(split(levels(bpm$type)[(which(apply(as.matrix(table(bpm$sc.id, bpm$type)), 1, ">", 0))-1) %% 11 + 1],
                                               which(apply(as.matrix(table(bpm$sc.id, bpm$type)), 1, ">", 0)) %/% 11),
                                         paste, collapse = ", ")),
             all = T, sort = F)[, c(colnames(bpm), "sc.comp")]

# plot superclusters by composition (values shown: flat-fielded values)
{
    pdf(paste0(fpath, "Superclusters-by-px-type.pdf"))
    par(mfrow = c(6,6), mar = c(2,2,1,1))
    for (i in unique(bpm$sc.id[!is.na(bpm$sc.id)][order(bpm$sc.comp)[!is.na(bpm$sc.id)]])) {
        # get range of cells to plot
        focus <- get.focus(bpm[which(bpm$sc.id == i),], 3)
        
        plot.bad.px(bpm, block = c("edge"), cex = 2, xlim = range(focus[,1]), ylim = range(focus[,2]), cex.axis = 0.7)
        text(focus, cex = 0.5, labels = round(ff[,, "160314"][focus]/1000,0))
        points(bpm[which(bpm$sc.id == i),1:2], pch = 0, cex = 2)
    }
    dev.off()
}

bpm <- classify.states(bpm)

table(bpm$cat, bpm$type, useNA = "ifany")

bpm[bpm$cat == "sc-no response",]

incl <- c("sc-dead", "sc-no-response", "sc-hot", "sc-bright", "c-dead", "c-no response", "c-hot", "c-v.bright", 
          "c-bright", "c-v.dim", "c-dim", "s-dead", "s-no response", "s-hot", "s-v.bright", "s-bright", 
          "s-v.dim", "s-dim")

cat.pch <- c(16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 1,1,1,1,1,1)
plot.bad.px(bpm[bpm$cat %in% incl,], pch = cat.pch[bpm$cat[bpm$cat %in% incl]])

####################################################################################################

# STATE TRANSITIONS                                                                             ####

# extend screen spots to include adjacent slightly dim pixels
bp2 <- lapply(bp, enhance.spots)

# cluster all images
cl <- lapply(bp2, cluster.by.type)

# identify superclusters
sc <- lapply(cl, superclusters)

# finally, classify states according to supercluster composition
bpm <- lapply(sc, classify.states)

bpm <- list()
bpm[["141009"]] <- classify.states(sc[["141009"]])
bpm[["141118"]] <- classify.states(sc[["141118"]])
bpm[["141217"]] <- classify.states(sc[["141217"]])
bpm[["150108"]] <- classify.states(sc[["150108"]])
bpm[["150113"]] <- classify.states(sc[["150113"]])
bpm[["150126"]] <- classify.states(sc[["150126"]])
bpm[["150529"]] <- classify.states(sc[["150529"]])  # failed - not sure why! problem with sc composition
bpm[["150730"]] <- classify.states(sc[["150730"]])
bpm[["150828"]] <- classify.states(sc[["150828"]])
bpm[["151015"]] <- classify.states(sc[["151015"]])
bpm[["160314"]] <- classify.states(sc[["160314"]])

# compare to get state transitions
incl <- c("sc-dead", "sc-no-response", "sc-hot", "sc-bright", "c-dead", "c-no response", "c-hot", "c-v.bright", 
          "c-bright", "c-v.dim", "c-dim", "s-dead", "s-no response", "s-hot", "s-v.bright", "s-bright", 
          "s-v.dim", "s-dim")

bp.match <- list()
bp.match[[1]] <- droplevels(merge(bpm[[1]][bpm[[1]]$cat %in% incl,], bpm[[2]][bpm[[2]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".141009", ".141118")))
bp.match[[2]] <- droplevels(merge(bpm[[2]][bpm[[2]]$cat %in% incl,], bpm[[3]][bpm[[3]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".141118", ".141217")))
bp.match[[3]] <- droplevels(merge(bpm[[3]][bpm[[3]]$cat %in% incl,], bpm[[4]][bpm[[4]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".141217", ".150108")))
bp.match[[4]] <- droplevels(merge(bpm[[4]][bpm[[4]]$cat %in% incl,], bpm[[5]][bpm[[5]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".150108", ".150113")))
bp.match[[5]] <- droplevels(merge(bpm[[5]][bpm[[5]]$cat %in% incl,], bpm[[6]][bpm[[6]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".150113", ".150126")))
bp.match[[6]] <- droplevels(merge(bpm[[6]][bpm[[6]]$cat %in% incl,], bpm[[7]][bpm[[7]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".150126", ".150730")))
bp.match[[7]] <- droplevels(merge(bpm[[7]][bpm[[7]]$cat %in% incl,], bpm[[8]][bpm[[8]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".150730", ".150828")))
bp.match[[8]] <- droplevels(merge(bpm[[8]][bpm[[8]]$cat %in% incl,], bpm[[9]][bpm[[9]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".150828", ".151015")))
bp.match[[9]] <- droplevels(merge(bpm[[9]][bpm[[9]]$cat %in% incl,], bpm[[10]][bpm[[10]]$cat %in% incl,], by = c("x", "y"), all = T, suffix = c(".151015", ".160314")))


table("141009" = bp.match[[1]]$cat.141009, "141118" = bp.match[[1]]$cat.141118, useNA = "ifany")
table("141118" = bp.match[[2]]$cat.141118, "141217" = bp.match[[2]]$cat.141217, useNA = "ifany")
table("141217" = bp.match[[3]]$cat.141217, "150108" = bp.match[[3]]$cat.150108, useNA = "ifany")
table("150108" = bp.match[[4]]$cat.150108, "150113" = bp.match[[4]]$cat.150113, useNA = "ifany")
table("150113" = bp.match[[5]]$cat.150113, "150126" = bp.match[[5]]$cat.150126, useNA = "ifany")
table("150126" = bp.match[[6]]$cat.150126, "150730" = bp.match[[6]]$cat.150730, useNA = "ifany")
table("150730" = bp.match[[7]]$cat.150730, "150828" = bp.match[[7]]$cat.150828, useNA = "ifany")
table("150828" = bp.match[[8]]$cat.150828, "151015" = bp.match[[8]]$cat.151015, useNA = "ifany")
table("151015" = bp.match[[9]]$cat.151015, "160314" = bp.match[[9]]$cat.160314, useNA = "ifany")


####################################################################################################

# FURTHER THRESHOLD REVISION?                                                                   ####

# revise inner threshold to bring into line with colours used in plotting?
# use mean + 1 or 2 SD to get cutoff? (or use MAD instead?)
{
    hist(res[,,"grey", "160314"], breaks = "fd", ylim = c(0,50), xlim = c(-10000,10000))
    
    abline(v = mean(res[,,"grey", "160314"]), col = "red")
    abline(v = mean(res[,,"grey", "160314"]) + (c(-1,1) * sd(res[,,"grey", "160314"])), col = "red", lty = 2)
    abline(v = mean(res[,,"grey", "160314"]) + (c(-2,2) * sd(res[,,"grey", "160314"])), col = "red", lty = 3)
    
    abline(v = median(res[,,"grey", "160314"]), col = "blue")
    abline(v = median(res[,,"grey", "160314"]) + (c(-1,1) * sd(res[,,"grey", "160314"])), col = "blue", lty = 2)
    abline(v = median(res[,,"grey", "160314"]) + (c(-1,1) * mad(res[,,"grey", "160314"])), col = "blue", lty = 3)
    
}

####################################################################################################

# COMPARE THIS BAD PIXEL MAP TO THE 'OFFICIAL' VERSION                                          ####

# REVISE THRESHOLDS: ADD ANOTHER STRIP AT INTERIOR?                                             ####

# NEXT STEP: LOOK FOR COLUMNS OF BRIGHT/DIM PIXELS                                              ####

# LINK COLUMNS WITH SUPERCLUSTERS?                                                              ####

# APPLY NEW METHOD TO OLD DATA SET                                                              ####