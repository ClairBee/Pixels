
# DETECTION OF SLIGHTLY BRIGHT (OR DIM) LINES

library("IO.Pixels"); library("CB.Misc")
if (getwd() != "/home/clair/Documents/Pixels") {setwd("/home/clair/Documents/Pixels")}
fpath <- "./Notes/Line-detection/Lines-detected/fig/"

load.pixel.means()

####################################################################################################

# THURSDAY MAY 12TH: DONE!

{
    # run function over all black images: bright & dim
    # run function over all grey images: bright & dim
    # run function over all white images: bright & dim
}

# for all lines found, plot
{
    #   - development over time at each power level
    #   - overplot convolutions over time: is difference between line & neighbours constant?
    #   - try to associate lines with superclusters
}

{
    # run over Jay's graduated images: are lines identified same at all settings?
    # is behaviour relative to neighbours same at all settings?
}

# apply function to as many old images as possible
# repeat the above for any & all lines found

# write all this up & send it to Julia.

#-------------------------------------------------------------------------------------------------

# possible further cleaning step... convolve with panel-edge kernel
# if that score is lower than bright kernel, the panel edge is bright.

#-------------------------------------------------------------------------------------------------

# also have a quick look at matching that diagonal pattern across images
#   - perhaps use differences across rows/columns? Same pattern of increases/decreases?
#   - may need to normalise (but prob. not in black images?)
#   - try to find a particularly egregious area & focus on details
#   - maybe pattern of gradients is more useful than actual change in value?

####################################################################################################

# RUN CONVOLUTIONS                                                                              ####

bright.b <- alply(pw.m[,,"black", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
                  sm.size = 11, min.length = 6, .dims = T)

bright.g <- alply(pw.m[,,"grey", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
                  sm.size = 11, min.length = 6, .dims = T) 

bright.w <- alply(pw.m[,,"white", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
                  sm.size = 11, min.length = 6, .dims = T) 

dim.b <- alply(pw.m[,,"black", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
               sm.size = 11, min.length = 6, dim.lines = T, .dims = T)

dim.g <- alply(pw.m[,,"grey", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
               sm.size = 11, min.length = 6, dim.lines = T, .dims = T)

dim.w <- alply(pw.m[,,"white", ], 3, find.lines, k.size = 5, threshold.at = 5500, 
               sm.size = 11, min.length = 6, dim.lines = T, .dims = T)

# ran same over rows: no bright rows

# plot all bright lines identified
{
    pdf(paste0(fpath, "bright-lines-identified.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, bright.w[[1]], add = F, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        for (i in 2:length(bright.w)) {
            image(1:1996, 1:1996, bright.w[[i]], add = T, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        }
        for (i in 1:length(bright.g)) {
            image(1:1996, 1:1996, bright.g[[i]], add = T, col = c(NA, "cyan3"), breaks = c(-0.5, 0.5, 99))
        }
        for (i in 1:length(bright.b)) {
            image(1:1996, 1:1996, bright.b[[i]], add = T, col = c(NA, "red"), breaks = c(-0.5, 0.5, 99))
        }
        dev.off()
    }
    pdf(paste0(fpath, "bright-lines-identified-with-panels.pdf"))
    draw.panels()
}

# plot all dim lines identified
{
    pdf(paste0(fpath, "dim-lines-identified.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, dim.w[[1]], add = F, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        for (i in 2:length(dim.w)) {
            image(1:1996, 1:1996, dim.w[[i]], add = T, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        }
        for (i in 1:length(dim.g)) {
            image(1:1996, 1:1996, dim.g[[i]], add = T, col = c(NA, "cyan3"), breaks = c(-0.5, 0.5, 99))
        }
        for (i in 1:length(dim.b)) {
            image(1:1996, 1:1996, dim.b[[i]], add = T, col = c(NA, "red"), breaks = c(-0.5, 0.5, 99))
        }
        dev.off()
    }
    pdf(paste0(fpath, "dim-lines-identified-with-panels.pdf"))
    draw.panels()
}

####################################################################################################

# SUMMARISE LINES FOUND IN ALL IMAGES                                                           ####

df.bright.b <- lapply(bright.b, summarise.lines)
df.bright.g <- lapply(bright.g, summarise.lines)
df.bright.w <- lapply(bright.w, summarise.lines)

# get master list of bright lines as CSV
{
    bright.cols <- ddply(rbind(rbind.fill(df.bright.b),
                               rbind.fill(df.bright.g),
                               rbind.fill(df.bright.w)),
                         .(col, panel), summarise, min.y = min(ymin), max.y = max(ymax))
    
    # tabulate appearances
    for (i in 1:length(df.bright.b)) {
        bright.cols[, paste0("imb", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(bright.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.bright.b[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
        bright.cols[, paste0("img", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(bright.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.bright.g[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
        bright.cols[, paste0("imw", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(bright.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.bright.w[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
    }
    
    # set output to create easy table
    bright.cols[, 5 + (0:11) * 3][!is.na(bright.cols[, 5 + (0:11) * 3])] <- "$\\bullet$"
    bright.cols[, 6 + (0:11) * 3][!is.na(bright.cols[, 6 + (0:11) * 3])] <- "\\textcolor{BlueGreen}{$\\bullet$}"
    bright.cols[, 7 + (0:11) * 3][!is.na(bright.cols[, 7 + (0:11) * 3])] <- "\\textcolor{gold}{$\\bullet$}"
    
    bright.cols[is.na(bright.cols)] <- ""
    
    # add flag to indicate whether at panel edge
    bright.cols$spedge <- ""
    bright.cols$spedge[(bright.cols$col+2) %% 128 < 3] <- "$\\mathbf{\\times}$" 
    bright.cols$spedge[(bright.cols$col+2) %% 128 > 126] <- "$\\mathbf{\\times}$" 
    bright.cols$spedge[bright.cols$col <= 5 ] <- "$\\mathbf{\\times}$" 
    bright.cols$spedge[bright.cols$col >= 1991 ] <- "$\\mathbf{\\times}$" 
    
    write.csv(bright.cols, paste0(fpath, "bright-columns.csv"), quote = F, row.names = F)
}

#-------------------------------------------------------------------------------------------

df.dim.b <- lapply(dim.b, summarise.lines)
df.dim.g <- lapply(dim.g, summarise.lines)
df.dim.w <- lapply(dim.w, summarise.lines)

# get master list of dim lines as CSV
{
    dim.cols <- ddply(rbind.fill(rbind.fill(df.dim.b),
                                 rbind.fill(df.dim.g),
                                 rbind.fill(df.dim.w)),
                      .(col, panel), summarise, min.y = min(ymin), max.y = max(ymax))
    
    # tabulate appearances
    for (i in 1:length(df.dim.b)) {
        dim.cols[, paste0("imb", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(dim.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.dim.b[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
        dim.cols[, paste0("img", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(dim.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.dim.g[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
        dim.cols[, paste0("imw", letters[i])] <- sapply(gsub("^\\s+|\\s+$", "", apply(dim.cols[,1:2], 1, paste, collapse = "")),
                                                           match,
                                                           gsub("^\\s+|\\s+$", "", apply(df.dim.w[[i]][,1:2], 1, paste, collapse = "")),
                                                           nomatch = "")
    }
    
    # set output to create easy table
    dim.cols[, 5 + (0:11) * 3][!is.na(dim.cols[, 5 + (0:11) * 3])] <- "$\\bullet$"
    dim.cols[, 6 + (0:11) * 3][!is.na(dim.cols[, 6 + (0:11) * 3])] <- "\\textcolor{BlueGreen}{$\\bullet$}"
    dim.cols[, 7 + (0:11) * 3][!is.na(dim.cols[, 7 + (0:11) * 3])] <- "\\textcolor{gold}{$\\bullet$}"
    
    dim.cols[is.na(dim.cols)] <- ""
    
    # add flag to indicate whether at panel edge
    dim.cols$spedge <- ""
    dim.cols$spedge[(dim.cols$col+2) %% 128 < 3] <- "$\\mathbf{\\times}$" 
    dim.cols$spedge[(dim.cols$col+2) %% 128 > 126] <- "$\\mathbf{\\times}$" 
    dim.cols$spedge[dim.cols$col <= 5 ] <- "$\\mathbf{\\times}$" 
    dim.cols$spedge[dim.cols$col >= 1991 ] <- "$\\mathbf{\\times}$" 
    
    write.csv(dim.cols, paste0(fpath, "dim-columns.csv"), quote = F, row.names = F)
}

####################################################################################################

# LINE BEHAVIOUR AT GRADUATED POWER SETTINGS                                                    ####

# import data
jay.load <- function(lvl) {
    
    lvl <- toString(lvl)
    jpath = paste0("/home/clair/Documents/Pixels/Other-data/Other-images/line_investigation/", lvl, "ua/")

    jay.tiffs <- list.files(jpath, pattern = "\\.tif$")
    
    # create array to hold loaded data
    ims <- array(dim = c(1996, 1996, length(jay.tiffs)))

    for (i in 1:length(jay.tiffs)) {
        tmp <- readTIFF(paste0(jpath, jay.tiffs[i]), as.is = T)
        
        # transpose & rotate to get image right way up
        ims[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
    }
    
    # get pixelwise mean
    m <- apply(ims, c(1, 2), mean)
    
    # assign to new object
    assign(paste0("ua", lvl), m, envir = .GlobalEnv)
}
lines <- lapply(c(20, 40, 60, 80, 100), jay.load)

lines <- lapply(list(ua20, ua40, ua60, ua80, ua100), find.lines, k.size = 5, threshold.at = 5500, sm.size = 11, min.length = 6)

lapply(lines, summarise.lines)

dim.lines <- lapply(list(ua20 = ua20, ua40 = ua40, ua60 = ua60, ua80 = ua80, ua100 = ua100), 
                    find.lines, dim.lines = T, k.size = 5, threshold.at = 5500, sm.size = 11, min.length = 6)
lapply(dim.lines, summarise.lines)

# plot all bright lines identified - basically, all panel edges
{
    pdf(paste0(fpath, "bright-lines-by-power-setting.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, lines[[5]], add = F, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        for (i in 4:1) {
            image(1:1996, 1:1996, lines[[i]], add = T, col = c(NA, c("orange", "red", "purple", "dodgerblue")[i]), breaks = c(-0.5, 0.5, 99))
        }
        dev.off()
    }
    pdf(paste0(fpath, "bright-lines-by-power-setting-with-panels.pdf"))
    draw.panels()
}

# plot all dim lines identified - basically, all panel edges
{
    pdf(paste0(fpath, "dim-lines-by-power-setting.pdf")); {
        par(mar = c(2,2,1,1))
        image(1:1996, 1:1996, dim.lines[[5]], add = F, col = c(NA, "gold"), breaks = c(-0.5, 0.5, 99))
        for (i in 4:1) {
            image(1:1996, 1:1996, dim.lines[[i]], add = T, col = c(NA, c("orange", "red", "purple", "dodgerblue")[i]), breaks = c(-0.5, 0.5, 99))
        }
        dev.off()
    }
    pdf(paste0(fpath, "dim-lines-by-power-setting-with-panels.pdf"))
    draw.panels()
}

# still same two columns identified.

#-----------------------------------------------------------------------------------------

plot.cols <- c("black", "deepskyblue", "blue", "purple", "magenta3", "red", "orange", "gold")

# transect plots
pdf(paste0(fpath, "bright-line-upper.pdf")); {
    par(mar = c(2,2,1,1))
    
    plot(ua100[427,], xlim = c(993, 1996), ylim = c(0,65535), type = "l", col = plot.cols[8], xlab = "", ylab = "")
    lines(pw.m[427,,"white", "160430"], col = plot.cols[7])
    lines(ua80[427,], col = plot.cols[6])
    lines(ua60[427,], col = plot.cols[5])
    lines(ua40[427,], col = plot.cols[4])
    lines(ua20[427,], col = plot.cols[3])
    lines(pw.m[427,,"grey", "160430"], col = plot.cols[2])
    lines(pw.m[427,,"black", "160430"], col = plot.cols[1])
    
    text(rep(1020, 8), 1000 * c(7, 15, 20, 30, 40, 45, 50.5, 58), 
         label = c("black", "grey", "ua20", "ua40", "ua60", "ua80", "white", "ua100"),
         col = plot.cols)
    
    dev.off()
} # upper

pdf(paste0(fpath, "bright-line-lower.pdf")); {
    par(mar = c(2,2,1,1))
    
    plot(ua100[809,], xlim = c(1, 992), ylim = c(0,65535), type = "l", col = plot.cols[8], xlab = "", ylab = "")
    lines(pw.m[809,,"white", "160430"], col = plot.cols[7])
    lines(ua80[809,], col = plot.cols[6])
    lines(ua60[809,], col = plot.cols[5])
    lines(ua40[809,],col = plot.cols[4])
    lines(ua20[809,], col = plot.cols[3])
    lines(pw.m[809,,"grey", "160430"], col = plot.cols[2])
    lines(pw.m[809,,"black", "160430"], col = plot.cols[1])
    
    text(rep(970, 8), 1000 * c(7, 15.5, 20, 30, 40, 46, 51, 58.5), 
         label = c("black", "grey", "ua20", "ua40", "ua60", "ua80", "white", "ua100"),
         col = plot.cols)
    
    dev.off()
} # lower


# transects with neighbours
{
    pdf(paste0(fpath, "upper-by-power-detail.pdf")); {
        par(mar = c(2,2,1,1), mfrow = c(4, 2))
        
        plot(ua100[427,], xlim = c(993, 1996), ylim = c(54000, 59000), type = "l", col = plot.cols[7], xlab = "", ylab = "")
        lines(ua100[426,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua100[428,], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 59000, "ua100")
        
        plot(pw.m[427,, "white", "160430"], xlim = c(993, 1996), ylim = c(45000, 50000), type = "l", col =  plot.cols[7], xlab = "", ylab = "")
        lines(pw.m[426,, "white", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[428,, "white", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 50000, "white")
        
        plot(ua80[427,], xlim = c(993, 1996), ylim = c(45000, 50000), type = "l", col = plot.cols[6], xlab = "", ylab = "")
        lines(ua80[426,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua80[428,], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 50000, "ua80")
        
        plot(ua60[427,], xlim = c(993, 1996), ylim = c(35000, 40000), type = "l", col = plot.cols[5], xlab = "", ylab = "")
        lines(ua60[426,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua60[428,], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 40000, "ua60")
        
        plot(ua40[427,], xlim = c(993, 1996), ylim = c(25000, 30000), type = "l", col = plot.cols[4], xlab = "", ylab = "")
        lines(ua40[426,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua40[428,], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 30000, "ua40")
        
        plot(ua20[427,], xlim = c(993, 1996), ylim = c(15000, 20000), type = "l", col = plot.cols[3], xlab = "", ylab = "")
        lines(ua20[426,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua20[428,], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 20000, "ua20")
        
        plot(pw.m[427,, "grey", "160430"], xlim = c(993, 1996), ylim = c(15000, 20000), type = "l", col = plot.cols[2], xlab = "", ylab = "")
        lines(pw.m[426,, "grey", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[428,, "grey", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 20000, "grey")
        
        plot(pw.m[427,, "black", "160430"], xlim = c(993, 1996), ylim = c(3000, 8000), type = "l", col = plot.cols[1], xlab = "", ylab = "")
        lines(pw.m[426,, "black", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[428,, "black", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(1010, 8000, "black")
        
        dev.off()
    }
    
    pdf(paste0(fpath, "lower-by-power-detail.pdf")); {
        par(mar = c(2,2,1,1), mfrow = c(4, 2))
        
        plot(ua100[809,], xlim = c(1, 992), ylim = c(54000, 59000), type = "l", col = plot.cols[8], xlab = "", ylab = "")
        lines(ua100[808,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua100[810,], col = adjustcolor("green3", alpha = 0.4))
        text(15, 59000, "ua100")
        
        plot(pw.m[809,, "white", "160430"], xlim = c(1, 992), ylim = c(45000, 50000), type = "l", col = plot.cols[7], xlab = "", ylab = "")
        lines(pw.m[808,, "white", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[810,, "white", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(15, 50000, "white")
        
        plot(ua80[809,], xlim = c(1, 992), ylim = c(45000, 50000), type = "l", col = plot.cols[6], xlab = "", ylab = "")
        lines(ua80[808,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua80[810,], col = adjustcolor("green3", alpha = 0.4))
        text(15, 50000, "ua80")
        
        plot(ua60[809,], xlim = c(1, 992), ylim = c(35000, 40000), type = "l", col = plot.cols[5], xlab = "", ylab = "")
        lines(ua60[808,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua60[810,], col = adjustcolor("green3", alpha = 0.4))
        text(15, 40000, "ua60")
        
        plot(ua40[809,], xlim = c(1, 992), ylim = c(25000, 30000), type = "l", col = plot.cols[4], xlab = "", ylab = "")
        lines(ua40[808,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua40[810,], col = adjustcolor("green3", alpha = 0.4))
        text(15, 30000, "ua40")
        
        plot(ua20[809,], xlim = c(1, 992), ylim = c(15000, 20000), type = "l", col = plot.cols[3], xlab = "", ylab = "")
        lines(ua20[808,], col = adjustcolor("cyan3", alpha = 0.4))
        lines(ua20[810,], col = adjustcolor("green3", alpha = 0.4))
        text(15, 20000, "ua20")
        
        plot(pw.m[809,, "grey", "160430"], xlim = c(1, 992), ylim = c(15000, 20000), type = "l", col = plot.cols[2], xlab = "", ylab = "")
        lines(pw.m[808,, "grey", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[810,, "grey", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(15, 20000, "grey")
        
        plot(pw.m[809,, "black", "160430"], xlim = c(1, 992), ylim = c(3000, 8000), type = "l", col = plot.cols[1], xlab = "", ylab = "")
        lines(pw.m[808,, "black", "160430"], col = adjustcolor("cyan3", alpha = 0.4))
        lines(pw.m[810,, "black", "160430"], col = adjustcolor("green3", alpha = 0.4))
        text(15, 8000, "black")
        
        dev.off()
    }
}

# convolution transects
conv <- lapply(list(ua100 = ua100, w = pw.m[,,"white", "160430"], ua80 = ua80, ua60 = ua60, ua40 = ua40,
                    ua20 = ua20, g = pw.m[,,"grey", "160430"]),  b = pw.m[,,"black", "160430"], convolve.lines, k.size = 5)

{
    pdf(paste0(fpath, "convolutions-by-power-upper.pdf"), height = 3, width = 7); {
        par(mar = c(2,2,1,1))
        plot(conv$ua100[427,], type = "l", col = plot.cols[8], xlim = c(993, 1996), ylim = c(-1000, 20000), xlab = "", ylab = "")
        for (i in 2:8) {
            lines(conv[[i]][427,], col = rev(plot.cols)[i])
        }
        dev.off()
    }

    pdf(paste0(fpath, "convolutions-by-power-lower.pdf"), height = 3, width = 7); {
        par(mar = c(2,2,1,1))
        plot(conv$ua100[809,], type = "l", col = plot.cols[8], xlim = c(1, 992), ylim = c(-1000, 20000), xlab = "", ylab = "")
        for (i in 2:8) {
            lines(conv[[i]][809,], col = rev(plot.cols)[i])
        }
        dev.off()
    }
    
    pdf(paste0(fpath, "convolutions-by-power-normal.pdf"), height = 3, width = 7); {
        par(mar = c(2,2,1,1))
        plot(conv$ua100[809,], type = "l", col = plot.cols[8], xlim = c(993,1996), ylim = c(-1000, 20000), xlab = "", ylab = "")
        for (i in 2:8) {
            lines(conv[[i]][809,], col = rev(plot.cols)[i])
        }
        dev.off()
    }   
    
}

# plot differences?
{
    plot(ua100[427,] - ua100[426,], xlim = c(993, 1996), ylim = c(-1000, 3000), type = "l", col = plot.cols[8], xlab = "", ylab = "")
    lines(pw.m[427,,"white", "160430"] - pw.m[426,,"white", "160430"], col = plot.cols[7])
    lines(ua80[427,] - ua80[426,], col = plot.cols[6])
    lines(ua60[427,] - ua60[426,], col = plot.cols[5])
    lines(ua40[427,] - ua40[426,], col = plot.cols[4])
    lines(ua20[427,] - ua20[426,], col = plot.cols[3])
    lines(pw.m[427,,"grey", "160430"] - pw.m[426,,"grey", "160430"], col = plot.cols[2])
    lines(pw.m[427,,"black", "160430"] - pw.m[426,,"black", "160430"], col = plot.cols[1])
    
    
    sd(unlist(lapply(list(ua100 = ua100, w = pw.m[,,"white", "160430"], ua80 = ua80, ua60 = ua60, ua40 = ua40,
                ua20 = ua20, g = pw.m[,,"grey", "160430"], b = pw.m[,,"black", "160430"]),
           function (x) median(c(x[427, 1196:1996] - x[426, 1196:1996], 
                                 x[427, 1196:1996] - x[428, 1196:1996])))))
    
    sd(unlist(lapply(list(ua100 = ua100, w = pw.m[,,"white", "160430"], ua80 = ua80, ua60 = ua60, ua40 = ua40,
                          ua20 = ua20, g = pw.m[,,"grey", "160430"], b = pw.m[,,"black", "160430"]),
                     function (x) median(c(x[427, 1:992] - x[426, 1:992], 
                                           x[427, 1:992] - x[428, 1:992])))))
    
    sd(unlist(lapply(list(ua100 = ua100, w = pw.m[,,"white", "160430"], ua80 = ua80, ua60 = ua60, ua40 = ua40,
                          ua20 = ua20, g = pw.m[,,"grey", "160430"], b = pw.m[,,"black", "160430"]),
                     function (x) median(c(x[809, 1:992] - x[808, 1:992], 
                                           x[809, 1:992] - x[810, 1:992])))))
    
    sd(unlist(lapply(list(ua100 = ua100, w = pw.m[,,"white", "160430"], ua80 = ua80, ua60 = ua60, ua40 = ua40,
                          ua20 = ua20, g = pw.m[,,"grey", "160430"], b = pw.m[,,"black", "160430"]),
                     function (x) median(c(x[809, 1196:1996] - x[808, 1196:1996], 
                                           x[809, 1196:1996] - x[810, 1196:1996])))))
}

####################################################################################################

# PLOTS OF BEHAVIOUR AROUND/AFTER SPIKE                                                         ####

# development of 'lift' around midline
{
    {
        plot(pw.m[427, , "black", "151015"], type = "l", xlim = c(992, 1996), ylim = c(4000,7000), col = adjustcolor("cyan3", alpha = 0.5), xlab = "", ylab = "")
        lines(pw.m[427, , "black", "160314"], col = adjustcolor("orange", alpha = 0.5))
        lines(pw.m[427, , "black", "160430"])
        
        abline(v = 992.5, lty = 4, col = "red")
    }
    
    pdf(paste0(fpath, "column-diffs-over-time-upper.pdf"), height = 3, width = 7); {
            par(mar = c(2,2,1,1))
            plot(pw.m[427, , "black", "151015"] - pw.m[426, , "black", "151015"], ylim = c(-500,1000), type = "l", xlim = c(993, 1996), col = adjustcolor("cyan3", alpha = 0.5), xlab = "", ylab = "")
            lines(pw.m[426, , "black", "160430"] - pw.m[425, , "black", "160430"], col = adjustcolor("gold", alpha = 0.1))
            lines(pw.m[428, , "black", "160430"] - pw.m[429, , "black", "160430"], col = adjustcolor("gold", alpha = 0.1))
            
            lines(pw.m[427, , "black", "151015"] - pw.m[428, , "black", "151015"], col = adjustcolor("cyan3", alpha = 0.5), xlab = "", ylab = "")
            lines(pw.m[427, , "black", "160314"] - pw.m[426, , "black", "160314"], col = adjustcolor("orange", alpha = 0.5))
            lines(pw.m[427, , "black", "160314"] - pw.m[428, , "black", "160314"], col = adjustcolor("orange", alpha = 0.5))
            lines(pw.m[427, , "black", "160430"] - pw.m[426, , "black", "160430"])
            lines(pw.m[427, , "black", "160430"] - pw.m[428, , "black", "160430"])
            abline(v = 992.5, lty = 4, col = "red")
            dev.off()
        }

    pdf(paste0(fpath, "column-diffs-over-time-lower.pdf"), height = 3, width = 7); {
            par(mar = c(2,2,1,1))
            plot(pw.m[809, , "black", "151015"] - pw.m[808, , "black", "151015"], ylim = c(-500,1000), type = "l", xlim = c(1, 992), col = adjustcolor("cyan3", alpha = 0.5), xlab = "", ylab = "")
            lines(pw.m[810, , "black", "160430"] - pw.m[811, , "black", "160430"], col = adjustcolor("gold", alpha = 0.1))
            lines(pw.m[808, , "black", "160430"] - pw.m[807, , "black", "160430"], col = adjustcolor("gold", alpha = 0.1))
            
            lines(pw.m[809, , "black", "151015"] - pw.m[810, , "black", "151015"],col = adjustcolor("cyan3", alpha = 0.5), xlab = "", ylab = "")
            lines(pw.m[809, , "black", "160314"] - pw.m[808, , "black", "160314"], col = adjustcolor("orange", alpha = 0.5))
            lines(pw.m[809, , "black", "160314"] - pw.m[810, , "black", "160314"], col = adjustcolor("orange", alpha = 0.5))
            lines(pw.m[809, , "black", "160430"] - pw.m[808, , "black", "160430"])
            lines(pw.m[809, , "black", "160430"] - pw.m[810, , "black", "160430"])
            
            abline(v = 992.5, lty = 4, col = "red")
            dev.off()
    }
}

plot(0, type = "n", xlim = c(1, 992), ylim = c(-500, 1000), xlab = "", ylab = "")

matplot(pw.m[809, 1:992, "black",] - pw.m[810, 1:992, "black",], type = "l", ylim = c(-500,1000))
matplot(pw.m[427, 993:1996, "black",] - pw.m[426, 993:1996, "black",], type = "l", ylim = c(-500,1000))

median(c(pw.m[809, 1:177, "black", "160430"] - pw.m[810, 1:177, "black", "160430"],
     pw.m[809, 1:177, "black", "160430"] - pw.m[808, 1:177, "black", "160430"]))
abline(h = 366, col = "red")
median(c(pw.m[809, 500:992, "black", "160430"] - pw.m[810,  500:992, "black", "160430"],
       pw.m[809, 500:992, "black", "160430"] - pw.m[808,  500:992, "black", "160430"]))

plot(pw.m[809, 1:992, "black", "160430"] - pw.m[810,  1:992, "black", "160430"], type = "l", ylim = c(-500, 1000), xlab = "", ylab = "")

plot(pw.m[809, 1:992, "black", "160430"], type = "l", ylim = c(4500, 6000), xlab = "", ylab = "")
lines(pw.m[810, 1:992, "black", "160430"], col = "red")

####################################################################################################

# SUPERCLUSTERS                                                                                 ####

# use bad pixel map from simplest parametric model
bp <- readRDS("./Models/Simple-parametric/Bad-px-new-thresholds.rds")

focus <- get.focus(data.frame(x = c(427, 427), y = c(1195, 1196)), surround = 10)

# upper line
{
    pdf(paste0(fpath, "supercluster-427.pdf")); {
        par(mar = c(2,2,1,1))
        plot(bp$"160314"[,1:2], xlim = c(422, 432), ylim = c(1195, 1205), pch = 0, asp = T, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"black", "160314"][focus]/1000,0))
        rect(426.5, 1195.5, 427.5, 1996.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        plot(bp$"160314"[,1:2], xlim = c(422, 432), ylim = c(1195, 1205), pch = 0, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"grey", "160314"][focus]/1000,0))
        rect(426.5, 1195.5, 427.5, 1996.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        plot(bp$"160314"[,1:2], xlim = c(422, 432), ylim = c(1195, 1205), pch = 0, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"white", "160314"][focus]/1000,0))
        rect(426.5, 1195.5, 427.5, 1996.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        dev.off()
    }
    
    focus <- get.focus(data.frame(x = c(809, 809), y = c(177, 178)), surround = 10)
    
    pdf(paste0(fpath, "supercluster-809.pdf")); {
        par(mar = c(2,2,1,1))
        plot(bp$"160314"[,1:2], xlim = range(focus[,1]), ylim = range(focus[,2]), pch = 0, asp = T, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"black", "160314"][focus]/1000,0))
        rect(808.5, 0.5, 809.5, 178.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        plot(bp$"160314"[,1:2], xlim = range(focus[,1]), ylim = range(focus[,2]), pch = 0, asp = T, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"grey", "160314"][focus]/1000,0))
        rect(808.5, 0.5, 809.5, 178.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        plot(bp$"160314"[,1:2], xlim = range(focus[,1]), ylim = range(focus[,2]), pch = 0, asp = T, cex = 5, lwd = 3, col = bp.colours(block = "edge")[bp$"160314"$type], xlab = "", ylab = "")
        text(focus, cex = 1.2, labels = round(pw.m[,,"white", "160314"][focus]/1000,0))
        rect(808.5, 0.5, 809.5, 178.5, border = "cyan2", lwd = 2, col = adjustcolor("cyan2", alpha = 0.1))
        
        dev.off()
    }
}


####################################################################################################

# RUN OVER OLD DATA                                                                             ####

# try running over old data

# MAY NEED TO REPOSITION CENTRE LINE

tmp.black <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapBlack.tif", as.is = T)
tmp.black <- t(tmp.black[nrow(tmp.black):1, , drop = FALSE])

tmp.grey <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapGrey.tif", as.is = T)
tmp.grey <- t(tmp.grey[nrow(tmp.grey):1, , drop = FALSE])

tmp.white <- readTIFF("~/Documents/Pixels/Other-data/Old-data/130613-bad-pixel-map/BadPixelMapWhite.tif", as.is = T)
tmp.white <- t(tmp.white[nrow(tmp.white):1, , drop = FALSE])

zz <- find.lines(tmp.black)

image(1:2000, 1:2000, zz, col = c(NA, "red"), xlim = c(970,1000), ylim = c(0,300))
# got one!

which(zz > 0, arr.ind = T)

pixel.image(tmp.black, xlim = c(970,1000), ylim = c(0,300))

o.plot(tmp.black[970:1000, 100]); o.plot(tmp.black[970:1000, 103], add = T); o.plot(tmp.black[970:1000, 97], add = T)

#-------------------------------------------------------------------------------
# dim lines?

qq <- find.lines(tmp.black, dim.lines = T)

table(qq) # two lines found
image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999))

summarise.lines(qq)

image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999), xlim = c(190, 220), ylim = c(900, 1300))
pixel.image(tmp.black, xlim = c(190, 220), ylim = c(900, 1300))
o.plot(tmp.black[190:230, 1100]); o.plot(tmp.black[190:230, 1101], add = T); o.plot(tmp.black[190:230, 1099], add = T);
o.plot(tmp.black[208, ], xlim = c(990, 1100))

image(1:2000, 1:2000, qq, col = c(NA, "red"), breaks = c(-0.5, 0.5, 999), xlim = c(1300, 1400), ylim = c(1000, 1100))
pixel.image(tmp.black,  xlim = c(1300, 1400), ylim = c(1000, 1100))
o.plot(tmp.black[1372, ], xlim = c(990, 1100))



# plots of line 1
{
    o.plot(tmp[208, ], xlim = c(990, 1300), ylim = c(0,50000))
    o.plot(tmp.grey[208, ], xlim = c(990, 1300), add = T, col = "grey")
    o.plot(tmp.white[208, ], xlim = c(990, 1300), add = T, col = "gold")
    abline(v = 1000.5, lwd = 2)                 # panel midline
    abline(v = 1247.5, col = "red", lty = 2)    # end of dim line
}

# plots of line 2
{
    o.plot(tmp[1372, ], xlim = c(990, 1100), ylim = c(0,50000))
    o.plot(tmp.grey[1372, ], xlim = c(990, 1100), add = T, col = "grey")
    o.plot(tmp.white[1372, ], xlim = c(990, 1100), add = T, col = "gold")
    abline(v = 1000.5, lwd = 2)                 # panel midline
    abline(v = 1071.5, col = "red", lty = 2)    # end of dim line
}



####################################################################################################

# SCRATCH/MISC                                                                                  ####

# plots of problematic areas
{
    # line segment retained in column 3
    {
        o.plot(pw.m[3, 1940:1970, "black", "160314"], ylim = c(0, 25000))
        o.plot(pw.m[4, 1940:1970, "black", "160314"], add = T, col = adjustcolor("green3", alpha = 0.3))
        o.plot(pw.m[2, 1940:1970, "black", "160314"], add = T, col = adjustcolor("gold", alpha = 0.3))
        
        o.plot(conv[3, 1940:1970], add = T, col = "cyan3")
        abline(h = 5500, col = "red", lty = 2)
        o.plot(sm[3, 1940:1970] * 2500, add = T, col = "blue")
        abline(h = 5.5 * 2500, col = "blue", lty = 2)
        
        pixel.image(pw.m[1:126, , "black", "160314"], xlim = c(0,30), ylim = c(1940,1970))
        rect(2.5, 1950.5, 3.5, 1965.5)
    }
    
    # line segment retained in column 234
    {
        o.plot(pw.m[234, 260:300, "black", "160314"], ylim = c(0, 25000))
        o.plot(pw.m[233, 260:300, "black", "160314"], add = T, col = adjustcolor("green3", alpha = 0.3))
        o.plot(pw.m[232, 260:300, "black", "160314"], add = T, col = adjustcolor("gold", alpha = 0.3))
        
        o.plot(conv[234, 260:300], add = T, col = "cyan3")
        abline(h = 5500, col = "red", lty = 2)
        o.plot(sm[234, 260:300] * 2500, add = T, col = "blue")
        abline(h = 5.5 * 2500, col = "blue", lty = 2)
        
        pixel.image(pw.m[200:300, , "black", "160314"], xlim = c(210,250), ylim = c(260,300))
        rect(233.5, 274.5, 234.5, 293.5)
    }
    
    # line identified in column 427
    {
        o.plot(pw.m[427,, "black", "160314"], ylim = c(-1000, 10000), xlim = c(1100,1996))
        o.plot(conv[427,], add = T, col = "green3")
        abline(h = 5500, col = "red")
        o.plot(sm[427,] * 900, add = T, col = "blue")
        abline(h = 5.5*900, col = "blue", lty = 2)
    }
    
    # line identified in column 809
    {
        o.plot(pw.m[809,, "black", "160314"], ylim = c(-1000, 10000), xlim = c(0,992))
        points(pw.m[810,, "black", "160314"], type = "l", col = "orange")
        points(pw.m[808,, "black", "160314"], type = "l", col = "gold")
        o.plot(conv[809,], add = T, col = "green3")
        abline(h = 5500, col = "red")
        o.plot(sm[809,] * 900, add = T, col = "blue")
        abline(h = 5.5*900, col = "blue", lty = 2)
    }
}
