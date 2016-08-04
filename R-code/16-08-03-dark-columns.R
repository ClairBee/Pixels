

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Image-plots/dark-columns/"

pw.m <- load.pixel.means()

####################################################################################################

# FIND ALL DARK LINES                                                                           ####

clump.lines <- function(px) {
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = c(2048, 2048))), dir = 4)
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))), 
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    df <- ddply(xy, .(id, x), summarise, 
                ymin = min(y), ymax = max(y), length = ymax - ymin + 1)
    df <- df[df$length > 200,]
    df
}

dpx <- apply(pw.m[,,"white",], 3, function(im) which (im < 20000, arr.ind = T))
dl <- lapply(dpx, clump.lines)

# sapply(dl, nrow)

invisible(lapply(names(dl)[-c(7:19)],
                 function(dt) {
                     apply(dl[[dt]], 1,
                           function(ll) {
                               vv <- pw.m[as.integer(ll["x"]),,, dt]
                               db <- pw.m[as.integer(ll["x"]), as.integer(ll["ymin"]):as.integer(ll["ymax"]), , dt]
                               neighb <- pw.m[as.integer(ll["x"])-1,,"black", dt]
                               neighb2 <- pw.m[as.integer(ll["x"])+1,,"black", dt]
                               sd.l <- apply(db, 2, sd, na.rm = T)
                               sd.neighb <- round(median(c(sd(neighb[1:1024], na.rm = T), sd(neighb[1024:2048], na.rm = T), sd(neighb2[1:1024], na.rm = T), sd(neighb2[1024:2048], na.rm = T))), 1)
                               xl <- c(0,1024) + (ll["ymin"] > 1024.5) * 1024
                               pdf(paste0(fpath, "dark-col-", dt, "-", ll["x"], ".pdf")); {
                                   par(mar = c(2,2,3,1))
                                   plot(vv[,"white"], type = "l", col = "gold", xlim = xl, ylim = range(pretty(db)), xlab = "", ylab = "",
                                        main = paste0(dt, " - ", ll["x"], ": ", paste(round(sd.l,1), collapse = ", "), " (", sd.neighb, ")"))
                                   lines(neighb, col = adjustcolor("cyan3", alpha = 0.4))
                                   lines(neighb2, col = adjustcolor("skyblue", alpha = 0.4))
                                   
                                   lines(vv[,"grey"], col = "green3")
                                   lines(vv[,"black"])
                                   dev.off()
                                   }
                               })
                     }))
                     
           
####################################################################################################

# PRODUCE SUMMARY STATISTICS FOR EACH COLUMN                                                    ####

# bgw sd along defect region
# sd along neighouring healthy line (2 columns away: checked, no overlap)
# end point: in visible region?
# difference between midline & endoipnt values
# difference between endpoint & 200px from endpoint

# will need to condense to individual lines: 131122, loan, MCT225
df <- invisible(sapply(c("130701", "140128", "loan"),
                       function(dt) {
                           rbind.fill(apply(dl[[dt]], 1,
                                 function(ll) {
                                     colm <- as.numeric(ll["x"])
                                     ymin <- as.numeric(ll["ymin"]); ymax <- as.numeric(ll["ymax"])
                                     
                                     db <- pw.m[colm, ymin:ymax, , dt]
                                     neighb <- pw.m[colm - 2, ymin:ymax, "black", dt]
                                     
                                     mid <- c(ymin, ymax)[which.min(abs(c(ymin, ymax)-1024.5))]
                                     end <- c(ymin, ymax)[which.max(abs(c(ymin, ymax)-1024.5))]
                                     fc <- end %in% range(which(!is.na(pw.m[1389,,"white",dt]), arr.ind = T))

                                     mid.v <- apply(pw.m[colm, c(mid + c(0:9) * ((end > 1024.5) * 2 - 1)), , dt], 2, median)
                                     end.v <- apply(pw.m[colm, c(end - c(0:9) * ((end > 1024.5) * 2 - 1)), , dt], 2, median)
                                     end200.v <- apply(pw.m[colm, c(end - c(100:109) * ((end > 1024.5) * 2 - 1)), , dt], 2, median)
                                     
                                     data.frame("acq" = dt, "col" = colm, "mid" = mid, "end" = end, "full" = fc, "black" = sd(db[,"black"]), 
                                                "grey" = sd(db[,"grey"]), "white" = sd(db[,"white"]), "neighb" = sd(neighb, na.rm = T),
                                                diff.b = end.v["black"] - mid.v["black"], diff.g = end.v["grey"] - mid.v["grey"], diff.w = end.v["white"] - mid.v["white"],
                                                diff200.b = end.v["black"] - end200.v["black"], diff200.g = end.v["grey"] - end200.v["grey"], diff200.w = end.v["white"] - end200.v["white"],
                                                stringsAsFactors = F)
                       }))}, simplify = F))

# run MCT225 separately because no midline

df$"MCT225" <- rbind.fill(apply(dl[["MCT225"]], 1,
                       function(ll) {
                           colm <- as.numeric(ll["x"])
                           ymin <- as.numeric(ll["ymin"]); ymax <- as.numeric(ll["ymax"])
                           
                           db <- pw.m[colm, ymin:ymax, , dt]
                           neighb <- pw.m[colm - 2, ymin:ymax, "black", dt]
                           
                           mid <- ymax
                           end <- ymin
                           fc <- end == 25
                           
                           mid.v <- apply(pw.m[colm, mid - c(0:9), , dt], 2, median)
                           end.v <- apply(pw.m[colm, end + c(0:9), , dt], 2, median)
                           end200.v <- apply(pw.m[colm, end - c(100:109), , dt], 2, median)
                           
                           data.frame("acq" = dt, "col" = colm, "mid" = mid, "end" = end, "full" = fc, "black" = sd(db[,"black"]), 
                                      "grey" = sd(db[,"grey"]), "white" = sd(db[,"white"]), "neighb" = sd(neighb, na.rm = T),
                                      diff.b = end.v["black"] - mid.v["black"], diff.g = end.v["grey"] - mid.v["grey"], diff.w = end.v["white"] - mid.v["white"],
                                      diff200.b = end.v["black"] - end200.v["black"], diff200.g = end.v["grey"] - end200.v["grey"], diff200.w = end.v["white"] - end200.v["white"],
                                      stringsAsFactors = F)
                           }))

df.a <- rbind.fill(df)
df.a[order(df.a$col, df.a$acq),]

plot(df.a$diff.b, df.a$diff200.b, pch = 20, col = c("black", "red")[df.a$full + 1])
abline(h = 0); abline(v = 0)

plot(df.a[,c("diff.b", "diff.g", "diff.w", "diff200.b", "diff200.g", "diff200.w")],
     pch = 20, col = c("black", "red")[df.a$full + 1], lower.panel = panel.cor)

####################################################################################################

# DARK ROW IN LOAN PANEL                                                                        ####

# row 77 - plots for 

jpeg("./01_Paper/fig/exploratory/row-defect-image.jpg"); {
    par(mar = c(2,2,1,1))
    pixel.image(pw.m[,,"black", "loan"], xlim = c(256,512), ylim = c(0,256))
    dev.off()
}

jpeg("./01_Paper/fig/exploratory/row-defect-transect.jpg", height = 240); {
    par(mar = c(2,2,1,1))
    plot(pw.m[,76,"black", "loan"], type = "l", col = adjustcolor("cyan3", alpha = 0.4), ylim = range(pw.m[,76:77,"black", "loan"], na.rm = T), xlab = "", ylab = "")
    lines(pw.m[,78,"black", "loan"], col = adjustcolor("green3", alpha = 0.4))
    lines(pw.m[,77,"black", "loan"])
    legend("bottomright", col = c("cyan3", "black", "green3"), legend = paste("Row ", c(76:78), sep = ""), lty = 1)
    dev.off()
}

jpeg("./01_Paper/fig/exploratory/row-defect-diffs.jpg", height = 240); {
    par(mar = c(2,2,1,1))
    plot(pw.m[,76,"black", "loan"] - pw.m[,75,"black", "loan"], ylim = c(-1000, 1000), type = "l", col = adjustcolor("cyan3", alpha = 0.4), xlab = "", ylab = "")
    lines(pw.m[,79,"black", "loan"] - pw.m[,78,"black", "loan"], col = adjustcolor("green3", alpha = 0.4))
    lines(pw.m[,77,"black", "loan"] - pw.m[,78,"black", "loan"])
    legend("topright", col = c("cyan3", "black", "green3"), lty = 1,
           legend = c("Row 76 - row 75", "Row 77 - row 76", "row 78 - row 79"))
    dev.off()
}

jpeg("./01_Paper/fig/exploratory/row-defect-powers.jpg", height = 240); {
    par(mar = c(2,2,1,1))
    plot(pw.m[,77,"white", "loan"] - pw.m[,76,"white", "loan"], type = "l", col = "gold", xlab = "", ylab = "")
    lines(pw.m[,77,"grey", "loan"] - pw.m[,78,"grey", "loan"], col = "green3")
    lines(pw.m[,77,"black", "loan"] - pw.m[,78,"black", "loan"])
    legend("topright", col = c("gold", "green3", "black"), lty = 1,
           legend = c("White", "Grey", "Black"))
    dev.off()
}

