

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Image-plots/dark-columns/"

pw.m <- load.pixel.means()

####################################################################################################

# find all dark lines

clump.lines <- function(px) {
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = c(2048, 2048))), dir = 4)
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))), 
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    df <- ddply(xy, .(id, x), summarise, 
                ymin = min(y), ymax = max(y), length = ymax - ymin + 1)
    df <- df[df$length > 10,]
    df
}

os <- pw.m[,,"white",] - pw.m[,,"black",]
dpx <- apply(os, 3, function(im) which(im < 15000, arr.ind = T))
dl <- lapply(dpx, clump.lines)

sapply(dl, nrow)

invisible(lapply(names(dl)[-c(7:19)],
                 function(dt) {
                     apply(dl[[dt]], 1,
                           function(ll) {
                               vv <- pw.m[as.integer(ll["x"]),(1:1024) + (ll["ymin"] > 1024.5) * 1024,, dt]
                               db <- pw.m[as.integer(ll["x"]), as.integer(ll["ymin"]):as.integer(ll["ymax"]), , dt]
                               neighb <- pw.m[as.integer(ll["x"])-1,(1:1024) + (ll["ymin"] > 1024.5) * 1024,"black", dt]
                               neighb2 <- pw.m[as.integer(ll["x"])+1,(1:1024) + (ll["ymin"] > 1024.5) * 1024,"black", dt]
                               sd.l <- apply(db, 2, sd, na.rm = T)
                               pdf(paste0(fpath, "dark-col-", dt, "-", ll["x"], ".pdf")); {
                                   par(mar = c(2,2,3,1))
                                   plot(vv[,"white"], type = "l", col = "gold", ylim = range(pretty(db)), xlab = "", ylab = "",
                                        main = paste0(dt, " - ", ll["x"], ": ", paste(round(sd.l,1), collapse = ", ")))
                                   lines(neighb, col = adjustcolor("cyan3", alpha = 0.4))
                                   lines(neighb2, col = adjustcolor("skyblue", alpha = 0.4))
                                   
                                   lines(vv[,"grey"], col = "green3")
                                   lines(vv[,"black"])
                                   dev.off()
                                   }
                               })
                     }))
                     
           