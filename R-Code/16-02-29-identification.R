
library("IO.Pixels"); library(raster)
load.images(150828, "black")
load.images(150828, "grey")
load.images(150828, "white")


# Various approaches to identification of dead pixels

#=============================================================================================
#
# approach based on manual, but without performing offset:
#   dead pixel: daily mean value is exactly 0
#   hot pixel: daily mean value is exactly 65535
#   underperforming bright pixel: value > 1.5x median
#   unerperforming dark pixel: value < 0.45x median
#   noise: sigma > 6x median sigma
#   globally non-uniform: value is > 2% from median value
#   locally non-uniform: value is > 1% from median of its 9x9 neighbours

get.bad.pixels <- function(data) {
    
    # create temporary list to hold output
    bp <- list()
    
    pw.mean <- pixelwise.mean(data)
    pw.sd <- pixelwise.sd(data)
    median.value <- median(data)
    median.sd <- median(pw.sd)
    
    
    # always use tryNULL when creating data frames, in case no rows are returned
    
    # order in which each type is added to df is order of priority when classifying
    # (eg. duplicates in later categories will be removed)
    
    # dead pixels: daily mean value is exactly 0
    bp$dead <- tryNULL(cbind(data.frame(which(pw.mean == 0.0, arr.ind = T)), 
                             type = "Dead"))
    
    # hot pixels: daily mean value is exactly 65535
    bp$hot <- tryNULL(cbind(data.frame(which(pw.mean == 65535.0, arr.ind = T)),
                            type = "Hot"))
    
    # underperforming bright pixels: value > 1.5x median
    bp$bright <- tryNULL(cbind(data.frame(unique(which(data > (1.5 * median.value), arr.ind = T)[,c(1:2)])),
                               type = "Too bright"))
    if (!is.null(bp$bright)) {colnames(bp$bright) <- c("row", "col", "type")}
    
    # underperforming dark pixels: value < 0.45x median
    bp$dim <- tryNULL(cbind(data.frame(unique(which(data < (0.45 * median.value), arr.ind = T)[,c(1:2)])),
                            type = "Too dark"))
    if (!is.null(bp$dim)) {colnames(bp$dim) <- c("row", "col", "type")}
    
    # noise: pixel sigma > 6x median sigma
    bp$noisy <- tryNULL(cbind(data.frame(which(pw.sd > (6 * median.sd), arr.ind = T)),
                              type = "Noisy"))
    
    # merge separate data frames into one
    df <- data.frame(do.call("rbind", bp))
    
    # remove duplicates and return df
    return(df[!duplicated(df[,c(1:2)]),])
}

system.time(bp <- get.bad.pixels(b.150828))     # 174.03 elapsed

table(bp$type)

pw.m.b.150828 <- pixelwise.mean(b.150828)

pixel.image(pw.m.b.150828, title = "Pixelwise mean value, black channel 15-08-28")
points(bp[,1], bp[,2], pch = 20, cex = 0.5, col = c("black", "white", "aquamarine", "grey", "darkgreen")[bp$type])


# how does bad pixel identification differ if we carry it out on a panel-by-panel basis?
p <- panel.edges()
panel.bp <- list()

for (i in 1:32) {
    j <- ceiling(i/16)
    k <- ((i-1) %% 16)+1
    
    panel <- b.150828[c(p$x[k]:(p$x[k+1]-1)),
                      c(p$y[j]:(p$y[j+1]-1)),]
    tmp <- get.bad.pixels(panel)
    panel.bp[[i]] <- rbind(panel = paste0(p$x[k], ":", (p$x[k+1]-1),", ", p$y[j],":",(p$y[j+1]-1)),
                           tmp)
}

panel.all <- do.call("rbind", panel.bp)

# save all as .csv to create neater tables later
write.csv(panel.all, "./Other-data/Panelwise-bad-pixels.csv", row.names = F)

bp$panel <- apply(cbind(cut(bp$col, p$y-1), "-", cut(bp$row, p$x-1)), 1, paste, collapse = "")
write.csv(bp$panel, "./Other-data/All-bad-pixels.csv", row.names = F)

#=============================================================================================
# local non-uniformity:
system.time({
local.median.9.9 <- as.matrix(focal(raster(pw.m.b), 
                                    w = matrix(c(rep(1,81)), nrow = 9), 
                                    fun = median, na.rm = T))
})      # 158.45 elapsed

tmp <- which(pw.m.b < (0.99 * local.median.9.9), arr.ind = T)
bad.pixels <- rbind(bad.pixels, 
                    data.frame(x = tmp[,1], y = tmp[,2], type = rep("Local low", nrow(tmp))))

tmp <- which(pw.m.b > (1.01 * local.median.9.9), arr.ind = T)
bad.pixels <- rbind(bad.pixels, 
                    data.frame(x = tmp[,1], y = tmp[,2], type = rep("Local high", nrow(tmp))))


# global non-uniformity:
# temporarily excluded: picks up 1/3 of all pixels each time under current method
# tmp <- which((pw.m.b > (1.02 * median.black)), arr.ind = T)       # 1341656
# tmp <- which((pw.m.b < (0.98 * median.black)), arr.ind = T)       # 1301899


# count of pixel types before removing duplicated
# Dead     Hot    Too bright   Too dark      Noisy  Local low   Local high 
# 5        125       1130          5         1319     821413      810268 

# remove duplicates from bad pixel list - should keep first instance as long as data isn't sorted
bad.pixels <- bad.pixels[!duplicated(bad.pixels[,c(1:2)]),]


# Dead     Hot    Too bright   Too dark      Noisy  Local low   Local high 
# 5        125       1005          0          788     821407       808417 