
add.scale <- function(z, zlim = range(z, na.rm = T), col = sd.colours(), breaks = sd.levels(z),
                      horiz = TRUE, ylim = NULL, xlim = NULL, ...) {
    if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
    }
    if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    poly <- vector(mode="list", length(col))
    for(i in seq(poly)){
        poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
    }
    xaxt <- ifelse(horiz, "s", "n")
    yaxt <- ifelse(horiz, "n", "s")
    if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
    if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
    if(missing(xlim)) xlim=XLIM
    if(missing(ylim)) ylim=YLIM
    plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
    for(i in seq(poly)){
        if(horiz){
            polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
        }
        if(!horiz){
            polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
        }
    }
}


layout(matrix(c(1,2), nrow = 2, ncol = 1), heights = c(2,1), widths = c(1,1))
# layout.show(2)

par(mar = c(2,2,3,1))
pixel.image(sc[,,"loan"], title = "Shading-corrected: loan")

par(mar = c(3,1,1))
add.scale(sc[,,"loan"], horiz = T)

image.scale(sc[,,"loan"], axis.pos = 1)
plot.new()
dev.off()
