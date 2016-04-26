
####################################################################################################
# WRITE UP IDENTIFICATION OF NON-RESPONSIVE PIXELS
# ADD NON-RESPONSIVE PIXELS TO SUPERCLUSTERS

library("IO.Pixels")

fpath <- "./Notes/Superclusters/fig/"
load.pixel.means()

bp <- setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                   rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                          factor("edge", levels = c("edge", "no.response", "dead", "hot", "dim spot", "v.bright", "bright"))),
               c("x", "y", "type"))

####################################################################################################
# IDENTIFY NON-RESPONSIVE PIXELS                                                                ####

# back-to-back histogram of black & grey pixels, showing limits of 'normal' behaviour
{
    hb <- hist(pw.m[,,"black", "160314"], breaks = "fd", plot = F)
    
    hg <- hist(pw.m[,,"grey", "160314"], breaks =  "fd", plot = F)
    hg.c <- hist(pw.m[11:1985, 11:1985, "grey", "160314"], breaks = "fd", plot = F)
    hb$counts <- -hb$counts
    
    hw <- hist(pw.m[,,"white", "160314"], breaks =  "fd", plot = F)
    hw.c <- hist(pw.m[11:1985, 11:1985, "white", "160314"], breaks = "fd", plot = F)

    JF <- JohnsonFit(pw.m[,, "black", "160314"])
    
    pdf(paste0(fpath, "grey-hist-no-response.pdf"), height = 4, width = 7)
    {
        par(mar = c(2,2,1,1))
        plot(hg, ylim = c(-30, 30), col = "lightseagreen", border = "lightseagreen", xlab = "", ylab = "", main = "")
        plot(hb, add = T)
        rect(qJohnson(0.005, JF), -50, qJohnson(0.995, JF), 50, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
        plot(hg.c, add = T, col = "blue", border = "blue")
        legend("bottomright", pch = 15, bty = "n", cex = 0.8, inset = 0.05,
               col = c("lightseagreen", "blue", "black", adjustcolor("skyblue", alpha = 0.3)),
               legend = c("All grey pixels", "Grey pixels, cropped", "Black pixels", "Normal black response"))
    }
    dev.off()
    
    
    pdf(paste0(fpath, "white-hist-no-response.pdf"), height = 4, width = 7)
    {
        par(mar = c(2,2,1,1))
        plot(hw, ylim = c(-30, 30), col = "lightseagreen", border = "lightseagreen", xlab = "", ylab = "", main = "")
        plot(hb, add = T)
        rect(qJohnson(0.005, JF), -50, qJohnson(0.995, JF), 50, col = adjustcolor("skyblue", alpha = 0.3), border = NA)
        plot(hw.c, add = T, col = "blue", border = "blue")
        legend("bottomright", pch = 15, bty = "n", cex = 0.8, inset = 0.05,
               col = c("lightseagreen", "blue", "black", adjustcolor("skyblue", alpha = 0.3)),
               legend = c("All white pixels", "White pixels, cropped", "Black pixels", "Normal black response"))
    }
    dev.off()
    
    
}

bp <- rbind(bp,
            setNames(data.frame(no.response(160314), type = "no.response"), c("x", "y", "type")))

table(bp$type)

# IDENTIFY HOT PIXELS                                                                           ####

# IDENTIFY DEAD PIXELS                                                                          ####

# IDENTIFY BRIGHT PIXELS                                                                        ####

# shades of bright/dim pixels: use midpoints as well as/instead of residuals?

# FIND SUPERCLUSTERS OF BRIGHT, HOT, NONRESPONSIVE                                              ####
