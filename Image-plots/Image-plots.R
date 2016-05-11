
library("IO.Pixels"); library("CB.Misc")

load.pixel.means()

####################################################################################################

# DIFFERENCES BETWEEN SUCCESSIVE ACQUISITIONS                                                   ####

diff <- pw.m[,,, 2:12] - pw.m[,,, 1:11]

for (col in dimnames(diff)[[3]]) {
    for (dt in dimnames(diff)[[4]]) {
        pdf(paste0("./Diffs-between-acquisitions/image-diffs-", col, "-", dt, ".pdf"))
            pixel.image(diff[ , , col, dt], title = paste0(col, " - ", dt))
        dev.off()
    }
}

####################################################################################################

# RAW IMAGES                                                                                    ####

