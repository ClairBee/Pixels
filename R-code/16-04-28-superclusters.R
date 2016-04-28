
library("IO.Pixels")

fpath <- "./Notes/Superclusters/fig/"
load.pixel.means()
res <- readRDS("./Models/Simple-parametric/Residuals-simple-parametric.rds")
bp <- readRDS("./Models/Simple-parametric/Bad-px-new-thresholds.rds")

####################################################################################################



####################################################################################################
# GET BAD PIXELS ACROSS ALL IMAGES                                                              ####
# approx. 8 minutes to run across 11 image sets
    {
        # bp <- lapply(dimnames(pw.m)[[4]], bad.px)
        # saveRDS(bp, "./Models/Simple-parametric/Bad-px-new-thresholds.rds")
    }

####################################################################################################

# CLUSTERS & SUPERCLUSTERS                                                                      ####
bpm <- enhance.spots(bp$"160314")


# REVISE THRESHOLDS: ADD ANOTHER STRIP AT INTERIOR?                                             ####