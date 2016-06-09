
# COMPARISON OF PROPOSED CLASSIFICATION METHODS & SELECTION OF FINAL APPROACH

# considering observable state spaces only

library("IO.Pixels"); library("CB.Misc")
fpath <- "./Notes/Final-classifications/fig/"

load.pixel.means()

# median-differenced images to find local non-uniformity
md.b <- readRDS("./Other-data/Median-diffs-black.rds")
md.g <- readRDS("./Other-data/Median-diffs-grey.rds")

####################################################################################################

# MEAN-VALUED CLASSIFICATION                                                                    ####
    # this is the 'standard' classification approach used throughout the project so far.
    # simplest approah - use as baseline against which to look for improvements.

# define possible bad pixel states 
Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "l.bright", "screen spot", "line.d", "edge", "v.dim", "dim", "l.dim")
Cat.cols <- c("purple", "black", "magenta3", "red", "orange", "gold", "yellow", "grey", "violet", NA, "green3", "lightskyblue", "grey")


# classify bad pixels
{
    bp <- apply(pw.m, 4,
                function(acq) rbind(data.frame(edge.px(acq), type = ordered("edge", levels = Cat)),
                                    data.frame(no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]), type = "no.resp"),
                                    data.frame(which(acq[,,"black"] == 65535, arr.ind = T), type = "hot"),
                                    data.frame(which(acq[,,"white"] == 0, arr.ind = T), type = "dead"),
                                        data.frame(which(acq[,,"white"] > median(acq[,,"white"]) + (65535 - median(acq[,,"white"])) / 2, arr.ind = T), type = "v.bright"),
                                        data.frame(which(acq[,,"white"] > median(acq[,,"white"]) + (65535 - median(acq[,,"white"])) / 4, arr.ind = T), type = "bright"),
                                        data.frame(which(acq[,,"white"] < median(acq[,,"white"]) / 2, arr.ind = T), type = "v.dim"),
                                        data.frame(which(acq[,,"white"] < median(acq[,,"white"]) * 0.75, arr.ind = T), type = "dim"),
                                    data.frame(which(acq[,,"grey"] > median(acq[,,"grey"]) + (65535 - median(acq[,,"grey"])) / 2, arr.ind = T), type = "v.bright"),
                                    data.frame(which(acq[,,"grey"] > median(acq[,,"grey"]) + (65535 - median(acq[,,"grey"])) / 4, arr.ind = T), type = "bright"),
                                    data.frame(which(acq[,,"grey"] < median(acq[,,"grey"]) / 2, arr.ind = T), type = "v.dim"),
                                    data.frame(which(acq[,,"grey"] < median(acq[,,"grey"]) * 0.75, arr.ind = T), type = "dim"),
                                        data.frame(which(acq[,,"black"] > median(acq[,,"black"]) + (65535 - median(acq[,,"black"])) / 2, arr.ind = T), type = "v.bright"),
                                        data.frame(which(acq[,,"black"] > median(acq[,,"black"]) + (65535 - median(acq[,,"black"])) / 4, arr.ind = T), type = "bright"),
                                        data.frame(which(acq[,,"black"] < median(acq[,,"black"]) / 2, arr.ind = T), type = "v.dim"),
                                        data.frame(which(acq[,,"black"] < median(acq[,,"black"]) * 0.75, arr.ind = T), type = "dim"),
                                    data.frame(which(find.lines(acq[, , "black"]) > 0, arr.ind = T), type = "line.b"),
                                    data.frame(which(find.lines(acq[, , "black"], dim.lines = T) > 0, arr.ind = T), type = "line.d"),
                                    data.frame(which(find.lines(acq[, , "black"], horizontal = T) > 0, arr.ind = T), type = "line.b"),
                                    data.frame(which(find.lines(acq[, , "black"], dim.lines = T, horizontal = T) > 0, arr.ind = T), type = "line.d"),
                                        data.frame(which(threshold(md.b[[x]], level = mad(pw.m[,,"black", x]) * 2) > 0, arr.ind = T), type = "l.bright")))
    
    bp <- apply(pw.m, 4,
                function(img) edge.px(img))                 
    # screen spots
    # locally dim
    # locally bright
    qq <- apply(pw.m, 4, apply, 3, 
                function(img) data.frame(which(img > median(img) + (65535 - median(img)) / 2, arr.ind = T)))
    
    qq <- data.frame(unique(rbind.fill(apply(pw.m[,,,"160430"], 3, function(img) data.frame(which(img > median(img) + (65535 - median(img)) / 2, arr.ind = T))))),
                     type = "v.bright")
    
    data.frame(which(acq[,,"white"] > median(acq[,,"white"]) + (65535 - median(acq[,,"white"])) / 2, arr.ind = T), type = "v.bright"),
    data.frame(which(acq[,,"white"] > median(acq[,,"white"]) + (65535 - median(acq[,,"white"])) / 4, arr.ind = T), type = "bright"),
    data.frame(which(acq[,,"white"] < median(acq[,,"white"]) / 2, arr.ind = T), type = "v.dim"),
    data.frame(which(acq[,,"white"] < median(acq[,,"white"]) * 0.75, arr.ind = T), type = "dim"),
}

# plot transitions as function of real time - are rates constant?
{}

# check stability of classes
{
    # mean transition rates with 95% confidence intervals?
}

# cluster pixels according to path taken?
{
    # filter out all pixels that are always healthy
    
    
    # plot trajectory through state space
    
    
    # cluster manually or by distance metric?
    
    
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
