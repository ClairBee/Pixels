
# read .xtek files to extract usage information

xtekct <- as.data.frame(t(as.matrix(read.csv("~/Documents/Pixels/Other-data/Old-data/Inside-Out_White_Check.xtekct", 
                                             header = T, sep = "=", stringsAsFactors = F))), stringsAsFactors = F)
    
####################################################################################################

# GET ALL FILES                                                                                 ####

list.xtec.files <- function(path.head, extensions = c(".xtec", ".xtekct")) {
    folders <- list.dirs(path.head, full.names = F, recursive = T)
    
    extensions <- paste("\\", extensions, sep = "")

    files <- list()
    for (ff in folders) {
        files[[ff]] <- unlist(lapply(extensions, 
                               function(x) list.files(paste0(path.head, ff), pattern = x, full.names = T)))
    }
    unname(unlist(files))
}

grab.xtec.data <- function(f.list, n = c(1:length(f.list))) {

    f.list <- f.list[n]
    
    xt.list <- lapply(f.list, 
                      function(ff) as.data.frame(t(as.matrix(read.csv(ff, 
                                                                      header = T, sep = "=", stringsAsFactors = F))), 
                                                 stringsAsFactors = F))
    rbind.fill(xt.list)
}

# run separately in case of crashes
fl <- list.xtec.files("./Other-data/Z/CT data/", extensions = c(".xtec", ".xtekct", ".xml"))

df <- grab.xtec.data(fl, n = c(3))


####################################################################################################

# AND NOW, IN PYTHON                                                                            ####

import os

for root, dirs, files in os.walk("./Other-data/Z/CT data/"):
    for file in files:
    if file.endswith(".xtec"):
    return(os.path.join(root, file))

