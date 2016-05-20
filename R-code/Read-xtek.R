
# read .xtek files to extract usage information

xtekct <- as.data.frame(t(as.matrix(read.csv("~/Documents/Pixels/Other-data/Old-data/Inside-Out_White_Check.xtekct", 
                                             header = T, sep = "=", stringsAsFactors = F))), stringsAsFactors = F)
    
    