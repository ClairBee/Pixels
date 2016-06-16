
# V1 ####

list.xtec.files <- function(path.head, extensions = c(".ctprofile")) {
  
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
fl <- list.xtec.files("Y:/CT data/")
df <- grab.xtec.data(fl, n = c(1:10))

length(list.dirs("Y:/CT data/", recursive = F))

##############################################################################################

# FINAL ####

# get list of ct profiles found
ctprofile.list <- function(folder.list) {
    unlist(lapply(folder.list, function(ff) list.files(ff, pattern = "\\.ctprofile", full.names = T)))
}


# read ct profiles into single data frame
get.ctprofiles <- function(ctp.list) {
  require(plyr); require(XML)
  xt.list <- lapply(ctp.list,
                    function(ctp) as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(ctp)),
                                                                      list(recursive=TRUE))), FUN = unlist),
                                                stringsAsFactors = F))
  xt.list <- cbind(rbind.fill(xt.list), filenm = ctp.list)
}

#---------------------------------------------------------------------------------------------

# run for Y:/CT data
{
  ff <- list.dirs("Y:/CT data", recursive = F)
  
  ct.list <- ctprofile.list(ff)
  
  zz <- get.ctprofiles(ct.list)
  
  write.csv(zz, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-1.csv")
}

#---------------------------------------------------------------------------------------------

# run for Y:/CT data exceptions
{
  ff <- list.dirs("Y:/CT data exceptions", recursive = F)
  
  ct.list <- ctprofile.list(ff)
  
  zz <- get.ctprofiles(ct.list)
  
  write.csv(zz, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-2.csv")
}

#---------------------------------------------------------------------------------------------

# run for X:/Nikon data
{
  ff <- list.dirs("X:/Nikon data", recursive = F)
  
  ct.list <- ctprofile.list(ff)
  
  zz <- get.ctprofiles(ct.list)
  
  write.csv(zz, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-3.csv")
}

#---------------------------------------------------------------------------------------------



# FIND ANY MISSED PROFILES ####

subfolders <- function(folder.list) {
  unlist(lapply(folder.list, list.dirs, full.names = T, recursive = F))
}

# run for Y:/CT data
{
  ff <- list.dirs("Y:/CT data", recursive = F)

  ct.list <- as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-1.csv")$filenm)
  ct.list <- paste(unlist(lapply(lapply(lapply(ct.list,
                                               strsplit, split = "]"), "[[", 1), "[[", 1)), "]", sep = "")
  
  # get CT profiles
  zz <- get.ctprofiles(ctprofile.list(subfolders(ff[!(ff %in% ct.list)])))
  
  write.csv(zz, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-1-2.csv")
}

# check subfolders
subfolder.ctprofiles <- function(path.head, df) {
    ff <- list.dirs(path.head, recursive = F)
    
    ct.list <- as.character(df$filenm)
    ct.list <- paste(unlist(lapply(lapply(lapply(ct.list,
                                                 strsplit, split = "]"), "[[", 1), "[[", 1)), "]", sep = "")
    get.ctprofiles(ctprofile.list(subfolders(ff[!(ff %in% ct.list)])))
}

# nothing found in Y:/CT data exceptions
{
  ff <- list.dirs("Y:/CT data exceptions", recursive = F)
  
  ct.list <- as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-2.csv")$filenm)
  ct.list <- paste(unlist(lapply(lapply(lapply(ct.list,
                                               strsplit, split = "]"), "[[", 1), "[[", 1)), "]", sep = "")
  
  # get CT profiles
  zz <- get.ctprofiles(ctprofile.list(subfolders(ff[!(ff %in% ct.list)])))
}

# nothing found in T:/CT data exceptions
{
  ff <- list.dirs("X:/Nikon data", recursive = F)
  
  ct.list <- as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-3.csv")$filenm)
  ct.list <- paste(unlist(lapply(lapply(lapply(ct.list,
                                               strsplit, split = "]"), "[[", 1), "[[", 1)), "]", sep = "")
  
  # get CT profiles
  zz <- get.ctprofiles(ctprofile.list(subfolders(ff[!(ff %in% ct.list)])))
  write.csv(zz, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-3-2.csv")
}

##############################################################################################

# LIST ALL UNMATCHED FILES                                                                ####

ff <- c(list.dirs("Y:/CT data", recursive = F),
        list.dirs("Y:/CT data exceptions", recursive = F),
        list.dirs("X:/Nikon data", recursive = F))

hh <- c(as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-1.csv")$filenm),
        as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-2.csv")$filenm),
        as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-3.csv")$filenm),
        as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-1-2.csv")$filenm),
        as.character(read.csv("C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles-3-2.csv")$filenm))

hh <- unique(lapply(lapply(strsplit(hh, "/"), "[", 1:3), function(s) paste(unlist(s), collapse = "/")))

qq <- ff[!(ff %in% hh)]
write.csv(qq, "C:/Users/coe-admin/Documents/R/CT profiles/Folders-without-ctprofile.csv")


length(unique(pp))
