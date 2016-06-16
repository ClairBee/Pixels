
# EXTRACTION CODE                                                                               ####

ctprofile.list <- function(folder.list) {
    unlist(lapply(folder.list, function(ff) list.files(ff, pattern = "\\.ctprofile", full.names = T)))
}

get.ctprofiles <- function(ctp.list) {
    require(plyr); require(XML)
    xt.list <- lapply(ctp.list,
                      function(ctp) as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(ctp)),
                                                                        list(recursive=TRUE))), FUN = unlist),
                                                  stringsAsFactors = F))
    xt.list <- cbind(rbind.fill(xt.list), filenm = ctp.list)
}

subfolders <- function(folder.list) {
    unlist(lapply(folder.list, list.dirs, full.names = T, recursive = F))
}

subfolder.ctprofiles <- function(path.head, df) {
    ff <- list.dirs(path.head, recursive = F)
    
    ct.list <- as.character(df$filenm)
    ct.list <- paste(unlist(lapply(lapply(lapply(ct.list,
                                                 strsplit, split = "]"), "[[", 1), "[[", 1)), "]", sep = "")
    get.ctprofiles(ctprofile.list(subfolders(ff[!(ff %in% ct.list)])))
}

# get CT profiles from top level
{
    df.y1 <- get.ctprofiles(ctprofile.list(list.dirs("Y:/CT data", recursive = F)))
    df.y2 <- get.ctprofiles(ctprofile.list(list.dirs("Y:/CT data exceptions", recursive = F)))
    df.x <- get.ctprofiles(ctprofile.list(list.dirs("X:/Nikon data", recursive = F)))

    found <- rbind.fill(list(df.y1, df.y2, df.x))
    write.csv(found, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles.csv")
}

# get CT profiles from first level of subfolders
{
    df2.y1 <- subfolder.ctprofiles("Y:/CT data", found)
    df2.y2 <- subfolder.ctprofiles("Y:/CT data exceptions", found)
    df2.x <- subfolder.ctprofiles("X:/Nikon data", found)
    
    found.sub <- rbind.fill(list(df2.y1, df2.y2, df2.x))
    write.csv(found.sub, "C:/Users/coe-admin/Documents/R/CT profiles/CT-profiles.csv")
}

####################################################################################################

# SUMMARISE DETECTOR USAGE                                                                      ####

xtek <- rbind.fill(lapply(list.files("./Other-data/Z/", pattern = "CT-profiles.+csv", full.names = T),
             read.csv, row.names = 1, stringsAsFactors = F))

write.csv(xtek, "./Other-data/xtek-usage.csv")

xt <- data.frame(acq.date = strptime(lapply(sapply(lapply(sapply(xtek$filenm, strsplit, "[[]"), 
                                                               "[", 2), strsplit, "[]]"), "[", 1),
                                     "%Y-%m-%d %H.%M.%S"),
                 kV = xtek$XraySettings.Settings.kV,
                 uA = xtek$XraySettings.Settings.uA,
                 exp.time = xtek$ImagingSettings..attrs.exposure,
                 n.proj = xtek$Projections,
                 frames.per.proj = xtek$FramesPerProjection,
                 filter.mat = xtek$XrayFilterMaterial,
                 filter.thk = xtek$XrayFilterThickness,
                 filenm = xtek$filenm,
                 stringsAsFactors = F)

####################################################################################################

# EMPTY FOLDERS                                                                                 ####

mm <- unlist(read.csv("./Other-data/Z/Folders-without-ctprofile.csv", row.names = 1, stringsAsF = F))

acq <- strptime(lapply(sapply(lapply(sapply(mm, strsplit, "[[]"), 
                              "[", 2), strsplit, "[]]"), "[", 1),
         "%Y-%m-%d %H.%M.%S")

dd <- as.Date(acq)

sort(dd)
