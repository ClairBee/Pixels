
####################################################################################################

# DATA FLOW IN XRD-1621 PANEL                                                                   ####

library(shape)

pdf("./01_Paper/fig/intro/data-readout.pdf"); {
    plot(0, type = "n", xlim = c(0,2048), ylim = c(0,2048), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    draw.panels()
    rect(0,0,2048,2048)
    
    text(x = rep(c(128 * 1:16) - 56, 2),
         y = c(rep(1024+512, 16), rep(1024-512, 16)),
         labels = c(1:32), cex = 0.6)
    
    b <- 20
    arr.col <- "cyan3"
    
    Arrows(x0 = c(128 * 1:16) - b, y0 = 1024 + b,
           x1 = c(128 * 1:16) - b, y1 = 2048 + b,
           col = arr.col, code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    Arrows(x0 = c(128 * 1:16) - b, y0 = 2048 + b,
           x1 = c(128 * 0:15) + b, y1 = 2048 + b,
           col = arr.col, code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    
    Arrows(x0 = c(128 * 0:15) + b, y0 = 1024 - b,
           x1 = c(128 * 0:15) + b, y1 = -b,
           col = arr.col, code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    Arrows(x0 = c(128 * 0:15) + b, y0 = -b,
           x1 = c(128 * 1:16) - b, y1 = -b,
           col = arr.col, code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    
    Arrows(x0 = 2048, y0 = 2048 + b * 2.5, x1 = 0, y1 = 2048 + b * 2.5,
           col = "green3", code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    Arrows(x0 = 2048, y0 = - b * 2.5, x1 = 0, y1 = - b * 2.5,
           col = "green3", code = 2, arr.type = "curved", arr.length = 0.2, arr.width = 0.1, arr.adj = 1)
    
    points(x = c(c(128 * 0:15) + b, c(128 * 1:16) - b), y = c(rep(2048 + b, 16), rep(-b, 16)),
           pch = 15, col = "red", cex = 0.6)
    
    dev.off()
}

crop.pdf("./01_Paper/fig/intro/data-readout.pdf")



####################################################################################################

# SUMMARY TABLE OF ALL IMAGES                                                                   ####

library("IO.Pixels")
# get xml data for new images
new.xml <- summarise.profiles()

df <- setNames(new.xml[new.xml$batch == "grey", c("date", "ImageNotes", "kV", "uA", "ExposureTime", "ImageDetector", "XraySource", "XraySourceTargetMaterial")],
               nm = c("code", "timestamp", "g.kV", "g.uA", "g.exp", "model", "source", "target"))
df$model[is.na(df$model)] <- df$model[df$code == "160430"]
df$source[is.na(df$source)] <- df$source[df$code == "160430"]
df$target[is.na(df$target)] <- df$target[df$code == "160430"]

df$w.kV <- new.xml$kV[new.xml$batch == "white"]
df$w.uA <- new.xml$uA[new.xml$batch == "white"]
df$w.exp <- new.xml$Exposure[new.xml$batch == "white"]

# get xml data for old images

old.dts <- c("130613", "130701", "131002", "131122", "140128", "140129")
old.xml <- rbind.fill(lapply(paste("./Other-data/Old-data/", old.dts, "/", sep = ""),
                  function(dt) {
                      as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(paste0(dt, "CalibrationParameters.xml"))),
                                                          list(recursive = TRUE))), 
                                           FUN = unlist), stringsAsFactors = F)
                  }))
zz <- data.frame(code = old.dts,
                 timestamp = NA,
                 g.kV = NA,
                 g.uA = NA,
                 g.exp = NA,
                 model = df$model[df$code == "160430"],
                 source = df$source[df$code == "160430"],
                 target = df$target[df$code == "160430"],
                 w.kV = as.numeric(old.xml$kV),
                 w.uA = as.numeric(old.xml$uA),
                 w.exp = as.numeric(old.xml$ImagingSettings..attrs.exposure),
                 stringsAsFactors = F)

# manually add timestamps based on tech report
zz$timestamp = c("13 June 2013 13:31:51", "01 July 2013 11:49:29", "02 October 2013 13:41:00",
                 "22 November 2013 10:54:30", "28 January 2014 11:48:00", "28 January 2014 15:14:02")

# combine into data frame & output to .csv for import into LaTeX
xml.summ <- rbind(zz, df)
rownames(xml.summ) <- xml.summ[,1]

# convert timestamps to POSIX format
xml.summ$timestamp <- strptime(xml.summ$timestamp, format = "%d %B %Y %H:%M:%S")

# add columns denoting history & panel ID
xml.summ$notes <- ""
xml.summ[c("140128", "141009"), "notes"] <- "Refurbished"
xml.summ[c("160314"), "notes"] <- "New software"

xml.summ$panel.id <- "WMG"
xml.summ["loan", "panel.id"] <- "WMG loan"
xml.summ["MCT225", "panel.id"] <- "Nikon"

# calculate power/exp usage
xml.summ$g.power <- round(xml.summ$g.kV * xml.summ$g.uA * xml.summ$g.exp / 1000 / 1000, 1)
xml.summ$w.power <- round(xml.summ$w.kV * xml.summ$w.uA * xml.summ$w.exp / 1000 / 1000, 1)

xml.summ[is.na(xml.summ)] <- "-"

# shorten detector names
xml.summ$model <- gsub("PerkinElmer", "", xml.summ$model)
# remove fullstops from column names, for transfer to LaTeX
colnames(xml.summ) <- gsub("\\.", "", colnames(xml.summ))

# flag to add horizontal lines in table
xml.summ$hlflag <- ""
xml.summ[c("141009", "loan", "MCT225"), "hlflag"] <- "x"

xml.summ$head <- ""
xml.summ["140129", "head"] <- "\\multirow{-6}{1em}{\\rotatebox{90}{\\centering{Pilot data}}}"
xml.summ["160705", "head"] <- "\\multirow{-14}{1em}{\\rotatebox{90}{\\centering{Main data}}}"

write.csv(xml.summ, "./01_Paper/fig/tech-data/Detector-details.csv", quote = F, row.names = F)

xml.summ <- read.csv("./01_Paper/fig/tech-data/Detector-details.csv", as.is = T)
row.names(xml.summ) <- xml.summ[,1]

