
library("IO.Pixels")

############################################################################################################
# Table of acquisitions                                                                                 ####

# extract as much as possible from image xml
df <- rbind(rbind.fill(lapply(list.dirs("./Other-data/Old-data", recursive = F, full.names = F)[1:6], 
                              # iterate over pilot set
                              function(fnm) {
                                  c.params <- xmlToList(xmlParse(paste0("./Other-data/Old-data/", fnm, "/CalibrationParameters.xml")))
                                  
                                  data.frame("code" = fnm, 
                                             "timestamp" = NA,
                                             "model" = "",
                                             "source" = "",
                                             "target" = "",
                                             "gkV" = NA,
                                             "guA" = NA,
                                             "gexp" = NA,
                                             "wkV" = as.numeric(c.params$kV),
                                             "wuA" = as.numeric(c.params$uA),
                                             "wexp" = as.numeric(c.params$ImagingSettings$.attrs["exposure"]),
                                             stringsAsFactors = F)
                              })),
            rbind.fill(lapply(list.dirs("./Image-data", recursive = F, full.names = F),
                              # iterate over experimental sets
                              function(fnm) {
                                  g.profile <- xmlToList(xmlParse(list.files(paste0("./Image-data/", fnm, "/grey"),
                                                                             pattern = "\\.xml", full.names = T)[1]))
                                  w.profile <- xmlToList(xmlParse(list.files(paste0("./Image-data/", fnm, "/white"),
                                                                             pattern = "\\.xml", full.names = T)[1]))
                                  
                                  data.frame("code" = fnm, 
                                             "timestamp" = paste(as.Date(w.profile$TimeAcquired), w.profile$Time),
                                             "model" = gsub("PerkinElmer", "", toString(w.profile$ImageDetector)),
                                             "source" = toString(w.profile$XraySource),
                                             "target" = toString(w.profile$XraySourceTargetMaterial),
                                             "gkV" = as.numeric(g.profile$kV),
                                             "guA" = as.numeric(g.profile$uA),
                                             "gexp" = as.numeric(g.profile$ExposureTime),
                                             "wkV" = as.numeric(w.profile$kV),
                                             "wuA" = as.numeric(w.profile$uA),
                                             "wexp" = as.numeric(w.profile$ExposureTime),
                                             stringsAsFactors = F)
                              })))
row.names(df) <- df$code

# manually add pilot timestamps based on tech report
df$timestamp[1:6] <- c("2013-06-13 13:31:51", "2013-07-01 11:49:29", "2013-10-02 13:41:00",
                       "2013-11-22 10:54:30", "2014-01-28 11:48:00", "2014-01-28 15:14:02")

# calculate relative grey & white power input
df$gpower <- round(df$gkV * df$guA * df$gexp / 1000 / 1000, 2)
df$wpower <- round(df$wkV * df$wuA * df$wexp / 1000 / 1000, 2)

df[is.na(df)] <- "-"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add formatting for LaTeX csvreader...

# flag to add horizontal lines in table
df$hlflag <- ""
df[c("141009", "loan", "MCT225"), "hlflag"] <- "x"

# describe data streams (with multirow formatting)
df$head <- ""
df["140129", "head"] <- "\\multirow{-6}{1em}{\\rotatebox{90}{\\centering{Pilot data}}}"
df["160705", "head"] <- "\\multirow{-14}{1em}{\\rotatebox{90}{\\centering{Main data}}}"

# panel ID
df$panel <- ""
df["140129", "panel"] <- "\\multirow{-6}{1em}{\\centering{A}}"
df["160705", "panel"] <- "\\multirow{-14}{1em}{\\centering{A}}"
df["loan3", "panel"] <- "\\multirow{-3}{1em}{\\centering{B}}"
df["MCT225", "panel"] <- "\\multirow{-1}{1em}{\\centering{C}}"

# model as multirow
df$mmodel <- ""
df["140129", "mmodel"] <- paste0("\\multirow{-6}{1em}{\\centering{", df["160705", "model"],"}}")
df["160705", "mmodel"] <- paste0("\\multirow{-14}{1em}{\\centering{", df["160705", "model"],"}}")
df["loan2", "mmodel"] <- paste0("\\centering{", df["loan3", "model"], "}")
df["MCT225", "mmodel"] <- paste0("\\centering{", df["MCT225", "model"], "}")

# notes (on refurbishments, essentially)
df$notes <- ""
df[c("140128", "141009"), "notes"] <- "Refurbished"
df[c("160314"), "notes"] <- "New software"
df[c("loan"), "notes"] <- "Loan panel"
df[c("MCT225"), "notes"] <- "Loan panel"

# relabel loan panels with datestamp
rl <- c("loan", "loan2", "loan3", "MCT225")
df[rl, "code"] <- gsub("-", "", sapply(df[rl, "timestamp"], substr, 3,10))

# remove row for 15-07-02
df <- df[df$code != "150702",]

write.csv(df, "./01_Paper/fig/tech-data/Detector-details.csv", quote = F, row.names = F)


############################################################################################################
# Charge transfer diagram                                                                               ####
