
# convert all columns of one type to another
DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)