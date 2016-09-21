
library("IO.Pixels"); library("CB.Misc")
pw.m <- load.objects("./02_Objects/images/", otype = "pwm")
md7 <- load.objects("./02_Objects/med-diffs/", otype = "md7")

fpath <- "./Image-plots/linear-res/"

.smoothScatter <- hijack(smoothScatter, nrpoints = 0, 
                         colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                         xlab = "", ylab = "")

####################################################################################################

# FUNCTIONS                                                                                     ####

quick.lm <- function(im, terms, midline = 1024.5, crop.data = T) {
    
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    
    if(!is.na(midline)) {
        df$upper <- df$y > midline
    }
    
    # fit linear model to central part of image only (excludes edge effects)
    if (crop.data) {
        df <- df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                             findInterval(df$y, c(40.5, 2008.5)) == 1, ]
    }
    
    lm(as.formula(terms), data = df)
}

lm.res <- function(im, terms, midline = 1024.5, crop.data = T) {
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    
    if(!is.na(midline)) {
        df$upper <- df$y > midline
    }
    
    # fit linear model to central part of image only (excludes edge effects)
    if (crop.data) {
        df <- df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                     findInterval(df$y, c(40.5, 2008.5)) == 1, ]
    }
    
    mod <- lm(as.formula(terms), data = df)
    
    fv <- predict(mod, df)
    array(df$g - fv, dim = dim(im[,,1]))
}

####################################################################################################

# FIT ALL MODELS                                                                                ####

# no interaction, data cropped
{
l.models <- invisible(apply(pw.m[,,,], 4, 
                            function(im) {
                                mm <- quick.lm(im, terms = "g ~ upper + (b + w)", crop.data = T)
                                c(c(coef(mm), "b:w" = 0)[1:5], rmse = summary(mm)$sigma, r2 = summary(mm)$r.squared)
                            }))

residuals <- invisible(apply(pw.m[,,,], 4,
                             function(im) {
                                 lm.res(im, terms = "g ~ upper + (b + w)", crop.data = T)
                             }))
                       
write.csv(data.frame(t(l.models)), paste0(fpath, "cropped-no-int.csv"), quote = F)
saveRDS(residuals, paste0(fpath, "res-cropped-no-int.rds"))
}

# no interaction, all data
{
    l.models <- invisible(apply(pw.m[,,,], 4, 
                                function(im) {
                                    mm <- quick.lm(im, terms = "g ~ upper + (b + w)", crop.data = F)
                                    c(c(coef(mm), "b:w" = 0)[1:5], rmse = summary(mm)$sigma, r2 = summary(mm)$r.squared)
                                }))
    
    residuals <- invisible(apply(pw.m[,,,], 4,
                                 function(im) {
                                     lm.res(im, terms = "g ~ upper + (b + w)", crop.data = F)
                                 }))
    
    write.csv(data.frame(t(l.models)), paste0(fpath, "uncropped-no-int.csv"), quote = F)
    saveRDS(residuals, paste0(fpath, "res-uncropped-no-int.rds"))
}

# with interaction, data cropped
{
    l.models <- invisible(apply(pw.m[,,,], 4, 
                                function(im) {
                                    mm <- quick.lm(im, terms = "g ~ upper + (b * w)", crop.data = T)
                                    c(c(coef(mm), "b:w" = 0)[1:5], rmse = summary(mm)$sigma, r2 = summary(mm)$r.squared)
                                }))
    
    residuals <- invisible(apply(pw.m[,,,19:22], 4,
                                 function(im) {
                                     lm.res(im, terms = "g ~ upper + (b * w)", crop.data = T)
                                 }))
    
    write.csv(data.frame(t(l.models)), paste0(fpath, "cropped-with-int.csv"), quote = F)
    saveRDS(residuals, paste0(fpath, "res-cropped-with-int.rds"))
}

# with interaction, all data
{
    l.models <- invisible(apply(pw.m[,,,], 4, 
                                function(im) {
                                    mm <- quick.lm(im, terms = "g ~ upper + (b * w)", crop.data = F)
                                    c(c(coef(mm), "b:w" = 0)[1:5], rmse = summary(mm)$sigma, r2 = summary(mm)$r.squared)
                                }))
    
    residuals <- invisible(apply(pw.m[,,,19:22], 4,
                                 function(im) {
                                     lm.res(im, terms = "g ~ upper + (b * w)", crop.data = F)
                                 }))
    
    write.csv(data.frame(t(l.models)), paste0(fpath, "uncropped-with-int.csv"), quote = F)
    saveRDS(residuals, paste0(fpath, "res-uncropped-with-int-1-5.rds"))
    saveRDS(residuals, paste0(fpath, "res-uncropped-with-int-6-11.rds"))
    saveRDS(residuals, paste0(fpath, "res-uncropped-with-int-12-18.rds"))
    saveRDS(residuals, paste0(fpath, "res-uncropped-with-int-19-22.rds"))
    
    r1 <- readRDS(paste0(fpath, "res-uncropped-with-int-1-5.rds"))
    r2 <- readRDS(paste0(fpath, "res-uncropped-with-int-6-11.rds"))
    r3 <- readRDS(paste0(fpath, "res-uncropped-with-int-12-18.rds"))
    r4 <- readRDS(paste0(fpath, "res-uncropped-with-int-19-22.rds"))
    
    rr <- array(abind(r1, r2, r3, r4, along = 2),
                dim = dim(pw.m[,,"black",]), dimnames = dimnames(pw.m[,,"black",]))
    
    saveRDS(rr, paste0(fpath, "res-uncropped-with-int.rds"))
}

####################################################################################################

# COMPARE MODELS                                                                                ####

models <- rbind.fill(list(data.frame(int = F, cropped = T,
                                read.csv(paste0(fpath, "cropped-no-int.csv"), row.names = NULL), stringsAsFactors = F),
               data.frame(int = F, cropped = F,
                                read.csv(paste0(fpath, "uncropped-no-int.csv"), row.names = NULL), stringsAsFactors = F),
               data.frame(int = T, cropped = T,
                                read.csv(paste0(fpath, "cropped-with-int.csv"), row.names = NULL), stringsAsFactors = F),
               data.frame(int = T, cropped = F,
                                read.csv(paste0(fpath, "uncropped-with-int.csv"), row.names = NULL), stringsAsFactors = F)))

# interaction makes very little difference between models. Cropping more so.
pdf(paste0(fpath, "linear-model-plots.pdf"), width = 4 * 4, height = 4 * 2); {
    par(mfrow = c(2,4))

    # RMSE by model
    {                       
        plot(models$rmse[models$int == T & models$cropped == T], 
             models$rmse[models$int == T & models$cropped == F], pch = 15,
             main = "RMSE - with interaction", xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
        
        plot(models$rmse[models$int == F & models$cropped == T], 
             models$rmse[models$int == F & models$cropped == F], pch = 15,
             main = "RMSE - without interaction", xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
        
        plot(models$rmse[models$int == T & models$cropped == T], 
             models$rmse[models$int == F & models$cropped == T], pch = 15,
             main = "RMSE - cropped data", xlab = "With interaction", ylab = "Without interaction")
        abline(0,1, col = "red")
        
        plot(models$rmse[models$int == T & models$cropped == F], 
             models$rmse[models$int == F & models$cropped == F], pch = 15,
             main = "RMSE - all data", xlab = "With interaction", ylab = "Without interaction")
        abline(0,1, col = "red")
    }
    
    # R2 by model
    {
        plot(models$r2[models$int == T & models$cropped == T], 
             models$r2[models$int == T & models$cropped == F], pch = 15,
             main = expression(paste(R^2, " - with interaction", collapse = "")), xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
        
        plot(models$r2[models$int == F & models$cropped == T], 
             models$r2[models$int == F & models$cropped == F], pch = 15,
             main = expression(paste(R^2, " - without interaction", collapse = "")), xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
        
        plot(models$r2[models$int == T & models$cropped == T], 
             models$r2[models$int == F & models$cropped == T], pch = 15,
             main = expression(paste(R^2, " - cropped data", collapse = "")), xlab = "With interaction", ylab = "Without interaction")
        abline(0,1, col = "red")
        
        plot(models$r2[models$int == T & models$cropped == F], 
             models$r2[models$int == F & models$cropped == F], pch = 15,
             main = expression(paste(R^2, " - all data", collapse = "")), xlab = "With interaction", ylab = "Without interaction")
        abline(0,1, col = "red")
    }
    dev.off()
}


####################################################################################################

# ROBUST MODEL FITTING                                                                          ####

quick.rlm <- function(im, terms, midline = 1024.5, crop.data = T) {
    
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    
    if(!is.na(midline)) {
        df$upper <- df$y > midline
    }
    
    # fit linear model to central part of image only (excludes edge effects)
    if (crop.data) {
        df <- df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                     findInterval(df$y, c(40.5, 2008.5)) == 1, ]
    }
    
    rlm(as.formula(terms), data = df)
}

r.models <- invisible(apply(pw.m[,,,], 4, 
                            function(im) {
                                mm <- quick.rlm(im, terms = "g ~ upper + (b + w)", crop.data = F)
                                c(c(coef(mm), "b:w" = 0)[1:5], rmse = summary(mm)$sigma, r2 = summary(mm)$r.squared)
                            }))

write.csv(data.frame(t(r.models)), paste0(fpath, "robust-uncropped-no-int.csv"), quote = F)

robust.models <- rbind.fill(list(data.frame(int = F, cropped = F,
                                     read.csv(paste0(fpath, "robust-uncropped-no-int.csv"), row.names = NULL), stringsAsFactors = F),
                          data.frame(int = F, cropped = T,
                                     read.csv(paste0(fpath, "robust-cropped-no-int.csv"), row.names = NULL), stringsAsFactors = F)))

# interaction makes very little difference between models. Cropping more so.
pdf(paste0(fpath, "robust-linear-model-plots.pdf"), width = 4 * 2, height = 4 * 2); {
    par(mfrow = c(2,2))
    
    # RMSE by model
    {                       
        plot(robust.models$rmse[!robust.models$int & robust.models$cropped], 
             robust.models$rmse[!robust.models$int & !robust.models$cropped], pch = 15,
             main = "RMSE - without interaction", xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
    }
    dev.off()
}

all.models <- rbind(data.frame(robust = F, models),
                    data.frame(robust = T, robust.models))
attach(all.models)

pdf(paste0(fpath, "robust-vs-cropped-linear-model-plots.pdf"), width = 4 * 2, height = 4 * 2); {
    par(mfrow = c(2,2))
    
    # RMSE by model
    {                       
        plot(rmse[!robust & !int & cropped], rmse[!robust & !int & !cropped], pch = 15,
             main = "RMSE without interaction - not robust", xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
        
        plot(rmse[robust & !int & cropped], rmse[robust & !int & !cropped], pch = 15,
             main = "RMSE without interaction - robust", xlab = "Cropped", ylab = "Uncropped")
        abline(0,1, col = "red")
  
        plot(rmse[!robust & !int & cropped], rmse[robust & !int & cropped], pch = 15,
             main = "RMSE without interaction - cropped", xlab = "Not robust", ylab = "Robust")
        abline(0,1, col = "red")
        
        plot(rmse[!robust & !int & !cropped], rmse[robust & !int & !cropped], pch = 15,
             main = "RMSE without interaction - uncropped", xlab = "Not robust", ylab = "Robust")
        abline(0,1, col = "red")
        
        plot(rmse[!robust & !int & cropped], rmse[robust & !int & cropped], pch = 15,
             main = "RMSE without interaction", xlab = "Not robust, cropped", ylab = "Robust, uncropped")
        abline(0,1, col = "red")
        }
    dev.off()
}

# robust model almost always has lower RMSE

####################################################################################################

# PROBLEMS WITH NON-ROBUST MODEL                                                                ####

# check homoscedasticity of original least-squares model
{
    # lowest observed RMSE in 130613
    mm <- quick.lm(pw.m[,,,"130613"], terms = "g ~ upper + (b + w)", crop.data = F)
    .smoothScatter(mm$fitted.values, mm$residuals)
    
    # robust RMSE higher than LS rmse in 131122
    
    # strong suggestion of heteroscedasticity
    # also large number of outliers
    # Robust method is definitely appropriate
}

# M-estimator is known not to be robust to outliers in the explanatory variables
# filter these out (using extreme/local thresholding) before fitting linear model.

####################################################################################################

# PRE-FILTER DATA                                                                               ####
dt <- "loan"

# filter out all locally/globally extreme pixels
df <- setNames(data.frame(melt(pw.m[, , "black", dt]), 
                          melt(pw.m[, , "grey", dt]), 
                          melt(pw.m[, , "white", dt]),
                          melt(md7[, , "black", dt]),
                          melt(md7[, , "grey", dt]))[, c("X1", "X2", "value", "value.1", "value.2", "value.3", "value.4")],
               nm = c("x", "y", "b", "g", "w", "b.res", "g.res"))
df$upper <- df$y > 1024.5
df <- df[!is.na(df$b),]

zz <- which(findInterval(df$b, asymmetric.mad(df$b)) == 1 & 
                      findInterval(df$g, asymmetric.mad(df$g)) == 1 & 
                      findInterval(df$b.res, asymmetric.mad(df$b.res)) == 1 & 
                      findInterval(df$g.res, asymmetric.mad(df$g.res)) == 1)

th.lm <- lm(g ~ upper + b + w, data = df[zz,])
th.rlm <- rlm(g ~ upper + b + w, data = df[zz,])

full.lm <- lm(g ~ upper + b + w, data = df)
full.rlm <- rlm(g ~ upper + b + w, data = df)

models <- cbind(do.call("rbind", lapply(list(full.lm, full.rlm, th.lm, th.rlm), coef)),
      do.call("rbind", lapply(lapply(list(full.lm, full.rlm, th.lm, th.rlm), summary), "[", "sigma")))
rownames(models) <- c("full.lm", "full.rlm", "th.lm", "th.rlm")

pdf(paste0(fpath, "Filtering-effect.pdf"), width = 5 * 4, height = 5 * 2); {
    par(mfrow = c(2,4))
    
    .smoothScatter(df$g, full.lm$fitted.values, main = "All data, least squares", xlab = "Observed", ylab = "Fitted",
                   xlim = c(0,65535), ylim = c(0,80000))
    abline(0,1, col = "darkred", lty = 2)
    text(1000,65000, adj = 0, paste0("Sigma ", round(summary(full.lm)$sigma,0)))
    text(1000,55000, adj = 0, paste0("RMSE ", round(sqrt(mean((df$g - predict(full.lm, df))^2)), 0)))
    
    .smoothScatter(df$g, full.rlm$fitted.values, main = "All data, M-estimate", xlab = "Observed", ylab = "Fitted",
                   xlim = c(0,65535), ylim = c(0,80000))
    abline(0,1, col = "darkred", lty = 2)
    text(1000,55000, adj = 0, paste0("Sigma ", round(summary(full.rlm)$sigma,0)))
    text(1000,65000, adj = 0, paste0("RMSE ", round(sqrt(mean((df$g - predict(full.rlm, df))^2)), 0)))
    
    .smoothScatter(df$g, predict(th.lm, df), main = "Filtered data, least squares", xlab = "Observed", ylab = "Fitted",
                   xlim = c(0,65535), ylim = c(0,80000))
    abline(0,1, col = "darkred", lty = 2)
    text(1000,65000, adj = 0, paste0("RMSE ", round(sqrt(mean((df$g - predict(th.lm, df))^2)), 0)))
    
    .smoothScatter(df$g, predict(th.rlm, df), main = "Filtered data, M-estimate", xlab = "Observed", ylab = "Fitted",
                   xlim = c(0,65535), ylim = c(0,80000))
    abline(0,1, col = "darkred", lty = 2)
    text(1000,65000, adj = 0, paste0("RMSE ", round(sqrt(mean((df$g - predict(th.rlm, df))^2)), 0)))
    
    .smoothScatter(df$g, full.lm$residuals, xlab = "Observed value", ylab = "Residual", ylim = c(-6000, 6000))
    abline(0, 0, col = "darkred", lty = 2)
    abline(h = c(2,-2) * summary(full.lm)$sigma, col = "cyan3")
    abline(h = asymmetric.mad(full.lm$residuals), col = "green3", lty = 2)
    text(1000,6000, adj = 0, paste0(sum(abs(full.lm$residuals) > 2 * summary(full.lm)$sigma), "px > 2sd"))
    text(1000,5000, adj = 0, paste0(sum(findInterval(full.lm$residuals, asymmetric.mad(full.lm$residuals)) %in% c(0,2)), " by MAD"))
    
    .smoothScatter(df$g, full.rlm$residuals, xlab = "Observed value", ylab = "Residual", ylim = c(-6000, 6000))
    abline(0, 0, col = "darkred", lty = 2)
    abline(h = c(2,-2) * summary(full.rlm)$sigma, col = "cyan3")
    abline(h = asymmetric.mad(full.rlm$residuals), col = "green3", lty = 2)
    text(1000,6000, adj = 0, paste0(sum(abs(full.rlm$residuals) > 2 * summary(full.rlm)$sigma), "px > 2sd"))
    text(1000,5000, adj = 0, paste0(sum(findInterval(full.rlm$residuals, asymmetric.mad(full.rlm$residuals)) %in% c(0,2)), " by MAD"))
    
    .smoothScatter(df$g, df$g - predict(th.lm, df), xlab = "Observed value", ylab = "Residual", ylim = c(-6000, 6000))
    abline(0, 0, col = "darkred", lty = 2)
    abline(h = c(2,-2) * summary(th.lm)$sigma, col = "cyan3")
    abline(h = asymmetric.mad(th.lm$residuals), col = "green3", lty = 2)
    text(1000,6000, adj = 0, paste0(sum(abs(th.lm$residuals) > 2 * summary(th.lm)$sigma), "px > 2sd"))
    text(1000,5000, adj = 0, paste0(sum(findInterval(th.lm$residuals, asymmetric.mad(th.lm$residuals)) %in% c(0,2)), " by MAD"))
    
    .smoothScatter(df$g, df$g - predict(th.rlm, df), xlab = "Observed value", ylab = "Residual", ylim = c(-6000, 6000))
    abline(0, 0, col = "darkred", lty = 2)
    abline(h = c(2,-2) * summary(th.rlm)$sigma, col = "cyan3")
    abline(h = asymmetric.mad(th.rlm$residuals), col = "green3", lty = 2)
    text(1000,6000, adj = 0, paste0(sum(abs(th.rlm$residuals) > 2 * summary(th.rlm)$sigma), "px > 2sd"))
    text(1000,5000, adj = 0, paste0(sum(findInterval(th.rlm$residuals, asymmetric.mad(th.rlm$residuals)) %in% c(0,2)), " by MAD"))
    
    dev.off()
}


