
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

attach(all.models)
plot(all.models$b[!int & !cropped & !robust], all.models$b[!int & !cropped & robust], pch = 20, 
     xlab = "robust", ylab = "least squares", main = "Coefficient b")
abline(0,1,col = "darkred", lty = 2)

plot(w[!int & !cropped & !robust], w[!int & !cropped & robust], pch = 20, 
     xlab = "robust", ylab = "least squares", main = "Coefficient w")
abline(0,1,col = "darkred", lty = 2)

plot(b[!int & !cropped & !robust], X.Intercept.[!int & !cropped & !robust], pch = 20, 
     xlab = "robust", ylab = "least squares", main = "Black coef vs intercept")
abline(1300,-3000,col = "darkred", lty = 2)

line(b[!int & !cropped & !robust], X.Intercept.[!int & !cropped & !robust])

# white coefficient fitted by both approaches is v similar.

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

####################################################################################################

# SIMULATE EFFECT OF INCREASING NUMBERS OF EXTREME PIXELS                                       ####

# use 141009 as base data
im <- pw.m[,,,"141009"]

df <- setNames(data.frame(melt(im[, , "black"]), 
                          melt(im[, , "grey"]), 
                          melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
               nm = c("x", "y", "b", "g", "w"))
df$upper <- df$y > 1024.5
df <- df[!is.na(df$g),]

# standard model
{
    rlm.141009 <- rlm(g ~ upper + b + w, data = df)
    lm.141009 <- lm(g ~ upper + b + w, data = df)
}
# fit model with 50000 extreme white values
{
    df$w50000 <- df$w; df$w50000[sample(1:nrow(df), 50000)] <- 65535
    
    rlm.w50000 <- rlm(g ~ upper + b + w50000, data = df)
    r.res.w50000 <- rlm.w50000$residuals
    .smoothScatter(df$g, rlm.w50000$fitted.values); abline(0,1, col = "darkred", lty=2)
    
    coef(rlm.w50000); summary(rlm.w50000)$sigma
        #  (Intercept)    upperTRUE            b       w50000          sigma
        # -350.5630576   11.7636932    0.6574950    0.3404426       59.28925
    
    # least squares fitting
    lm.w50000 <- lm(g ~ upper + b + w50000, data = df)
    ls.res.w50000 <- lm.w50000$residuals
    .smoothScatter(df$g, lm.w50000$fitted.values); abline(0,1, col = "darkred", lty=2)
    
    coef(lm.w50000); summary(lm.w50000)$sigma
        #  (Intercept)    upperTRUE            b       w50000          sigma
        # 8477.6804073   85.9265596    0.5364828    0.1716512       462.0374
}

# fit model with 50000 extreme black values - effect is worse
{
    df$b50000 <- df$b; df$b50000[sample(1:nrow(df), 50000)] <- 65535

    rlm.b50000 <- rlm(g ~ upper + b50000 + w, data = df)
    r.res.b50000 <- rlm.b50000$residuals
    .smoothScatter(df$g, rlm.b50000$fitted.values, main = "Robust model, 50000 extreme black values")
    abline(0,1, col = "darkred", lty=2)
    
    coef(rlm.b50000); summary(rlm.b50000)$sigma
    #  (Intercept)    upperTRUE            w       b50000          sigma
    # 4937.9656502   27.6781716    0.3038765    0.0005827       147.7138
    
    # least squares fitting
    lm.b50000 <- lm(g ~ upper + b50000 + w, data = df)
    ls.res.b50000 <- lm.b50000$residuals
    .smoothScatter(df$g, lm.b50000$fitted.values, main = "Least squares model, 50000 extreme black values")
    abline(0,1, col = "darkred", lty=2)
    
    coef(lm.b50000); summary(lm.b50000)$sigma
    #  (Intercept)    upperTRUE            b       b50000          sigma
    # 17610.660067   159.714954     0.397766    0.0000005   642.7129
}

# fit model with 50000 increased black values (+12000)
{
    df$bx50000 <- df$b; 
    s <- sample(1:nrow(df), 50000); df$bx50000[s] <- df$b[s] + 12000
    
    rlm.bx50000 <- rlm(g ~ upper + bx50000 + w, data = df)
    r.res.bx50000 <- rlm.bx50000$residuals
    .smoothScatter(df$g, rlm.bx50000$fitted.values, main = "Robust model, 50000 black values + 12000")
    abline(0,1, col = "darkred", lty=2)
    
    round(coef(rlm.bx50000), 7); summary(rlm.bx50000)$sigma
    #  (Intercept)    upperTRUE            w      bx50000         sigma
    # 4792.2745541   27.2047898    0.3049917    0.0173753      145.7641
    
    # least squares fitting
    lm.bx50000 <- lm(g ~ upper + bx50000 + w, data = df)
    ls.res.bx50000 <- lm.bx50000$residuals
    .smoothScatter(df$g, lm.bx50000$fitted.values, main = "Least squares model, 50000 black values + 12000")
    abline(0,1, col = "darkred", lty=2)
    
    round(coef(lm.bx50000), 7); summary(lm.bx50000)$sigma
    #  (Intercept)    upperTRUE            w      bx50000       sigma
    # 877.9295939   20.1988309    0.3179752    0.0694438     296.8281
}

# fit model with 5000 extreme white values
{
    df$w5000 <- df$w; df$w5000[sample(1:nrow(df), 5000)] <- 65535
    
    rlm.w5000 <- rlm(g ~ upper + b + w5000, data = df)
    r.res.w5000 <- rlm.w5000$residuals
    .smoothScatter(df$g, rlm.w5000$fitted.values, main = "Robust model, 5000 extreme white values")
    abline(0,1, col = "darkred", lty=2)
    
    coef(rlm.w5000); summary(rlm.w5000)$sigma
    #  (Intercept)    upperTRUE            b        w5000          sigma
    # -797.1793570    8.6306691    0.6788071    0.3473058       57.25373
    
    # least squares fitting
    lm.w5000 <- lm(g ~ upper + b + w5000, data = df)
    ls.res.w5000 <- lm.w5000$residuals
    .smoothScatter(df$g, lm.w5000$fitted.values, main = "Least squares model, 5000 extreme white values")
    abline(0,1, col = "darkred", lty=2)
    
    coef(lm.w5000); summary(lm.w5000)$sigma
    #  (Intercept)    upperTRUE            b       w5000          sigma
    #  871.2298981   22.8992816    0.6554789   0.3154575       205.4668
}

# fit model with 5000 extreme black values
{
    df$b5000 <- df$w; df$b5000[sample(1:nrow(df), 5000)] <- 65535
    
    rlm.b5000 <- rlm(g ~ upper + b5000 + w, data = df)
    r.res.b5000 <- rlm.b5000$residuals
    .smoothScatter(df$g, rlm.b5000$fitted.values, main = "Robust model, 5000 extreme black values")
    abline(0,1, col = "darkred", lty=2)
    
    round(coef(rlm.b5000), 7); summary(rlm.b5000)$sigma
    #  (Intercept)    upperTRUE            w        b5000          sigma
    # 4943.3638726   27.6964262    0.3038848   -0.0000463       147.7879
    
    # least squares fitting
    lm.b5000 <- lm(g ~ upper + b5000 + w, data = df)
    ls.res.b5000 <- lm.b5000$residuals
    .smoothScatter(df$g, lm.b5000$fitted.values, main = "Least squares model, 5000 extreme black values")
    abline(0,1, col = "darkred", lty=2)
    
    round(coef(lm.b5000), 7); summary(lm.b5000)$sigma
    #  (Intercept)    upperTRUE            b        b5000          sigma
    # 4426.3433383   21.5855028    0.3146186   -0.0000687       312.5519
}

# fit model with thresholded & trimmed data
{
    df$th <- findInterval(df$b, asymmetric.mad(df$b)) * findInterval(df$g, asymmetric.mad(df$g))
    lm.th <- lm(g ~ upper + b + w, data = df[df$th == 1,])
    rlm.th <- rlm(g ~ upper + b + w, data = df[df$th == 1,])
}

zz <- cbind(rbind(coef(lm.141009), coef(rlm.141009),
                  coef(lm.w5000), coef(rlm.w5000), coef(lm.b5000), coef(rlm.b5000),
                  coef(lm.w50000), coef(rlm.w50000), coef(lm.b50000), coef(rlm.b50000),
                  coef(lm.bx50000), coef(rlm.bx50000)),
            c(summary(lm.141009)$sigma, summary(rlm.141009)$sigma,
              summary(lm.w5000)$sigma, summary(rlm.w5000)$sigma, summary(lm.b5000)$sigma, summary(rlm.b5000)$sigma,
              summary(lm.w50000)$sigma, summary(rlm.w50000)$sigma, summary(lm.b50000)$sigma, summary(rlm.b50000)$sigma,
              summary(lm.bx50000)$sigma, summary(rlm.bx50000)$sigma))

pp <- cbind(rbind(coef(lm.th), coef(rlm.th)), 
            c(summary(lm.th)$sigma, summary(rlm.th)$sigma),
            c(F, T), c("-", "-"))
      
models <- rbind(models, setNames(data.frame(pp), nm = colnames(models)))
rownames(models)[13:14] <- c("lm.th", "rlm.th")


models <- setNames(data.frame(zz), nm = c("os", "u", "b", "w", "sigma"))
rownames(models) <- c("lm", "rlm", "lm.w5000", "rlm.w5000", "lm.b5000", "rlm.b5000",
                      "lm.w50000", "rlm.w50000", "lm.b50000", "rlm.b50000", "lm.bx50000", "rlm.bx50000")
models$robust <- rep(c(F, T), 6)
models$var <- c("-", "-", "w", "w", "b", "b", "w", "w", "b", "b", "b", "b", "th", "th")

write.csv(models, paste0(fpath, "models-with-extreme-values.csv"), quote = F)

attach(models)

format(models, scientific = F)

plot(os, b, pch = c(4, 20, 0, 1)[as.factor(var)], col = c("black", "green3")[(robust == T)+1])
# changes to black values have strong effect on coefficient of b - absorbed by change in offset.

plot(b, w, pch = c(4, 20, 0, 1)[as.factor(var)], col = c("black", "green3")[(robust == T)+1])
# coefficient of w is fairly consistent - dominates the model

plot(w, sigma, pch = c(4, 20, 0, 1)[as.factor(var)], col = c("black", "green3")[(robust == T)+1])
# robust modelling always has lower sigma
# threshold the data to be extre sure (since we're already doing so)

####################################################################################################

# TO FILTER OR NOT TO FILTER?                                                                   ####
fpath <- "./Image-plots/linear-res/"
acq.names <- c("130613", "130701", "131002", "131122", "140128", "140129", "141009", "141118", "141217",
               "150108", "150113", "150126", "150529", "150730", "150828", "151015", "160314", "160430",
               "160705", "loan2", "loan", "MCT225")

invisible(lapply(acq.names,
                 function(dt) {
                     im <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
                     df <- setNames(data.frame(melt(im[, , "black", dt]), 
                                               melt(im[, , "grey", dt]), 
                                               melt(im[, , "white", dt]))[, c("X1", "X2", "value", "value.1", "value.2")],
                                    nm = c("x", "y", "b", "g", "w"))
                     df$upper <- df$y > 1024.5
                     df <- df[!is.na(df$g),]
                     df$th <- findInterval(df$b, asymmetric.mad(df$b)) * findInterval(df$g, asymmetric.mad(df$g))
                     
                     lm.full <- lm(g ~ upper + b + w, data = df)
                     rlm.full <- rlm(g ~ upper + b + w, data = df)
                     lm.th <- lm(g ~ upper + b + w, data = df[df$th == 1,])
                     rlm.th <- rlm(g ~ upper + b + w, data = df[df$th == 1,])
                     
                     coeffs <- rbind(coef(lm.full), coef(rlm.full), coef(lm.th), coef(rlm.th))
                     sig <- c(summary(lm.full)$sigma, 
                              summary(rlm.full)$sigma, 
                              summary(lm.th)$sigma,
                              summary(rlm.th)$sigma)
                     px.removed <- rep(c(0, sum(df$th != 1)), each = 2)
                     npx <- c(sum(findInterval(lm.full$residuals, asymmetric.mad(lm.full$residuals)) != 1),
                              sum(findInterval(rlm.full$residuals, asymmetric.mad(rlm.full$residuals)) != 1),
                              sum(findInterval(lm.th$residuals, asymmetric.mad(lm.th$residuals)) != 1),
                              sum(findInterval(rlm.th$residuals, asymmetric.mad(rlm.th$residuals)) != 1))
                     
                     df$res.lm.full <- df$g - predict(lm.full, df)
                     df$res.rlm.full <- df$g - predict(rlm.full, df)
                     df$res.lm.th <- df$g - predict(lm.th, df)
                     df$res.rlm.th <- df$g - predict(rlm.th, df)
                     
                     zz <- setNames(data.frame(model = c("lm.full", "rlm.full", "lm.th", "rlm.th"), 
                                               robust = rep(c(F, T), 2),
                                               filtered = rep(c(F, T), each = 2),
                                               coeffs, sig, px.removed, npx),
                                    c("model", "robust", "filtered", "os", "u", "b", "w", "sigma", "omitted", "nonlinear"))
                     
                     write.csv(zz, paste0(fpath, "th-vs-full-", dt, ".csv"), quote = F, row.names = F)
                     
                     bmp(paste0(fpath, "th-vs-full-px-", dt, ".bmp"), height = 960 * 2, width = 960 * 2); {
                         par(mfrow = c(2,2), mar = c(2,2,3,1))
                         pixel.plot(df[df$th != 1,], col = "skyblue", main = paste0("Least squares, all data - ", dt))
                         points(df[findInterval(df$res.lm.full, asymmetric.mad(df$res.lm.full)) != 1,],
                                pch = 15, cex = 0.4)
                         
                         pixel.plot(df[df$th != 1,], col = "skyblue", main = paste0("Robust, all data - ", dt))
                         points(df[findInterval(df$res.rlm.full, asymmetric.mad(df$res.rlm.full)) != 1,],
                                pch = 15, cex = 0.4)
                         
                         pixel.plot(df[df$th != 1,], col = "skyblue", main = paste0("Least squares, filtered"))
                         points(df[findInterval(df$res.lm.th, asymmetric.mad(df$res.lm.th)) != 1 & df$th == 1,],
                                pch = 15, cex = 0.4)
                         
                         pixel.plot(df[df$th != 1,], col = "skyblue", main = paste0("Robust, filtered"))
                         points(df[findInterval(df$res.rlm.th, asymmetric.mad(df$res.rlm.th)) != 1 & df$th == 1, ],
                                pch = 15, cex = 0.4)
                         dev.off()
                     }
                     
}))

qq <- sapply(acq.names, 
             function(dt) {
                 mm <- read.csv(paste0(fpath,  "th-vs-full-", dt, ".csv"), row.names = 1)
}, simplify = F)

qq.all <- rbind.fill(qq)
rownames(qq.all) <- paste(rep(acq.names, each = 4), rep(row.names(qq[[1]]),22), sep = ".")

write.csv(qq.all, paste0(fpath, "all-models-th-robust.csv"), quote = F, row.names = T)


attach(qq.all)

# check # pixels identified with thresholding & without
# orange = LS fitted, black = robust
plot(nonlinear[filtered], nonlinear[!filtered], pch = 20, xlab = "Filtered", ylab = "Unfiltered", 
     main = "Nonlinear px identified", col = c("orange", "black")[robust+1])
abline(0,1,col = "red", lty = 2)
# unfair test: if px were high, values were removed. Large areas of v. high-valued px in corners

plot(sigma[filtered], sigma[!filtered], pch = 20, xlab = "Filtered", ylab = "Unfiltered", 
     main = "RMSE", col = c("orange", "black")[robust+1])
abline(0,1,col = "red", lty = 2)

# RMSE almost exactly the same in robust fitting. Higher without filtering in least-squares fitting.


# what about increasing large # pixels by an offset?
dt <- "141009"

im <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
df <- setNames(data.frame(melt(im[, , "black", dt]), 
                          melt(im[, , "grey", dt]), 
                          melt(im[, , "white", dt]))[, c("X1", "X2", "value", "value.1", "value.2")],
               nm = c("x", "y", "b", "g", "w"))
df$upper <- df$y > 1024.5

# additional offset, uniform between 10000 & 30000, truncated at max GV
df$b[s] <- df$b[s] + a; df$g[s] <- df$g[s] + a; df$w[s] <- df$w[s] + a
df[df > 65535] <- 65535

rlm.full <- rlm(g ~ upper + b + w, data = df)

coef(rlm.full)

.smoothScatter(df$g, df$w)

####################################################################################################

# FIT ALL LINEAR MODELS, STORE RESIDUALS                                                        ####

acq.names <- c("130613", "130701", "131002", "131122", "140128", "140129", "141009", "141118", "141217",
               "150108", "150113", "150126", "150529", "150730", "150828", "151015", "160314", "160430",
               "160705", "loan2", "loan", "MCT225")

linear.models <- list()

# get residuals for each image in turn
invisible(lapply(acq.names,
                 function(dt) {
                     im <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = dt)
                     df <- setNames(data.frame(melt(im[, , "black", dt]), 
                                               melt(im[, , "grey", dt]), 
                                               melt(im[, , "white", dt]))[, c("X1", "X2", "value", "value.1", "value.2")],
                                    nm = c("x", "y", "b", "g", "w"))
                     df$upper <- df$y > 1024.5

                     model <- rlm(g ~ upper + b + w, data = df)
                     res <- array(df$g - predict(model, newdata = df), dim = dim(im[,,1, dt]))
                     saveRDS(res, paste0("./02_Objects/linear-res/l.res-", dt, ".rds"))
                     
                     linear.models[[dt]] <<- c(coef(model), sigma = summary(model)$sigma, 
                                              rmse = sqrt(sum(res^2, na.rm = T) / sum(!is.na(res))))
                 }))

all.models <- do.call("rbind", linear.models)

write.csv(all.models, "./02_Objects/linear-res/models-fitted.csv", quote = F, row.names = T)

####################################################################################################

# LINEAR RESPONSE VS SHADING CORRECTION?                                                        ####

im <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = "loan")[,,,1]
l.res <- load.objects("./02_Objects/linear-res/", otype = "l-res", acq.list = "loan")[,,1]
sc <- shading.corrected(im)

sc <- sc[!is.na(l.res)]; l.res <- l.res[!is.na(l.res)]
.smoothScatter(sc, l.res, xlab = "Shading-corrected value", ylab = "Linear residuals")
abline(line(sc, l.res), col = "darkred", lty = 2)

####################################################################################################

# FITTED RESPONSE VS EXPECTED?                                                                  ####

l.res <- load.objects("./02_Objects/linear-res/", otype = "l-res", acq.list = c("141009", "loan", "160430"))
pw.m <- load.objects("./02_Objects/images/", otype = "pwm", acq.list = c("141009", "loan", "160430"))

l.model <- pw.m[,,"grey",] - l.res
exp <- array(apply(pw.m, 4, function(im) 0.75 * im[,,"black"] + 0.25 * im[,,"white"]), dim = dim(pw.m[,,"white",]), dimnames = dimnames(pw.m[,,"grey", ]))

# fitted vs expected
{
    .smoothScatter(l.model[,,"160430"], exp[,,"160430"], xlab = "Fitted", ylab = "Expected", main = "160430")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(l.model[,,"141009"], exp[,,"141009"], xlab = "Fitted", ylab = "Expected", main = "141009")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(l.model[,,"loan"], exp[,,"loan"], xlab = "Fitted", ylab = "Expected", main = "loan")
    abline(0,1,col = "darkred", lty = 2)
}

# expected vs actual
{
    .smoothScatter(pw.m[,,"grey","160430"], exp[,,"160430"], xlab = "Observed", ylab = "Expected", main = "160430")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(pw.m[,,"grey","141009"], exp[,,"141009"], xlab = "Observed", ylab = "Expected", main = "141009")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(pw.m[,,"grey","loan"], exp[,,"loan"], xlab = "Observed", ylab = "Expected", main = "loan")
    abline(0,1,col = "darkred", lty = 2)
}

# fitted vs actual
{
    .smoothScatter(pw.m[,,"grey","160430"], l.model[,,"160430"], xlab = "Observed", ylab = "Fitted", main = "160430")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(pw.m[,,"grey","141009"], l.model[,,"141009"], xlab = "Observed", ylab = "Fitted", main = "141009")
    abline(0,1,col = "darkred", lty = 2)
    
    .smoothScatter(pw.m[,,"grey","loan"], l.model[,,"loan"], xlab = "Observed", ylab = "Fitted", main = "loan")
    abline(0,1,col = "darkred", lty = 2)
}

####################################################################################################

# COMPARISON OF ASYMMETRIC MAD BY IMAGE                                                         ####

m.b <- apply(pw.m[,,"black",], 3, function(im) abs(modal.density(im) - asymmetric.mad(im, n = 1)))
m.g <- apply(pw.m[,,"grey",], 3, function(im) abs(modal.density(im) - asymmetric.mad(im, n = 1)))
m.w <- apply(pw.m[,,"white",], 3, function(im) abs(modal.density(im) - asymmetric.mad(im, n = 1)))

all.mad <- data.frame(b.l = m.b[1,], b.u = m.b[2,],
                      g.l = m.g[1,], g.u = m.g[2,],
                      w.l = m.w[1,], w.u = m.w[2,])

plot(m.b[1,], m.b[2,], pch = 20, xlim = c(0,10000), ylim = c(0,10000))
points(m.g[1,], m.g[2,], pch = 20, col = "green3")
points(m.w[1,], m.w[2,], pch = 20, col = "gold")


invisible(lapply(dimnames(pw.m)[[4]],
               function(dt) {
                   bmp(paste0("./Image-plots/Shading-corrections/shading-correction-", dt, ".bmp"),
                       width = 2048, height = 2048, pointsize = 28); {
                       pixel.image(shading.corrected(pw.m[,,,dt]), 
                                   title = paste0(dt, " - shading corrected"))
                           dev.off()
                   }
               }))