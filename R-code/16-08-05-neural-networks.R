

library("IO.Pixels"); library("CB.Misc")
pw.m <- load.pixel.means()

norm <- pw.m[,,"black","160430"]

norm <- 1 - norm / mean(norm, na.rm = T)
pixel.image(norm, xlim = c(900,1100), ylim = c(0,256))

hist(norm, breaks = "fd")

pixel.plot(which(norm < -1, arr.ind = T), cex = 0.4)

tmp <- norm
tmp[tmp < -1] <- NA
pixel.image(tmp, xlim = c(900,1100), ylim = c(0,256), break.levels = sd.levels(norm))

####################################################################################################

# NEURALNET EXAMPLE                                                                             ####

# https://www.r-bloggers.com/fitting-a-neural-network-in-r-neuralnet-package/

set.seed(500)
library(MASS)
data <- Boston

# check for missing data
apply(data,2,function(x) sum(is.na(x)))

# randomly split data into training & test sets
index <- sample(1:nrow(data), round(0.75 * nrow(data)))
train <- data[index,]
test <- data[-index,]

# fit linear model on training set, test on test set
lm.fit <- glm(medv ~., data = train)
summary(lm.fit)
pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$medv)^2)/nrow(test)

                        #---------------------------------------------------#
                        ###          PREP & TRAIN NEURAL NETWORK          ###
                        #---------------------------------------------------#

# normalise data by min-max scaling
maxs <- apply(data, 2, max) 
mins <- apply(data, 2, min)

scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))

train_ <- scaled[index,]
test_ <- scaled[-index,]

library(neuralnet)
n <- names(train_)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)

# NOTES:

# For some reason the formula y~. is not accepted in the neuralnet() function.
# You need to first write the formula and then pass it as an argument in the fitting function.

# The hidden argument accepts a vector with the number of neurons for each hidden layer
# The argument linear.output is used to specify whether we want to do 
# regression (linear.output=TRUE) or classification (linear.output=FALSE)

plot(nn)

# use neural network to predict values
pr.nn <- compute(nn,test_[,1:13])

# un-normalise data
pr.nn_ <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
test.r <- (test_$medv)*(max(data$medv)-min(data$medv))+min(data$medv)

MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test_)

# compare the two MSEs

round(c(MSE.lm, MSE.nn), 1)


                #---------------------------------------------------#
                ###              EVALUATE RESULTS                 ###
                #---------------------------------------------------#


par(mfrow = c(1,2))

plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(test$medv,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)

par(mfrow = c(1,1))


plot(test$medv,pr.nn_,col='red',main='Real vs predicted',pch=18,cex=0.7)
points(test$medv,pr.lm,col='blue',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend=c('NN','LM'),pch=18,col=c('red','blue'))


# fast cross-validation

library(boot)
set.seed(200)
lm.fit <- glm(medv~.,data=data)
cv.glm(data,lm.fit,K=10)$delta[1]

# cross-validation of neural network
set.seed(450)
cv.error <- NULL
k <- 10

library(plyr) 
pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
    index <- sample(1:nrow(data),round(0.9*nrow(data)))
    train.cv <- scaled[index,]
    test.cv <- scaled[-index,]
    
    nn <- neuralnet(f,data=train.cv,hidden=c(5,2),linear.output=T)
    
    pr.nn <- compute(nn,test.cv[,1:13])
    pr.nn <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
    
    test.cv.r <- (test.cv$medv)*(max(data$medv)-min(data$medv))+min(data$medv)
    
    cv.error[i] <- sum((test.cv.r - pr.nn)^2)/nrow(test.cv)
    
    pbar$step()
}


boxplot(cv.error,xlab='MSE CV',col='cyan3',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)

