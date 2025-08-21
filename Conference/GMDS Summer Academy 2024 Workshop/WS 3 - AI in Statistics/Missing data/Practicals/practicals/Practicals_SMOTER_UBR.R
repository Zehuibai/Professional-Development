################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Implementation of SMOTER and UBR for MD           ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################


library(UBL)
library(randomForest)
library(plyr)
library(lsmeans)



############################################################################################
##################################### Test group ###########################################
############################################################################################

rm(list = ls())

all=read.csv("C:\\Users\\206586\\OneDrive - BioNTech SE\\HH\\Development\\GMDS 2024\\Simulations\\sep24\\simu_score_all.csv")

tapply(all$true, all$treat, summary)
tapply(all$score, all$treat, summary)


grp <- 1 

train.org <- subset(all, treat==grp & misstypen==9)
miss      <- subset(all, treat==grp & (misstypen==1 | misstypen==2))

tgt <- which(colnames(all) == "score")


boxplot(train.org$score)
summary(train.org$score)

# using a relevance function provided by the user 
rel <- matrix(0, ncol = 3, nrow = 0)
rel <- rbind(rel, c(0,  1,   0))
rel <- rbind(rel, c(3,  1,   0))
rel <- rbind(rel, c(5,  0.9, 0))
rel <- rbind(rel, c(7,  0.8, 0))
rel <- rbind(rel, c(9,  0.7, 0))
rel <- rbind(rel, c(10, 0.6, 0))
rel <- rbind(rel, c(12, 0.5, 0))
rel <- rbind(rel, c(13, 0,   0))
rel <- rbind(rel, c(15, 0,   0))
rel <- rbind(rel, c(17, 0,   0))
rel <- rbind(rel, c(20, 0,   0))

plot(rel[,1],rel[,2], ylim=c(-0.2,1.2),type = "l", xlab="Target variable", ylab="Relevance value")
abline(h=c(0.50),lty=3,col="red")

train <- SmoteRegress(score ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7,
                      rel = rel, thr.rel = 0.55,  k = 5,
                      train.org, dist = "Euclidean", C.perc = list(2,1))
boxplot(train.org$score)
boxplot(train$score)

hist(train.org$score)
hist(train$score)


plot(train.org$cov2,train.org$score, xlab="Covariate", ylab="Target variable")
plot(train$cov2,train$score, xlab="Covariate", ylab="Target variable")
summary(train.org$score)
summary(train$score)


# Estimation of parameters used for obtaining the relevance function.
control.parms <- phi.control(train[,tgt], method="extremes", extr.type = "both", coef=0.2)
control.parms$control.pts
#  Options for methods: "extremes" or "range".
# If the selected method is "extremes", the distribution of the target variable values is used to assign more importance to the most extreme values according to the boxplot. 
# If the selected method is "range", a matrix should be provided defining the important and unimportant values


### Relevance function.
# This function allows to obtain the relevance function values on a set of target variable values given the interpolating points.
y.phi <- phi(sort(train[,tgt]), control.parms=control.parms)
plot(sort(train[,tgt]), y.phi, ylim=c(0,1), type = "l", xlab = "Target variable", ylab = "Relevance value")



# the boundaries of the domain considered
minds <- min(train[,tgt])
maxds <- max(train[,tgt])


### A 3-column matrix with interpolating points for the off-diagonal cases (i.e., y != y.pred), provided by the user. 
# The first column has the y value, the second column the y.pred value and the third column has the corresponding utility value.
# At least, the off diagonal domain boundary points (i.e., points (minds, maxds, util) and (maxds, minds, util) ) must be provided in this matrix. 
# Moreover, the points provided through this parameter must be in [minds, maxds] range.
# build m.pts to include at least (minds, maxds) and (maxds, minds) points

m.pts <- matrix(c(minds, maxds, -1, maxds, minds, -1), byrow=TRUE, ncol=3)

miss$score <- 15 # Due to the bug in UtilOptimRegress function (doesn't support missing value), assign a value to missing values. The value of chosen here (e.g., value=1) doesn't impact the prediction results. 

pred.res <- UtilOptimRegress(score~cov1+cov2+cov3+cov4+cov5+cov6+cov7, train=train, test=miss,
                             type = "utility", strat = "interpol",
                             strat.parms=list(method = "bilinear"),
                             control.parms = control.parms,
                             m.pts = m.pts, minds = minds, maxds = maxds)


miss.pred <- miss
miss.pred$ubl.pred <- pred.res$optim

reg.true.pred <- lm(pred.res$optim ~ miss$true) # this is function from stats package, not from caret
plot(miss$true, miss.pred$ubl.pred, xlim=c(0,20), ylim=c(0,20), xlab="True value", ylab="UBL Prediction")
abline(coef = c(0,1), col="blue")
abline(reg.true.pred, col="red")
summary(reg.true.pred)


rf <- randomForest(score ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7, data=train)
miss.pred$rf.pred <- predict(rf, miss.pred)
reg.true.pred.rf <- lm(miss.pred$rf.pred ~ miss$true)
plot(miss.pred$true, miss.pred$rf.pred, xlim=c(0,20), ylim=c(0,20), xlab="True value", ylab="RF SMOTER Prediction")
abline(coef = c(0,1), col="blue")
abline(reg.true.pred.rf, col="red")
summary(reg.true.pred.rf)


rf.rog <- randomForest(score ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7, data=train.org)
miss.pred$rf.org.pred <- predict(rf.rog, miss.pred)
reg.true.pred.rf.org <- lm(miss.pred$rf.org.pred ~ miss$true)
plot(miss.pred$true, miss.pred$rf.org.pred, xlim=c(0,20), ylim=c(0,20), xlab="True value", ylab="RF Prediction")
abline(coef = c(0,1), col="blue")
abline(reg.true.pred.rf.org, col="red")
summary(reg.true.pred.rf.org)


train.org$ubl.pred <- NA 
train.org$rf.pred <- NA
train.org$rf.org.pred <- NA

impt.test <- rbind(train.org, miss.pred)


impt.test$impt.ubl <- with(impt.test, ifelse(misstype=="MCAR"|misstype=="MAR"|misstype=="MNAR", impt.test$ubl.pred, impt.test$score))
impt.test$impt.rf <- with(impt.test, ifelse(misstype=="MCAR"|misstype=="MAR"|misstype=="MNAR", impt.test$rf.pred, impt.test$score))
impt.test$impt.rf.org <- with(impt.test, ifelse(misstype=="MCAR"|misstype=="MAR"|misstype=="MNAR", impt.test$rf.org.pred, impt.test$score))


summary(impt.test$true)
t.test(impt.test$true, mu = 12, alternative = "two.sided")

summary(train.org$score)
t.test(train.org$score, mu = 12, alternative = "two.sided")

summary(impt.test$impt.ubl)
t.test(impt.test$impt.ubl, mu = 12, alternative = "two.sided")

summary(impt.test$impt.rf)
t.test(impt.test$impt.rf, mu = 12, alternative = "two.sided")

summary(impt.test$impt.rf.org)
t.test(impt.test$impt.rf.org, mu = 12, alternative = "two.sided")


############################## End of Test Group ##############################




############################################################################################
#####################################  Control Group  ######################################
############################################################################################




############################## End of Control Group ##############################




############################# Analysis of data with both groups combined ##############################



############################## End of Analysis ##############################

