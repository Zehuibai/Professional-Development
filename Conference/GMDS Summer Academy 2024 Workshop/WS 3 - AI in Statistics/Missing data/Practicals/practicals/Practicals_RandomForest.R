################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Implementation of Random Forest for MD            ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################

#install.packages("randomForest")
library(randomForest)


############################################################################################
##################################### Test group ###########################################
############################################################################################

rm(list = ls())

all=read.csv("C:\\Users\\206586\\OneDrive - BioNTech SE\\HH\\Development\\GMDS 2024\\Simulations\\sep24\\simu_score_all.csv")

tapply(all$true, all$treat, summary)
tapply(all$score, all$treat, summary)


grp <- 1 # Test Group #
train <- subset(all, treat==grp & misstypen==9)
miss  <- subset(all, treat==grp & (misstypen==1 | misstypen==2))


summary(train$score)
t.test(train$score, mu = 12, alternative = "two.sided")
boxplot(train$score)


miss.pred <- miss #copy the data with missing value, prepare for prediction #


### Using default setting in random forest: 
  # ntree = 500
  # mtry = floor(ncol(x)/3) for regression 
  # replace = TRUE
  # sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)) 
  # nodesize = 5 for regression. Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). 
  # maxnodes = null. Maximum number of terminal nodes trees in the forest can have. If not given, trees are grown to the maximum possible (subject to limits by nodesize). If set larger than maximum possible, a warning is issued. 



rf <- randomForest(score ~ cov1+cov2+cov3+cov4+cov5+cov6+cov7, data=train)

print(rf)
rf$mse
#mean(rf$mse)
sqrt(mean(rf$mse))  # RMSE or MAE #
rf$rsq
rf$importance
plot(rf)

# use model rf to predict the missing data 
miss.pred$rf.pred <- predict(rf, miss.pred)

# add the prediction to the data frame
reg.true.pred.rf <- lm(miss.pred$rf.pred ~ miss$true)

# check predicted values vs true values (in simulation dataset)
reg.true.pred.rf <- lm(miss.pred$rf.pred ~ miss$true)
plot(miss.pred$true, miss.pred$rf.pred, xlim=c(0,20), ylim=c(0,20), xlab="True value", ylab="RF SMOTER Prediction")
abline(coef = c(0,1), col="blue")
abline(reg.true.pred.rf, col="red")
summary(reg.true.pred.rf)


# combine trainning data and data with missing value ##

train$rf.pred <- NA
impt.test <- rbind(train, miss.pred)

impt.test$impt.rf <- with(impt.test, ifelse(misstype=="MCAR"|misstype=="MAR"|misstype=="MNAR", impt.test$rf.pred, impt.test$score))


summary(impt.test$true)
t.test(impt.test$true, mu = 12, alternative = "two.sided")

summary(train$score)
t.test(train$score, mu = 12, alternative = "two.sided")

summary(impt.test$impt.rf)
t.test(impt.test$impt.rf, mu = 12, alternative = "two.sided")


############################## End of Test Group ##############################




############################################################################################
#####################################  Control Group  ######################################
############################################################################################



############################## End of Test Group ##############################


############################# Analysis of data with both groups combined ##############################




############################## End of Analysis ##############################


