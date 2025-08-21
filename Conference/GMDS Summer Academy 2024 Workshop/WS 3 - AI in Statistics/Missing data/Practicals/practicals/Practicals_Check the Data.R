################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Check and understand the simulated data           ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################

#install.packages("plyr")
# plyr: Tools for Splitting, Applying and Combining Data
library(plyr)

#install.packages("lsmeans")
#Obtain least-squares means for linear, generalized linear, and mixed models. Compute contrasts or linear functions of
#  least-squares means, and comparisons of slopes.
library(lsmeans)
library(ggplot2)


rm(list=ls())

### read in data
all=read.csv("...\\simu_score_all.csv")


tapply(all$true, all$treat, summary)
tapply(all$score, all$treat, summary)

freq <- table(all$treat, all$misstype)
freq
100*freq/300 


test <- subset(all, (all$treat==1), )
ctrl <- subset(all, (all$treat==2), )

sd(test$true)  
sd(test$score, na.rm = TRUE)  ### MNAR also leads to underestimation of variability ###

sd(ctrl$true)
sd(ctrl$score, na.rm = TRUE)


ggplot(data = test, aes(true, id, color = misstype), coord_cartesian(xlim = c(0, 20))) +
  geom_point() +
  scale_color_manual(values = c("MNAR"="red", "MCAR"="green", "MAR"="blue", "Non-Missing"="black")) +
  labs(y= "Patient ID", x = "Outcome Variable")

ggplot(data = ctrl, aes(true, id, color = misstype), coord_cartesian(xlim = c(0, 20))) +
  geom_point() +
  scale_color_manual(values = c("MNAR"="red", "MCAR"="green", "MAR"="blue", "Non-Missing"="black")) +
  labs(y= "Patient ID", x = "Outcome Variable")




lm <- lm(data=all, true ~ treat+cov1)
summary(lm)
lsm. <- lsmeans(lm, "treat")  # one should not save this directly as df, as this matrix is input of contrast function 
lsm <- as.data.frame(lsm.)
print(lsm)
dif. <- contrast(lsm., "pairwise")
dif <- as.data.frame(confint(dif.))
print(dif)


lm <- lm(data=all, score ~ treat+cov1)
summary(lm)
lsm. <- lsmeans(lm, "treat")   
lsm <- as.data.frame(lsm.)
print(lsm)
dif. <- contrast(lsm., "pairwise")
dif <- as.data.frame(confint(dif.))
print(dif)



