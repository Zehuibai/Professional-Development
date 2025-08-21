################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Implementation of Multiple Imputation under MNAR  ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################



#install.packages("mice")
# mice: Multivariate Imputation by Chained Equations
library(mice)

#install.packages("plyr")
# plyr: Tools for Splitting, Applying and Combining Data
library(plyr)

#install.packages("lsmeans")
#Obtain least-squares means for linear, generalized linear, and mixed models. Compute contrasts or linear functions of
#  least-squares means, and comparisons of slopes.
library(lsmeans)



rm(list=ls())

### read in data
all=read.csv("C:\\Users\\206586\\OneDrive - BioNTech SE\\HH\\Development\\GMDS 2024\\Simulations\\sep24\\simu_score_all.csv")


### split by treatment group, as the MI is performed by treatment group (data distribution may be different between treatment groups)
### keep the relevant variables: outcome varialbe = score (with missing data), covariates = cov1-cov7 (no missing data)
test <- subset(all, (all$treat==1), select=c("score", "cov1","cov2","cov3","cov4","cov5","cov6","cov7"))
ctrl <- subset(all, (all$treat==2), select=c("score", "cov1","cov2","cov3","cov4","cov5","cov6","cov7"))



all.vars <- c("score", "cov1","cov2","cov3","cov4","cov5","cov6","cov7")
n.vars <- length(all.vars)
pred.matrix <- matrix(0, n.vars, n.vars,
                      dimnames = list(all.vars,all.vars))
pred.matrix["score",c("cov1","cov2","cov3","cov4","cov5","cov6","cov7")] <- 1 
# one could define a matrix that in each block one variable is dependent and others are predictors (in the cases where predictors with MD), 
#   but in this case, it doesn't matter, as only outcome variable with MD.
pred.matrix





################################ m imputations ############################################
### Method: 
# method = "mnar.norm"
# This function imputes data that are thought to be Missing Not at Random (MNAR) by the NARFCS method. 
#  *** NARFCS = not-at-random fully conditional specification ***
# The NARFCS procedure (Tompsett et al, 2018) generalises the so-called delta-adjustment
# sensitivity analysis method of Van Buuren, Boshuizen & Knook (1999) to the case with multiple
# incomplete variables within the FCS framework. In practical terms, the NARFCS procedure shifts the imputations drawn
# at each iteration of mice by a user-specified quantity that can vary across subjects, 
# to reflect systematic departures of the missing data from the data distribution imputed under MAR.


################################################# Test Group ###################################################
### penalize Test group by borrowing information from Control group, i.e., introducing a delta-adjustment ######

### ums = "unidentifiable model specification"
# A string containing the specification of the unidentifiable part of the imputation model, that is, the desired  
# delta-adjustment (offset) as a function of other variables and values for the corresponding deltas (sensitivity parameters).

m <- 50

mnar.blot <- list(score = list(ums = "-2"))
test.mnar <- mice(test, m=m,
                  predictorMatrix=pred.matrix,
                  method = "mnar.norm",
                  blots = mnar.blot, seed = 369)
# check the imputed values
test.mnar$imp$score  


############################################### Control Group ################################################
###############  Control Group still assume MAR ############### 
ctrl.mar <- mice(ctrl, m=m, 
                 predictorMatrix=pred.matrix,
                 method = "pmm",
                 seed=123)

# check the imputed values
ctrl.mar$imp$score      



### Get the imputed data sets 
#imp.test.mnar1 <- complete(test.mnar, action=1L)  # the first imputed data set only, action=1L is default 
imp.test.mnar <- complete(test.mnar, action="long")  # all imputed data sets stacked vertically, action="broad" produces a data set with where imputed data sets are stacked horizontally 
imp.test.mnar$treat <- 1   # re-assign treatment group

imp.ctrl.mar <- complete(ctrl.mar, action="long")
imp.ctrl.mar$treat <- 2   # re-assign treatment group

# combine the data from two group, which will be used for the analysis by imputation ID 
imp.mar <- rbind(imp.test.mnar, imp.ctrl.mar)



################################ Analysis of m imputed datasets separately ############################################

lsm.all <- data.frame() 
dif.all <- data.frame() 
for(i in 1:m){
  print(i)
  dat.in <- imp.mar[imp.mar$.imp==i,]
  lm <- lm(data=dat.in, score ~ treat+cov1)
  summary(lm)
  print(lm)
  lsm. <- lsmeans(lm, "treat")  # one should not save this directly as df, as this matrix is input of contrast function 
  lsm <- as.data.frame(lsm.)
  dif. <- contrast(lsm., "pairwise")
  dif <- as.data.frame(confint(dif.))
  lsm$impid <- i
  dif$impid <- i
  lsm.all <- rbind(lsm.all, lsm)
  dif.all <- rbind(dif.all, dif)
}


################################ Pooling the results together ############################################
#### Required reporting matrices: LS means per group and their 95% CI, difference in LS means and its 95% CI, 

### Pooling the LS means per group and get their 95% CIs

# Test group
dat1 <- lsm.all[lsm.all$treat==1,]
pool1 <- pool.scalar(Q=dat1[["lsmean"]], U=dat1[["SE"]], rule = c("rubin1987"))
pool1$t.stats <- qt(p = .975, df = pool1$df)
pool1$ci.lower <- pool1$qbar - pool1$t.stats*pool1$ubar
pool1$ci.upper <- pool1$qbar + pool1$t.stats*pool1$ubar
print(pool1$qbar)
print(pool1$ubar)
print(pool1$ci.lower)
print(pool1$ci.upper)


# Control group
dat2 <- lsm.all[lsm.all$treat==2,]
pool2 <- pool.scalar(Q=dat2[["lsmean"]], U=dat2[["SE"]], rule = c("rubin1987"))
pool2$t.stats <- qt(p = .975, df = pool2$df)
pool2$ci.lower <- pool2$qbar - pool2$t.stats*pool2$ubar
pool2$ci.upper <- pool2$qbar + pool2$t.stats*pool2$ubar
print(pool2$qbar)
print(pool2$ubar)
print(pool2$ci.lower)
print(pool2$ci.upper)

### Pooling the difference and get its 95% CI 
pool.dif <- pool.scalar(Q=dif.all[["estimate"]], U=dif.all[["SE"]], rule = c("rubin1987"))
pool.dif$t.stats <- qt(p = .975, df = pool.dif$df)
pool.dif$ci.lower <- pool.dif$qbar - pool.dif$t.stats*pool.dif$ubar
pool.dif$ci.upper <- pool.dif$qbar + pool.dif$t.stats*pool.dif$ubar
print(pool.dif$qbar)
print(pool.dif$ubar)
print(pool.dif$ci.lower)
print(pool.dif$ci.upper)


########################################### End of MNAR imputation ###########################################

