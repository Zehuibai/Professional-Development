################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Implementation of Multiple Imputation under MAR   ##########
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
pred.matrix
pred.matrix["score",c("cov1","cov2","cov3","cov4","cov5","cov6","cov7")] <- 1 
  # one could define a matrix that in each block one variable is dependent and others are predictors (in the cases where predictors with MD), 
  #   but in this case, it doesn't matter, as only outcome variable with MD.
pred.matrix





################################ m imputations ############################################
### Method: 
#   Drawing imputations from the observed values by predictive-mean matching (PMM):
#     For a given subject with missing data on the variable in question, PMM identifies those subjects with no missing data on the variable
#     in question whose linear predictors (created using the regression coefficients from the fitted imputation model) are close to the 
#     linear predictor of the given subject (created using the regression coefficients sampled from the appropriate posterior distribution, 
#     as described above). Of those subjects who are close, one subject is selected at random and the observed value of the given variable for
#     that randomly selected subject is used as the imputed value of the variable for the subject with missing data.


m <- 50

test.mar <- mice(test, m=m,
            predictorMatrix=pred.matrix,
            method = "pmm",
            seed=123)

#plot(test.mar)

# check the imputed values
test.mar$imp$score      


ctrl.mar <- mice(ctrl, m=m, 
                predictorMatrix=pred.matrix,
                method = "pmm",
                seed=123)

# check the imputed values
ctrl.mar$imp$score      


### Get the imputed data sets 
#imp.test.mar1 <- complete(test.mar, action=1L)  # the first imputed data set only, action=1L is default 
imp.test.mar <- complete(test.mar, action="long")  # all imputed data sets stacked vertically, action="broad" produces a data set with where imputed data sets are stacked horizontally 
imp.test.mar$treat <- 1   # re-assign treatment group

imp.ctrl.mar <- complete(ctrl.mar, action="long")
imp.ctrl.mar$treat <- 2   # re-assign treatment group

# combine the data from two group, which will be used for the analysis by imputation ID 
imp.mar <- rbind(imp.test.mar, imp.ctrl.mar)



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


### just to check the averaging of point estimate and SE - they ARE the same!                     
#aggregate(x = lsm.all[["lsmean"]], by = list(lsm.all$treat), FUN = mean)
#aggregate(x = lsm.all[["SE"]], by = list(lsm.all$treat), FUN = mean)


################################ MNAR ############################################

