################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Implementation of Multiple Imputation under MNAR  ##########
########## Author:                                   ##########
########## Date:                                          ##########
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



### Get the imputed data sets 


# combine the data from two group, which will be used for the analysis by imputation ID 




################################ Analysis of m imputed datasets separately ############################################



################################ Pooling the results together ############################################
#### Required reporting matrices: LS means per group and their 95% CI, difference in LS means and its 95% CI, 

### Pooling the LS means per group and get their 95% CIs

# Test group



# Control group


### Pooling the difference and get its 95% CI 



########################################### End of MNAR imputation ###########################################

