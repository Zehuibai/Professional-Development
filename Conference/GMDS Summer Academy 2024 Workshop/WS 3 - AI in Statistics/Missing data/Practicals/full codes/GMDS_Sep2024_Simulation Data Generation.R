################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Generate correlated simulation data               ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################



#install.packages("BinNor")
#install.packages("dplyr")
#help(package = "dplyr")
#install.packages("tidyverse")
#install.packages("ggplot2")

library(BinNor)
library(dplyr)
library(tidyverse)
library(ggplot2)

rm(list = ls())

# MCAR flag is independent from the outcome and the co-variates (Corr.Coeff.=0.02)
# MAR flag is correlated with the first co-variate only (Corr.Coeff.=0.4) and independent from the outcome and the other covariates 
# MNAR low flag is negative correlated with outcome (high level correlation) and negative correlated with the 2-7th co-variate (medium level correlation)
# MNAR high flag is positive correlated with outcome (high level correlation) and the 2-7th co-variate (medium level correlation)
# The first co-variate is correlated with MAR flag only
# The 2-7 co-variates are correlated with the outcome (high level correlation), therefore they are also correlated with each other (medium level correlation)


cc.low = 0.05
mnarl.cov =  -0.3
mnarl.out =  -0.55
#cov.out = 0.55
cov.out = 0.25
cov.cov = 0.25

no.rows=300
no.bin=2
no.nor=8
mean.vec.nor1=c(12,10,10,10,10,10,10,10)   # test group target variable mean=10, covariates mean=10
mean.vec.nor2=c(10,10,10,10,10,10,10,10)   # control group target variable mean=8, covariates mean=8
var.nor=c(10,10,10,10,10,10,10,10)


corr.vec=c(cc.low,   cc.low,   cc.low,   cc.low,   cc.low,   cc.low,   cc.low,   cc.low,   cc.low,
           mnarl.out,cc.low,   mnarl.cov,mnarl.cov,mnarl.cov,mnarl.cov,mnarl.cov,mnarl.cov,
           cc.low,   cov.out,  cov.out,  cov.out,  cov.out,  cov.out,  cov.out,
           cc.low,   cc.low,   cc.low,   cc.low,   cc.low,   cc.low,  
           cov.cov,  cov.cov,  cov.cov,  cov.cov,  cov.cov, 
           cov.cov,  cov.cov,  cov.cov,  cov.cov,
           cov.cov,  cov.cov,  cov.cov,
           cov.cov,  cov.cov,
           cov.cov);


prop.vec.bin1=c(0.10,0.30)  # Test Group: actual total will be between 30% and 40% due to overlapping
prop.vec.bin2=c(0.15,0.10)  # ctrlerence Group: actual total will be 15% MAR, the 10% MNAR will be re-set to be non-missing


simu.num=1

# generate TEST group data 
mytotdf <- data.frame()  # the first (blank) data set to be appended
for (i in 1:simu.num){
  set.seed(i)
  cmat = lower.tri.to.corr.mat(corr.vec,10)
  sigma.star=compute.sigma.star(no.bin=no.bin, no.nor=no.nor, prop.vec.bin=prop.vec.bin1, corr.mat=cmat)
  
  mydata=jointly.generate.binary.normal(no.rows, no.bin, no.nor, prop.vec.bin1,
                                        mean.vec.nor1, var.nor, sigma_star=sigma.star$sigma_star,
                                        continue.with.warning=TRUE)
  mydf <- as.data.frame(mydata)
  # add simulation id
  mydf$simu.id <- i
  # append data sets
  mytotdf <- rbind(mytotdf, mydf)
}
test.grp <-mytotdf 
test.grp$treat <- 1




# generate control group data 
mytotdf <- data.frame()  # the first (blank) data set to be appended
for (i in 1:simu.num){
  set.seed(i)
  cmat = lower.tri.to.corr.mat(corr.vec,10)
  sigma.star=compute.sigma.star(no.bin=no.bin, no.nor=no.nor, prop.vec.bin=prop.vec.bin2, corr.mat=cmat)
  
  mydata=jointly.generate.binary.normal(no.rows, no.bin, no.nor, prop.vec.bin2,
                                        mean.vec.nor2, var.nor, sigma_star=sigma.star$sigma_star,
                                        continue.with.warning=TRUE)
  mydf <- as.data.frame(mydata)
  # add simulation id
  mydf$simu.id <- i
  # append data sets
  mytotdf <- rbind(mytotdf, mydf)
}
ctrl.grp <-mytotdf 
ctrl.grp$treat <- 2
ctrl.grp$V2 <- 0   # No MNAR in control group


# append two groups 
df <- rbind(test.grp, ctrl.grp)
# add Subject ID variable
df$id <- seq.int(nrow(df))

## Rounding outcome variable to 2 dp 
df <- df %>% mutate(V3 = round(V3, 2)) 



names(df)[1] <- "marfl"
names(df)[2] <- "mnarfl"
names(df)[3] <- "true"
names(df)[4] <- "cov1"
names(df)[5] <- "cov2"
names(df)[6] <- "cov3"
names(df)[7] <- "cov4"
names(df)[8] <- "cov5"
names(df)[9] <- "cov6"
names(df)[10] <- "cov7"


aggregate(x = df$true,               # Specify data column
          by = list(df$treat),       # Specify group indicator
          FUN = mean)                # Specify function (i.e. mean)


### there should be no/very less negative value given the distribution 
### But to cope with the bug in UtilOptimRegress fuction, set all negative values to 0 
df$true <- pmax(df$true, 0)
df$cov1 <- pmax(df$cov1, 0)
df$cov2 <- pmax(df$cov2, 0)
df$cov3 <- pmax(df$cov3, 0)
df$cov4 <- pmax(df$cov4, 0)
df$cov5 <- pmax(df$cov5, 0)
df$cov6 <- pmax(df$cov6, 0)
df$cov7 <- pmax(df$cov7, 0)


aggregate(x = df$true,               # Specify data column
          by = list(df$treat),       # Specify group indicator
          FUN = mean)                # Specify function (i.e. mean)




#attach(df)

# check MD proportions   
table(df$treat, df$mnarfl)
table(df$treat, df$marfl)

# check flags   
table(df$mnarfl)
table(df$marfl)

# If MNAR data point is >=10 (on the right tail), then it is not MNAR
df <- within(df, mnarfl[mnarfl == 1 & true >= 12] <- 0)

# check to make sure no redundant missingness 
table(df$marfl, df$mnarfl)


# If one record is both MNAR and MAR, then it is MAR only
df <- within(df, mnarfl[mnarfl == 1 & marfl == 1] <- 0)


# check to make sure no redundant missingness 
table(df$marfl, df$mnarfl)


# generate outcome variable score which with missing values
df$score <- with(df, ifelse(mnarfl==1 | marfl==1, NA, df$true))


# check MD proportions   
table(df$treat, df$mnarfl)
table(df$treat, df$marfl)

# generate one category variable for missingness

df<-df%>%mutate(misstypen = case_when(
  marfl  ==1 ~ 1,
  mnarfl ==1 ~ 2,
  marfl==0 & mnarfl==0 ~ 9
))
df<-df%>%mutate(misstype = case_when(
  marfl  ==1 ~ 'MAR',
  mnarfl ==1 ~ 'MNAR',
  marfl==0 & mnarfl==0 ~ 'Non-Missing'
))

table(df$misstype)

tapply(df$true, df$treat, summary)
tapply(df$score, df$treat, summary)


write.csv(df, file.path("C:\\Users\\206586\\OneDrive - BioNTech SE\\HH\\Development\\GMDS 2024\\Simulations\\sep24\\simu_score_all.csv"))

########################## End of data generation ###############################
