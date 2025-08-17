# Title     : Quiz 2
# Objective : To run code for the 2nd week quiz of Statisticis for Genomics
# Created by: valengo
# Created on: 14/05/21

packages <- c("devtools", "broom", "limma", "edge")
BiocManager::install(packages, update = TRUE, ask = FALSE)
renv::snapshot()

con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

# What percentage of variation is explained by the 1st principal component in the data set if you:
# 1. Do no transformations?
svd1 <- svd(edata)
svd1_explained <- svd1$d^2/sum(svd1$d^2)
round(svd1_explained[1], 2)

# 2. log2(data + 1) transform?
edata_log2 <- log2(edata + 1)
svd2 <- svd(edata_log2)
svd2_explained <- svd2$d^2/sum(svd2$d^2)
round(svd2_explained[1], 2)

# 3. log2(data + 1) transform and subtract row means?
edata_log2 <- log2(edata + 1)
edata_centered <- edata_log2 - rowMeans(edata_log2)
svd3 <- svd(edata_centered)
svd3_explained <- svd3$d^2/sum(svd3$d^2)
round(svd3_explained[1], 2)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

log_edata <- log2(edata + 1)
edata_centered <- log_edata - rowMeans(log_edata)
set.seed(333)
kmeans1 <- kmeans(t(log_edata), centers=2)

svd1 <- svd(edata_centered)
# What is the correlation between the first singular vector and the sample clustering indicator?
round(cor(svd1$v[, 1], kmeans1$cluster), 2)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
edata <- Biobase::exprs(bm)
pdata_bm <- Biobase::pData(bm)

# Fit a linear model relating the first gene’s counts to the number of technical replicates,
# treating the number of replicates as a factor
lm1 <- lm(edata[1,] ~ as.factor(pdata_bm$num.tech.reps))
broom::tidy(lm1)
# Plot the data for this gene versus the covariate
plot(edata[1, ], pdata_bm$num.tech.reps, col=1)
abline(lm1$coeff[1],lm1$coeff[2], col=2,lwd=3)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
edata <- Biobase::exprs(bm)
pdata_bm <- Biobase::pData(bm)

# Fit a linear model relating the first gene’s counts to the age of the person and the sex of the samples
lm1 <- lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
broom::tidy(lm1)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

log_edata <- log2(edata + 1)
# Fit a regression model to each sample using population as the outcome
mod <- model.matrix(~ pdata$population)
fit <- lm.fit(mod, t(edata))
dim(fit$residuals)
dim(fit$effects)
dim(fit$coefficients)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

log_edata <- log2(edata + 1)
# Fit a regression model to each sample using population as the outcome
mod <- model.matrix(~ pdata$population)
fit <- lm.fit(mod, t(edata))
fit$effects[1:5, 1:3]
dim(fit$effects)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
edata <- Biobase::exprs(bm)
pdata_bm <- Biobase::pData(bm)

# Fit many regression models to the expression data where age is the outcome variable using the lmFit
table(pdata_bm$age,useNA="ifany")
pdata_bm_sub <- subset(pdata_bm,(!is.na(pdata_bm$age)))
edata_bm_sub <- as.matrix(edata[, row.names(pdata_bm_sub)])

mod <- model.matrix(~ pdata_bm$age)
fit <- limma::lmFit(edata_bm_sub, mod)
fit$coefficients[1000,]
plot(pdata_bm_sub$age, edata_bm_sub[1000,])
abline(fit$coefficients[1000,], col=2, lwd=3)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
edata <- Biobase::exprs(bm)
pdata_bm <- Biobase::pData(bm)

pdata_bm_sub <- subset(pdata_bm,(!is.na(pdata_bm$age)))
edata_bm_sub <- as.matrix(edata[, row.names(pdata_bm_sub)])

# Fit many regression models to the expression data where age is the outcome variable and
# tissue.type is an adjustment variable
mod_adj <- model.matrix(~pdata_bm_sub$age + pdata_bm_sub$tissue.type)
fit <- limma::lmFit(edata_bm_sub, mod_adj)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
edata <- Biobase::exprs(bm)
pdata_bm <- Biobase::pData(bm)

set.seed(33353)
pdata_bm_sub <- subset(pdata_bm,(!is.na(pdata_bm$age)))
edata_bm_sub <- as.matrix(edata[, row.names(pdata_bm_sub)])
log_edata <- log2(edata_bm_sub + 1)
pre_edata <- log_edata[rowMeans(log_edata) > 1, ]

mod <- model.matrix(~pdata_bm_sub$age, data=pdata_bm_sub)
mod0 <- model.matrix(~1, data=pdata_bm_sub)
sva1 <- sva::sva(pre_edata, mod, mod0, n.sv=1)
round(cor(sva1$sv, pdata_bm_sub$age), 2)
round(cor(sva1$sv, as.numeric(pdata_bm_sub$gender)), 2)
round(cor(sva1$sv, as.numeric(pdata_bm_sub$race)), 2)

