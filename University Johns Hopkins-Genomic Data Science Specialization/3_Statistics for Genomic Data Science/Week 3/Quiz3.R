# Title     : Quiz 3
# Objective : To run code for the 3rd week quiz of Statisticis for Genomics
# Created by: valengo
# Created on: 16/05/21

packages <- c("snpStats", "MASS", "DESeq2")
BiocManager::install(packages, update = TRUE, ask = FALSE)
renv::snapshot()

library(broom)

library(snpStats)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata <- sub.10@.Data
status <- subject.support$cc

# Fit a linear model and a logistic regression model to the data for the 3rd SNP
snp3 <- as.numeric(snpdata[,3])
snp3[snp3==0] <- NA

lm1 <- lm(status ~ snp3)
broom::tidy(lm1)

glm1 <- glm(status ~ snp3, family="binomial")
broom::tidy(glm1)

# ---------------------------------------------------------------- #
# convention: 0 = missing, 1 = "A/A", 2 = "A/B" or "B/A" and 3 = "B/B"
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata <- sub.10@.Data
status <- subject.support$cc

snp10 <- as.numeric(snpdata[,10])
snp10[snp10==0] <- NA

# Fit a logistic regression model on a recessive (need 2 copies of minor
# allele to confer risk: 0 - 0 - 2) and additive scale (0 - 1 - 2) for the 10th SNP
adtv_glm <- glm(status ~ snp10, family = "binomial")

snp10_rec <- (snp10 == 3)
rec_glm <- glm(status ~ snp10_rec, family = "binomial")

# Make a table of the fitted values versus the case/control status
broom::tidy(adtv_glm)
broom::tidy(rec_glm)

# ---------------------------------------------------------------- #
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata <- sub.10@.Data
status <- subject.support$cc

# Fit an additive logistic regression model to each SNP
n_rows <- dim(snpdata)[2]
effects <- matrix(NA, nrow=n_rows, ncol=1)

for(i in 1:n_rows) {
  snp <- as.numeric(snpdata[,i])
  snp[snp==0] <- NA
  effects[i,] <- broom::tidy(glm(status ~ snp,family="binomial"))$statistic[2]
}
# What is the average effect size (z-value of your logistic regression)?
# What is the max? What is the minimum?
mean(effects)
min(effects)
max(effects)

# ---------------------------------------------------------------- #
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata <- sub.10@.Data
status <- subject.support$cc

# Fit an additive logistic regression model to each SNP and square de coefficients
n_rows <- dim(snpdata)[2]
effects <- matrix(NA, nrow=n_rows, ncol=1)

for(i in 1:n_rows) {
  snp <- as.numeric(snpdata[,i])
  snp[snp==0] <- NA
  effects[i,] <- broom::tidy(glm(status ~ snp,family="binomial"))$statistic[2]
}
effects_coeff <- effects ^ 2

glm_all_adj <- snp.rhs.tests(status ~ 1,snp.data=sub.10)
chi_squared <- chi.squared(glm_all_adj)

cor(effects_coeff, chi_squared)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

edata <- log2(as.matrix(edata) + 1)
fstats <- genefilter::rowFtests(edata, pdata$study)
tstats <- genefilter::rowttests(edata, pdata$study)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
edata <- edata[rowMeans(edata) > 100,]
fdata <- Biobase::fData(mp)

# First test for differences between the studies using the DESeq2 package using the DESeq function
dds <- DESeq2::DESeqDataSetFromMatrix(countData = edata,
                                      colData = pdata,
                                      design = ~ study)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)
colnames(res)

edata <- log2(as.matrix(edata) + 1)
mod <- model.matrix(~ pdata$study)
fit_limma <- limma::lmFit(edata,mod)
ebayes_limma <- limma::eBayes(fit_limma)
top <- limma::topTable(ebayes_limma, number=dim(edata)[1], sort.by="none")
top <- top[match(rownames(res), rownames(top)),]
colnames(top)

round(cor(res$stat, top$t), 2)
DESeq2::plotMA(res)

par(mfcol=c(1,2))
plot(log2(res$baseMean+1), res$log2FoldChange, ylim=c(-5,5), pch=16)
plot(top$AveExpr, top$logFC, ylim=c(-5,5),, pch=16)

# ---------------------------------------------------------------- #
# Apply the Benjamni-Hochberg correction to the P-values from the two previous analyses
res_bh <- p.adjust(res$pvalue, method = "BH")
sum(res_bh < 0.05)

top_bh <- p.adjust(top$P.Value, method = "BH")
sum(top_bh < 0.05)