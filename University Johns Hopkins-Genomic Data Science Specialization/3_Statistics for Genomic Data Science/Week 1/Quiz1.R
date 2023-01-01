# Title     : Quiz 1
# Objective : To run code for the 1st week quiz of Statisticis for Genomics
# Created by: valengo
# Created on: 05/05/21
if (!requireNamespace("../renv", quietly = TRUE))
  install.packages("../renv")

# renv::restore()
renv::init(bare = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("Biobase", "GenomicRanges", "SummarizedExperiment", "dplyr", "DESeq2", "rafalib", "rmarkdown")
BiocManager::install(packages, update = TRUE, ask = FALSE)

renv::snapshot()

Biobase::data(sample.ExpressionSet, package = "Biobase")
se <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)
genomic_table <- SummarizedExperiment::assay(se)
pheno <- SummarizedExperiment::colData(se)
features <- SummarizedExperiment::rowData(se)
SummarizedExperiment::rowRanges(se)

# Load the Bottomly and the Bodymap data sets with the following code
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot <- bottomly.eset
pdata_bot <- Biobase::pData(bot)
table(pdata_bot$strain)
table(pdata_bot$strain, pdata_bot$experiment.number)
table(pdata_bot$strain, pdata_bot$num.tech.reps)
table(pdata_bot$num.tech.reps)

con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
pdata_bm <- Biobase::pData(bm)
table(pdata_bm$gender, pdata_bm$race)
table(pdata_bm$gender, pdata_bm$tissue.type)
table(pdata_bm$num.tech.reps)

# Load the Bottomly data
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset

edata <- Biobase::exprs(bm)
row_sums <- rowSums(edata)
edata <- edata[order(-row_sums),]
index <- 1:500
heatmap(edata[index,],Rowv=NA,Colv=NA) # 16864261,  12245926

edata <- Biobase::exprs(bm)
row_sums < rowSums(edata)
index <- which(rank(row_sums) < 500 )
heatmap(edata[index,],Colv=NA)

edata <- Biobase::exprs(bm)
row_sums <- rowSums(edata)
edata <- edata[order(row_sums),]
index <- which(rank(-row_sums) < 500 )
heatmap(edata[index,],Rowv=NA,Colv=NA) # 0, 0, 0, 0

edata <- Biobase::exprs(bm)
row_sums <- rowSums(edata)
index <- which(rank(-row_sums) < 500 )
heatmap(edata[index,],Rowv=NA,Colv=NA) # 307296, 202806, 319412

# Load the Bodymap data using the following code:
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm <- bodymap.eset
pdata <- Biobase::pData(bm)
edata <- Biobase::exprs(bm)

mm <- log2(edata[,1]+1) - log2(edata[,2]+1)
aa <- log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

edata_rlog <-DESeq2::rlog(edata)
mm <- (edata_rlog[,1]) - (edata_rlog[,2])
aa <- (edata_rlog[,1]) + (edata_rlog[,2])
plot(aa,mm,col=2)

# Load the Montgomery and Pickrell eSet:
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

# Clustering with no changes to the data
dist1 <- dist(t(edata)) # by default calculates the distance between rows
hclust1 <- hclust(dist1)
plot(hclust1, hang=-1) # plot dendogram
rafalib::myplclust(hclust1, lab.col = as.numeric(pdata$study))


# After filtering all genes with rowMeans less than 100
filtered <- edata[rowMeans(edata) > 99,]
dist1 <- dist(t(filtered))
hclust1 <- hclust(dist1)
plot(hclust1, hang=-1) # plot dendogram
rafalib::myplclust(hclust1, lab.col = as.numeric(pdata$study))

# After taking the log2 transform of the data without filtering
log_data <- log2(edata + 1)
dist1 <- dist(t(log_data))
hclust1 <- hclust(dist1)
plot(hclust1, hang=-1) # plot dendogram
rafalib::myplclust(hclust1, lab.col = as.numeric(pdata$study))

# Load the Montgomery and Pickrell eSet:
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp <- montpick.eset
pdata <- Biobase::pData(mp)
edata <- as.data.frame(Biobase::exprs(mp))
fdata <- Biobase::fData(mp)

# Cluster the samples using k-means clustering after applying the log transform
log_edata <- log2(edata + 1)
# Set a seed for reproducible results
set.seed(1235)
# Choose two clusters
kmeans1 <- kmeans(t(log_edata),centers=2)
matplot(t(kmeans1$centers),col=1:2,type="l",lwd=2)
table(kmeans1$cluster)

# cuttree
dist1 <- dist(t(log_edata)) # by default calculates the distance between rows
hclust1 <- hclust(dist1)
cutree1 <- cutree(hclust1, k = 2)
table(cutree1)