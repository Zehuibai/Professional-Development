# Title     : Quiz 4
# Objective : To run code for the 4th week quiz of Statisticis for Genomics
# Created by: valengo
# Created on: 31/05/21

packages <- c("goseq", "MatrixEQTL", "org.Mm.eg.db")
BiocManager::install(packages, update = TRUE, ask = FALSE)
renv::snapshot()

sup <- geneLenDataBase::supportedGenomes()

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot <- bottomly.eset
pdata_bot <- Biobase::pData(bot)
fdata_bot <- Biobase::featureData(bot)
edata <- Biobase::exprs(bot)
fdata_bot <- fdata_bot[rowMeans(edata) > 5]
edata <- edata[rowMeans(edata) > 5, ]
edata <- log2(edata+1)

# perform a differential expression analysis using limma with only the strain variable as an outcome
mod <- model.matrix(~ pdata_bot$strain)
fit_limma <- limma::lmFit(edata,mod)
ebayes_limma <- limma::eBayes(fit_limma)
top <- limma::topTable(ebayes_limma, number=dim(edata)[1], sort.by="none")

top_bh <- p.adjust(top$P.Value, method = "BH")
sum(top_bh < 0.05)

# ---------------------------------------------------------------- #
genes <- as.integer(top$adj.P.Val < 0.05)
not_na <- !is.na(genes)
names(genes) <- rownames(edata)
genes <- genes[not_na]

pwf <-goseq::nullp(genes,"mm9","ensGene")
head(pwf)

GO.wall <- goseq::goseq(pwf,"mm9","ensGene")
head(GO.wall)

# ---------------------------------------------------------------- #
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot <- bottomly.eset
pdata_bot <- Biobase::pData(bot)
fdata_bot <- Biobase::featureData(bot)
edata <- Biobase::exprs(bot)
fdata_bot <- fdata_bot[rowMeans(edata) > 5]
edata <- edata[rowMeans(edata) > 5, ]
edata <- log2(edata+1)

mod_adj <- model.matrix(~ pdata_bot$strain + as.factor(pdata_bot$lane.number))
fit_adj <- limma::lmFit(edata, mod_adj)
ebayes_limma_adj <- limma::eBayes(fit_adj)
top_adj <- limma::topTable(ebayes_limma_adj, number=dim(edata)[1], sort.by="none")

genes_adj <- as.integer(top_adj$adj.P.Val < 0.05)
not_na <- !is.na(genes_adj)
names(genes_adj) <- rownames(edata)
genes_adj <- genes_adj[not_na]

pwf_adj <-goseq::nullp(genes_adj,"mm9","ensGene")
head(pwf)

GO.wall_adj <- goseq::goseq(pwf_adj,"mm9","ensGene")
