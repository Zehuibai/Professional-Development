**Authors:** Hugo Tavares, Georg Zeller

----

These are extra materials used as a complement to 
[Data Carpentry in R](http://www.datacarpentry.org/R-ecology-lesson/) 
courses, and thus assume that some of those lessons were covered beforehand. 

These lessons are under active development and may change over time. 

The lessons are modular so can be taught in different order than shown here 
(apart from the introduction, which should always be the first):

* [Introduction](01_rnaseq_intro.html) to the dataset.
* [Basic exploratory analysis](02_rnaseq_exploratory.html)
to understand some properties of expression data.
* [Using principal component analysis (PCA)](03_rnaseq_pca.html)
to explore transcriptome-wide effects of our experimental design.
* Exploring gene expression patterns:
    * [Identifying candidate genes](04a_explore_test_results.html) from a differential analysis test.
    * [Using clustering](04b_rnaseq_clustering.html) to partition genes into groups.


### Important note

There are many dedicated packages to deal with RNAseq data, mostly 
within the [Bioconductor](https://bioconductor.org/) package repository. 

**This lesson is not about analysing RNAseq data** (that would be a topic for a whole 
course!), but rather to show you how the data manipulation principles learned 
so far can be applied to explore these kind of data. 

If you are doing RNAseq analysis, you should use 
[dedicated packages and workflows](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html), 
which implement models to account for particular features of these data.

If you are interested, you can see [how the data for this lesson was pre-processed](https://github.com/tavareshugo/data-carpentry-rnaseq/blob/master/prepare_fission_data.R) 
using the `DESeq2` package.

