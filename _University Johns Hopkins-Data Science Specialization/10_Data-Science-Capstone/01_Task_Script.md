Task 1: Getting and Cleaning the Data
================
Mark Blackmore
2017-11-27

-   [Download and explore the data](#download-and-explore-the-data)
-   [Sample the data and save the sample](#sample-the-data-and-save-the-sample)
-   [Clean the sample data](#clean-the-sample-data)
-   [Initial Exploratory Data Analysis](#initial-exploratory-data-analysis)
-   [Session info](#session-info)

Download and explore the data
-----------------------------

Create a data directory

``` r
if (!file.exists("data")) {
  dir.create("data")
}
```

Download the data

``` r
# url <- "https://d396qusza40orc.cloudfront.net/dsscapstone/dataset/Coursera-SwiftKey.zip"
# download(url, dest="dataset.zip", mode="wb") 
# unzip ("dataset.zip", exdir = "./data")
```

English Repository Files

``` r
blogs_file   <- "./data/final/en_US/en_US.blogs.txt"
news_file    <- "./data/final/en_US/en_US.news.txt"
twitter_file <- "./data/final/en_US/en_US.twitter.txt"  
```

File Sizes (Mb)

``` r
blogs_size   <- file.size(blogs_file) / (2^20)
news_size    <- file.size(news_file) / (2^20)
twitter_size <- file.size(twitter_file) / (2^20)
```

Read the data files

``` r
blogs   <- read_lines(blogs_file)
news    <- read_lines(news_file)
twitter <- read_lines(twitter_file) 
```

Number of Lines per file

``` r
blogs_lines   <- length(blogs)
news_lines    <- length(news)
twitter_lines <- length(twitter)
total_lines   <- blogs_lines + news_lines + twitter_lines
```

Distibution of characters per line, by file

``` r
blogs_nchar   <- nchar(blogs)
news_nchar    <- nchar(news)
twitter_nchar <- nchar(twitter)

boxplot(blogs_nchar, news_nchar, twitter_nchar, log = "y",
        names = c("blogs", "news", "twitter"),
        ylab = "log(Number of Characters)", xlab = "File Name") 
        title("Comparing Distributions of Chracters per Line")
```

![](01_Task_Script_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

Total characters per file

``` r
blogs_nchar_sum   <- sum(blogs_nchar)
news_nchar_sum    <- sum(news_nchar)
twitter_nchar_sum <- sum(twitter_nchar)
```

Total words per file

``` r
blogs_words <- wordcount(blogs, sep = " ")
news_words  <- wordcount(news,  sep = " ")
twitter_words <- wordcount(news, sep = " ")
```

Create summary of repo stats

``` r
repo_summary <- data.frame(f_names = c("blogs", "news", "twitter"),
                           f_size  = c(blogs_size, news_size, twitter_size),
                           f_lines = c(blogs_lines, news_lines, twitter_lines),
                           n_char =  c(blogs_nchar_sum, news_nchar_sum, twitter_nchar_sum),
                           n_words = c(blogs_words, news_words, twitter_words))
repo_summary <- repo_summary %>% mutate(pct_n_char = round(n_char/sum(n_char), 2))
repo_summary <- repo_summary %>% mutate(pct_lines = round(f_lines/sum(f_lines), 2))
repo_summary <- repo_summary %>% mutate(pct_words = round(n_words/sum(n_words), 2))
kable(repo_summary)
```

| f\_names |   f\_size|  f\_lines|    n\_char|  n\_words|  pct\_n\_char|  pct\_lines|  pct\_words|
|:---------|---------:|---------:|----------:|---------:|-------------:|-----------:|-----------:|
| blogs    |  200.4242|    899288|  206824505|  37334131|          0.36|        0.21|        0.35|
| news     |  196.2775|   1010242|  203223154|  34372528|          0.36|        0.24|        0.32|
| twitter  |  159.3641|   2360148|  162096031|  34372528|          0.28|        0.55|        0.32|

Sample the data and save the sample
-----------------------------------

Compute sample sizes in terms of lines

``` r
sample_pct = 0.05
set.seed(1001)
blogs_size   <- blogs_lines * sample_pct
news_size    <- news_lines * sample_pct
twitter_size <- twitter_lines * sample_pct
```

Create samples

``` r
blogs_sample   <- sample(blogs, blogs_size)
news_sample    <- sample(news, news_size)
twitter_sample <- sample(twitter, twitter_size)
repo_sample    <- c(blogs_sample, news_sample, twitter_sample)
```

Save sample

``` r
writeLines(repo_sample, "./data/final/en_US/en_US.repo_sample.txt")
saveRDS(repo_sample, file = "./data/final/en_US/repo_sample.rds" )
```

Clean the sample data
---------------------

Use `tm` to create and clean the corpus

``` r
clean_sample <- Corpus(VectorSource(repo_sample))
print(as.character(clean_sample[[1]]))
```

    ## [1] "Love the use of onomatopoeia, and I wish they made dumplings with just scallions and cabbage. But there are plenty of places on 8th where you can buy dumplings (veggie or otherwise) for 4 for a dollar or so. Several places that look like they don't sell anything (just a white booth with a see-through window) actually sell delicious stuff to take home a cook. It's inexpensive and usually quite good. Enjoy."

Remove URL's
Source: [R and Data Mining](%22http://www.rdatamining.com/books/rdm/faq/removeurlsfromtext%22)

``` r
removeURL <- function(x) gsub("http[^[:space:]]*", "", x)
clean_sample <- tm_map(clean_sample, content_transformer(removeURL))

# Remove anything other than English letters or space
removeNumPunct <- function(x) gsub("[^[:alpha:][:space:]]*", "", x)
clean_sample <- tm_map(clean_sample, content_transformer(removeNumPunct))
```

Transform sample to all lower case

``` r
clean_sample <- tm_map(clean_sample, content_transformer(tolower))
```

Create profanity filter
Source: [List-of-Dirty-Naughty-Obscene-and-Otherwise-Bad-Words](%22List-of-Dirty-Naughty-Obscene-and-Otherwise-Bad-Words/en%22)

``` r
profanity <- read.table("./data/final/en_US/profanity.txt", header = FALSE, sep ="\n")
```

Remove profanity

``` r
clean_sample <- tm_map(clean_sample, removeWords, profanity[,1])
```

Remove stopwords

``` r
clean_sample <- tm_map(clean_sample, removeWords, stopwords("english"))
clean_sample <- tm_map(clean_sample, removeWords, stopwords("SMART"))
print(as.character(clean_sample[[1]]))
```

    ## [1] "love    onomatopoeia     made dumplings   scallions  cabbage    plenty  places      buy dumplings veggie       dollar    places     dont sell    white booth   seethrough window  sell delicious stuff   home  cook  inexpensive    good enjoy"

Remove Whitespace

``` r
clean_sample <- tm_map(clean_sample, stripWhitespace)
print(as.character(clean_sample[[1]]))
```

    ## [1] "love onomatopoeia made dumplings scallions cabbage plenty places buy dumplings veggie dollar places dont sell white booth seethrough window sell delicious stuff home cook inexpensive good enjoy"

Save clean corpus

``` r
saveRDS(clean_sample, file = "./data/final/en_US/clean_sample.rds" )
```

Initial Exploratory Data Analysis
---------------------------------

Convert to document term matrix

``` r
docterm_corpus <- DocumentTermMatrix(clean_sample)
dim(docterm_corpus)
```

    ## [1] 213483 141426

``` r
new_docterm_corpus <- removeSparseTerms(docterm_corpus,sparse = 0.993)
dim(new_docterm_corpus)
```

    ## [1] 213483    105

Find frequent terms

``` r
colS <- colSums(as.matrix(new_docterm_corpus))
length(colS)
```

    ## [1] 105

``` r
doc_features <- data.table(name = attributes(colS)$names, count = colS)
```

Most frequent and least frequent words

``` r
doc_features[order(-count)][1:10] #top 10 most frequent words
```

    ##       name count
    ##  1:   time 10625
    ##  2:   good  8913
    ##  3:   dont  8841
    ##  4:    day  8488
    ##  5:   love  8005
    ##  6: people  7808
    ##  7:   back  6948
    ##  8:   make  6422
    ##  9:  great  6089
    ## 10:   year  5835

``` r
doc_features[order(count)][1:10] #least 10 frequent words
```

    ##         name count
    ##  1:  started  1622
    ##  2:     stop  1627
    ##  3:  twitter  1661
    ##  4:   person  1662
    ##  5:     guys  1676
    ##  6: tomorrow  1677
    ##  7:      win  1694
    ##  8:     head  1695
    ##  9:     open  1703
    ## 10: business  1713

Plot most frequent terms

``` r
ggplot(doc_features[count>5000],aes(name, count)) +
  geom_bar(stat = "identity",fill='lightblue',color='black') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_economist() + scale_color_economist() 
```

![](01_Task_Script_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png)

Create word cloud

``` r
wordcloud(names(colS), colS, min.freq = 500, 
          colors = brewer.pal(6, 'Dark2'), random.order = FALSE)  
```

![](01_Task_Script_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-26-1.png)

``` r
wordcloud(names(colS), colS, min.freq = 2000, 
          colors = brewer.pal(6, 'Dark2'), random.order = FALSE)  
```

![](01_Task_Script_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-26-2.png)

``` r
end <- Sys.time()
(ellapsed <- end - start)
```

    ## Time difference of 2.203485 mins

------------------------------------------------------------------------

Session info
------------

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 15063)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2        ngram_3.0.3         wordcloud_2.5      
    ##  [4] RColorBrewer_1.1-2  ggthemes_3.4.0      data.table_1.10.4-2
    ##  [7] dtplyr_0.0.2        dplyr_0.7.4         purrr_0.2.3        
    ## [10] readr_1.1.1         tidyr_0.7.1         tibble_1.3.4       
    ## [13] ggplot2_2.2.1       tidyverse_1.1.1     knitr_1.17         
    ## [16] tm_0.7-1            NLP_0.1-11          downloader_0.4     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] slam_0.1-40      reshape2_1.4.2   haven_1.1.0      lattice_0.20-35 
    ##  [5] colorspace_1.3-2 htmltools_0.3.6  yaml_2.1.14      rlang_0.1.2     
    ##  [9] foreign_0.8-69   glue_1.1.1       modelr_0.1.1     readxl_1.0.0    
    ## [13] bindr_0.1        plyr_1.8.4       stringr_1.2.0    munsell_0.4.3   
    ## [17] gtable_0.2.0     cellranger_1.1.0 rvest_0.3.2      psych_1.7.8     
    ## [21] evaluate_0.10.1  labeling_0.3     forcats_0.2.0    parallel_3.4.2  
    ## [25] highr_0.6        broom_0.4.2      Rcpp_0.12.13     backports_1.1.1 
    ## [29] scales_0.5.0     jsonlite_1.5     mnormt_1.5-5     hms_0.3         
    ## [33] digest_0.6.12    stringi_1.1.5    grid_3.4.2       rprojroot_1.2   
    ## [37] tools_3.4.2      magrittr_1.5     lazyeval_0.2.0   pkgconfig_2.0.1 
    ## [41] xml2_1.1.1       lubridate_1.6.0  assertthat_0.2.0 rmarkdown_1.6   
    ## [45] httr_1.3.1       R6_2.2.2         nlme_3.1-131     compiler_3.4.2
