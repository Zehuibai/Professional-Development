Task 04A: Fast Ngram Files
================
Mark Blackmore
2017-12-07

-   [Load the Data](#load-the-data)
-   [Sample the data](#sample-the-data)
-   [Clean the sample data](#clean-the-sample-data)
-   [Create all n-grams](#create-all-n-grams)
-   [Reduce n-grams files](#reduce-n-grams-files)
-   [What does the distribution on ngrams look like?](#what-does-the-distribution-on-ngrams-look-like)
-   [Separate words](#separate-words)
-   [Session info](#session-info)

Load the Data
-------------

English Repository Files

``` r
blogs_file   <- "./data/final/en_US/en_US.blogs.txt"
news_file    <- "./data/final/en_US/en_US.news.txt"
twitter_file <- "./data/final/en_US/en_US.twitter.txt"  
```

Read the data files

``` r
blogs   <- readLines(blogs_file, skipNul = TRUE)
news    <- readLines(news_file,  skipNul = TRUE)
twitter <- readLines(twitter_file, skipNul = TRUE)
```

Create dataframes

``` r
blogs   <- data_frame(text = blogs)
news    <- data_frame(text = news)
twitter <- data_frame(text = twitter)
```

Sample the data
---------------

``` r
set.seed(1001)
sample_pct <- 0.10

blogs_sample <- blogs %>%
  sample_n(., nrow(blogs)*sample_pct)
news_sample <- news %>%
  sample_n(., nrow(news)*sample_pct)
twitter_sample <- twitter %>%
  sample_n(., nrow(twitter)*sample_pct)
```

Create tidy repository

``` r
repo_sample <- bind_rows(mutate(blogs_sample, source = "blogs"),
                         mutate(news_sample,  source = "news"),
                         mutate(twitter_sample, source = "twitter")) 
repo_sample$source <- as.factor(repo_sample$source)
```

Clean up

``` r
rm(list = c("blogs", "blogs_file", "blogs_sample","news", "news_file",     
            "news_sample", "sample_pct", "twitter","twitter_file", 
            "twitter_sample"))
```

Clean the sample data
---------------------

Create filters: non-alphanumeric's, url's, repeated letters(+3x)

``` r
replace_reg <- "[^[:alpha:][:space:]]*"
replace_url <- "http[^[:space:]]*"
replace_aaa <- "\\b(?=\\w*(\\w)\\1)\\w+\\b"  
```

Clean the sample. Cleaning is separted from tidying so `unnest_tokens` function can be used for words, and ngrams.

``` r
clean_sample <-  repo_sample %>%
  mutate(text = str_replace_all(text, replace_reg, "")) %>%
  mutate(text = str_replace_all(text, replace_url, "")) %>%
  mutate(text = str_replace_all(text, replace_aaa, "")) %>% 
  mutate(text = iconv(text, "ASCII//TRANSLIT"))

rm(list = c("repo_sample"))
```

Create all n-grams
------------------

Bigrams

``` r
bigram_repo <- clean_sample  %>%
  unnest_tokens(bigram, text, token = "ngrams", n = 2)
```

Trigrams

``` r
trigram_repo <- clean_sample  %>%
  unnest_tokens(trigram, text, token = "ngrams", n = 3)
```

Quadgrams

``` r
quadgram_repo <- clean_sample  %>%
  unnest_tokens(quadgram, text, token = "ngrams", n = 4)
```

Quintgrams

``` r
quintgram_repo <- clean_sample  %>%
  unnest_tokens(quintgram, text, token = "ngrams", n = 5)
```

Sextgrams

``` r
sextgram_repo <- clean_sample  %>%
  unnest_tokens(sextgram, text, token = "ngrams", n = 6)
```

Reduce n-grams files
--------------------

Bigrams

``` r
bigram_cover <- bigram_repo %>%
  count(bigram) %>%  
  filter(n > 10) %>%
  arrange(desc(n))  
rm(list = c("bigram_repo"))
```

Trigrams

``` r
trigram_cover <- trigram_repo %>%
  count(trigram) %>%  
  filter(n > 10) %>%
  arrange(desc(n))  
rm(list = c("trigram_repo"))
```

Quadgrams

``` r
quadgram_cover <- quadgram_repo %>%
  count(quadgram) %>%  
  filter(n > 10) %>%
  arrange(desc(n))  
rm(list = c("quadgram_repo"))
```

Quintgrams

``` r
quintgram_cover <- quintgram_repo %>%
  count(quintgram) %>%  
  filter(n > 10) %>%
  arrange(desc(n))  
rm(list = c("quintgram_repo"))
```

Sextgrams

``` r
sextgram_cover <- sextgram_repo %>%
  count(sextgram) %>%  
  filter(n > 10) %>%
  arrange(desc(n))  
rm(list = c("sextgram_repo"))
```

What does the distribution on ngrams look like?
-----------------------------------------------

``` r
disty <- data_frame(ngram = c(rep("bigrams",  nrow(bigram_cover)),
                             rep("trigrams",  nrow(trigram_cover)),
                             rep("quadgrams", nrow(quadgram_cover)),
                             rep("quintgrams", nrow(quintgram_cover)),
                             rep("sextgrams",  nrow(sextgram_cover))),
                    number = c(bigram_cover$n,  trigram_cover$n, 
                              quadgram_cover$n, quintgram_cover$n,
                              sextgram_cover$n))
disty
```

    ## # A tibble: 150,271 x 2
    ##      ngram number
    ##      <chr>  <int>
    ##  1 bigrams  44296
    ##  2 bigrams  41480
    ##  3 bigrams  25015
    ##  4 bigrams  20646
    ##  5 bigrams  19953
    ##  6 bigrams  16338
    ##  7 bigrams  16290
    ##  8 bigrams  14462
    ##  9 bigrams  14357
    ## 10 bigrams  12060
    ## # ... with 150,261 more rows

``` r
disty$ngram <- as.factor(disty$ngram)
ggplot(data = disty, aes(y = number, x = reorder(ngram, -number))) +
  geom_boxplot() + scale_y_log10() +
  xlab("ngram")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-1.png)

``` r
ggsave("./ngram_match/www/ngrams.png")
```

    ## Saving 7 x 5 in image

``` r
sextgram_cover %>%
  top_n(15, n) %>%
  mutate(sextgram = reorder(sextgram, n)) %>%
  ggplot(aes(sextgram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip() +
  ggtitle("Sextgrams")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-2.png)

``` r
ggsave("./ngram_match/www/sextgrams.png")
```

    ## Saving 7 x 5 in image

``` r
quintgram_cover %>%
  top_n(15, n) %>%
  mutate(quintgram = reorder(quintgram, n)) %>%
  ggplot(aes(quintgram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip() +
  ggtitle("Quintgrams")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-3.png)

``` r
ggsave("./ngram_match/www/quintgrams.png")
```

    ## Saving 7 x 5 in image

``` r
quadgram_cover %>%
  top_n(15, n) %>%
  mutate(quadgram = reorder(quadgram, n)) %>%
  ggplot(aes(quadgram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip() +
  ggtitle("Quadgrams")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-4.png)

``` r
ggsave("./ngram_match/www/quadgrams.png")
```

    ## Saving 7 x 5 in image

``` r
trigram_cover %>%
  top_n(15, n) %>%
  mutate(trigram = reorder(trigram, n)) %>%
  ggplot(aes(trigram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip() +
  ggtitle("Trigrams")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-5.png)

``` r
ggsave("./ngram_match/www/trigrams.png")
```

    ## Saving 7 x 5 in image

``` r
bigram_cover %>%
  top_n(15, n) %>%
  mutate(bigram = reorder(bigram, n)) %>%
  ggplot(aes(bigram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip() +
  ggtitle("Bigrams")
```

![](04A_Task_Script_files/figure-markdown_github-ascii_identifiers/DistyPlot-6.png)

``` r
ggsave("./ngram_match/www/bigrams.png")
```

    ## Saving 7 x 5 in image

Separate words
--------------

``` r
bi_words <- bigram_cover %>%
  separate(bigram, c("word1", "word2"), sep = " ")
bi_words
```

    ## # A tibble: 87,239 x 3
    ##    word1 word2     n
    ##  * <chr> <chr> <int>
    ##  1    of   the 44296
    ##  2    in   the 41480
    ##  3    to   the 25015
    ##  4   for   the 20646
    ##  5    on   the 19953
    ##  6     c     e 16338
    ##  7    to    be 16290
    ##  8    at   the 14462
    ##  9   and   the 14357
    ## 10    in     a 12060
    ## # ... with 87,229 more rows

``` r
tri_words <- trigram_cover %>%
  separate(trigram, c("word1", "word2", "word3"), sep = " ")
tri_words
```

    ## # A tibble: 51,539 x 4
    ##     word1 word2 word3     n
    ##  *  <chr> <chr> <chr> <int>
    ##  1    one    of   the  3543
    ##  2      a   lot    of  2929
    ##  3 thanks   for   the  2501
    ##  4     to    be     a  1888
    ##  5  going    to    be  1825
    ##  6    the    of   the  1814
    ##  7      i  want    to  1629
    ##  8    the   end    of  1522
    ##  9    out    of   the  1520
    ## 10     it   was     a  1511
    ## # ... with 51,529 more rows

``` r
quad_words <- quadgram_cover %>%
  separate(quadgram, c("word1", "word2", "word3", "word4"), sep = " ")
quad_words
```

    ## # A tibble: 9,941 x 5
    ##    word1 word2 word3 word4     n
    ##  * <chr> <chr> <chr> <chr> <int>
    ##  1   the   end    of   the   806
    ##  2    at   the   end    of   656
    ##  3   the  rest    of   the   651
    ##  4   for   the first  time   613
    ##  5    at   the  same  time   506
    ##  6    is going    to    be   482
    ##  7    is   one    of   the   446
    ##  8   one    of   the  most   435
    ##  9  when    it comes    to   419
    ## 10 going    to    be     a   376
    ## # ... with 9,931 more rows

``` r
quint_words <- quintgram_cover %>%
  separate(quintgram, c("word1", "word2", "word3", "word4", "word5"), sep = " ")
quint_words
```

    ## # A tibble: 1,317 x 6
    ##    word1 word2 word3 word4 word5     n
    ##  * <chr> <chr> <chr> <chr> <chr> <int>
    ##  1    at   the   end    of   the   371
    ##  2   for   the first  time    in   145
    ##  3   the   end    of   the   day   128
    ##  4   for   the  rest    of   the   124
    ##  5    by   the   end    of   the   122
    ##  6 thank   you    so  much   for   108
    ##  7 there   are     a   lot    of   108
    ##  8   its going    to    be     a    98
    ##  9   let    me  know    if   you    94
    ## 10    is going    to    be     a    93
    ## # ... with 1,307 more rows

``` r
sext_words <- sextgram_cover %>%
  separate(sextgram, c("word1", "word2", "word3", "word4", "word5", "word6"), sep = " ")
sext_words
```

    ## # A tibble: 235 x 7
    ##    word1    word2    word3 word4 word5 word6     n
    ##  * <chr>    <chr>    <chr> <chr> <chr> <chr> <int>
    ##  1    at      the      end    of   the   day   107
    ##  2     c        c        c     c     c     c    63
    ##  3    on      the    other  side    of   the    57
    ##  4    me       me       me    me    me    me    40
    ##  5  this       is    going    to    be     a    36
    ##  6  cake     cake     cake  cake  cake  cake    34
    ##  7     i     just finished     a    mi   run    34
    ##  8  just finished        a    mi   run  with    33
    ##  9 thank      you       so  much   for   the    28
    ## 10   let       me     know  what   you think    27
    ## # ... with 225 more rows

Save data for the Shiny App

``` r
saveRDS(bi_words, "./ngram_match/app_data/bi_words_fast.rds")
saveRDS(tri_words, "./ngram_match/app_data/tri_words_fast.rds")
saveRDS(quad_words,"./ngram_match/app_data/quad_words_fast.rds")
saveRDS(quint_words,"./ngram_match/app_data/quint_words_fast.rds")
saveRDS(sext_words,"./ngram_match/app_data/sext_words_fast.rds")
```

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
    ##  [1] bindrcpp_0.2    knitr_1.17      stringr_1.2.0   dplyr_0.7.4    
    ##  [5] purrr_0.2.3     readr_1.1.1     tidyr_0.7.1     tibble_1.3.4   
    ##  [9] ggplot2_2.2.1   tidyverse_1.1.1 tidytext_0.1.4 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.13      cellranger_1.1.0  compiler_3.4.2   
    ##  [4] plyr_1.8.4        bindr_0.1         forcats_0.2.0    
    ##  [7] tokenizers_0.1.4  tools_3.4.2       digest_0.6.12    
    ## [10] lubridate_1.6.0   jsonlite_1.5      gtable_0.2.0     
    ## [13] evaluate_0.10.1   nlme_3.1-131      lattice_0.20-35  
    ## [16] pkgconfig_2.0.1   rlang_0.1.2       Matrix_1.2-11    
    ## [19] psych_1.7.8       yaml_2.1.14       parallel_3.4.2   
    ## [22] haven_1.1.0       xml2_1.1.1        httr_1.3.1       
    ## [25] janeaustenr_0.1.5 hms_0.3           tidyselect_0.2.2 
    ## [28] rprojroot_1.2     grid_3.4.2        glue_1.1.1       
    ## [31] R6_2.2.2          readxl_1.0.0      foreign_0.8-69   
    ## [34] rmarkdown_1.6     modelr_0.1.1      reshape2_1.4.2   
    ## [37] magrittr_1.5      scales_0.5.0      backports_1.1.1  
    ## [40] SnowballC_0.5.1   htmltools_0.3.6   rvest_0.3.2      
    ## [43] assertthat_0.2.0  mnormt_1.5-5      colorspace_1.3-2 
    ## [46] labeling_0.3      stringi_1.1.5     lazyeval_0.2.0   
    ## [49] munsell_0.4.3     broom_0.4.2
