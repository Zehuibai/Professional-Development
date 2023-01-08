Task 05: Prediction Model
================
Mark Blackmore
2017-11-27

Setup

``` r
go <- Sys.time()
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})
```

Load Training Data

``` r
bi_words <- readRDS("./clean_repos/bi_words_fast.rds")
tri_words  <- readRDS("./clean_repos/tri_words_fast.rds")
quad_words <- readRDS("./clean_repos/quad_words_fast.rds")
```

Create Ngram Matching Functions

``` r
bigram <- function(input_words){
                    num <- length(input_words)
                    filter(bi_words, 
                          word1==input_words[num]) %>% 
                    top_n(1, n) %>%
                    filter(row_number() == 1L) %>%
                    select(num_range("word", 2)) %>%
                    as.character() -> out
                    ifelse(out =="character(0)", "?", return(out))
}

trigram <- function(input_words){
                    num <- length(input_words)
                    filter(tri_words, 
                            word1==input_words[num-1], 
                            word2==input_words[num])  %>% 
                    top_n(1, n) %>%
                    filter(row_number() == 1L) %>%
                    select(num_range("word", 3)) %>%
                    as.character() -> out
                    ifelse(out=="character(0)", bigram(input_words), return(out))
}

quadgram <- function(input_words){
                    num <- length(input_words)
                    filter(quad_words, 
                            word1==input_words[num-2], 
                            word2==input_words[num-1], 
                            word3==input_words[num])  %>% 
                    top_n(1, n) %>%
                    filter(row_number() == 1L) %>%
                    select(num_range("word", 4)) %>%
                    as.character() -> out
                    ifelse(out=="character(0)", trigram(input_words), return(out))
}
```

Create User Input and Data Cleaning Function; Calls the matching functions

``` r
ngrams <- function(input){
  # Create a dataframe
  input <- data_frame(text = input)
  # Clean the Inpput
  replace_reg <- "[^[:alpha:][:space:]]*"
  input <- input %>%
    mutate(text = str_replace_all(text, replace_reg, ""))
  # Find word count, separate words, lower case
  input_count <- str_count(input, boundary("word"))
  input_words <- unlist(str_split(input, boundary("word")))
  input_words <- tolower(input_words)
  # Call the matching functions
  out <- ifelse(input_count == 1, bigram(input_words), 
              ifelse (input_count == 2, trigram(input_words), quadgram(input_words)))
  # Output
  return(out)
}
```

User Input and Program Ouput

``` r
input <- "In case of a"
ngrams(input)
```

    ## [1] "new"

Time it

``` r
stop <- Sys.time()
(how_long <- stop - go)
```

    ## Time difference of 3.601568 secs

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
    ## [1] bindrcpp_0.2    stringr_1.2.0   dplyr_0.7.4     purrr_0.2.3    
    ## [5] readr_1.1.1     tidyr_0.7.1     tibble_1.3.4    ggplot2_2.2.1  
    ## [9] tidyverse_1.1.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.13     cellranger_1.1.0 compiler_3.4.2   plyr_1.8.4      
    ##  [5] bindr_0.1        forcats_0.2.0    tools_3.4.2      digest_0.6.12   
    ##  [9] lubridate_1.6.0  jsonlite_1.5     evaluate_0.10.1  nlme_3.1-131    
    ## [13] gtable_0.2.0     lattice_0.20-35  pkgconfig_2.0.1  rlang_0.1.2     
    ## [17] psych_1.7.8      yaml_2.1.14      parallel_3.4.2   haven_1.1.0     
    ## [21] xml2_1.1.1       httr_1.3.1       knitr_1.17       hms_0.3         
    ## [25] rprojroot_1.2    grid_3.4.2       glue_1.1.1       R6_2.2.2        
    ## [29] readxl_1.0.0     foreign_0.8-69   rmarkdown_1.6    modelr_0.1.1    
    ## [33] reshape2_1.4.2   magrittr_1.5     backports_1.1.1  scales_0.5.0    
    ## [37] htmltools_0.3.6  rvest_0.3.2      assertthat_0.2.0 mnormt_1.5-5    
    ## [41] colorspace_1.3-2 stringi_1.1.5    lazyeval_0.2.0   munsell_0.4.3   
    ## [45] broom_0.4.2
