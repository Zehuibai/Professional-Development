Exploring the `tm` Package
================
Mark Blackmore
2017-11-28

-   [Introduction](#introduction)
-   [Example with `ovid`](#example-with-ovid)
-   [Example with `Reuters-21578`](#example-with-reuters-21578)
-   [Session info](#session-info)

Introduction
------------

This document works through several examples in the article:
[Text mining infrastucture in R](http://www.jstatsoft.org/v25/i05/)

This is one of the references cited in the course document: *Task 0: Understanding the Problem*

``` r
library(tm)
```

    ## Loading required package: NLP

Example with `ovid`
-------------------

### Find the sample corpus

``` r
txt <- system.file("texts", "txt", package = "tm")
```

### Load the corpus

``` r
(ovid <- Corpus(DirSource(txt, encoding = "UTF-8"),
                readerControl = list(reader = readPlain,
                                     language = "lat",
                                     load = TRUE)))
```

    ## <<SimpleCorpus>>
    ## Metadata:  corpus specific: 1, document level (indexed): 0
    ## Content:  documents: 5

### Examine metadata

Show first document and update the author's name

``` r
meta(ovid[[1]])
```

    ##   author       : character(0)
    ##   datetimestamp: 2017-11-29 00:50:06
    ##   description  : character(0)
    ##   heading      : character(0)
    ##   id           : ovid_1.txt
    ##   language     : lat
    ##   origin       : character(0)

``` r
meta(ovid[[1]])$author <- "Publius Ovidius Naso"
meta(ovid[[1]])
```

    ##   author       : character(0)
    ##   datetimestamp: 2017-11-29 00:50:06
    ##   description  : character(0)
    ##   heading      : character(0)
    ##   id           : content
    ##   language     : lat
    ##   origin       : character(0)

``` r
# same result
ovid[[1]][2]
```

    ## $meta
    ##   author       : character(0)
    ##   datetimestamp: 2017-11-29 00:50:06
    ##   description  : character(0)
    ##   heading      : character(0)
    ##   id           : content
    ##   language     : lat
    ##   origin       : character(0)

### Examine the first document's text

``` r
ovid[[1]][1]
```

    ## $content
    ## [1] "    Si quis in hoc artem populo non novit amandi,\n         hoc legat et lecto carmine doctus amet.\n    arte citae veloque rates remoque moventur,\n         arte leves currus: arte regendus amor.\n\n    curribus Automedon lentisque erat aptus habenis,\n         Tiphys in Haemonia puppe magister erat:\n    me Venus artificem tenero praefecit Amori;\n         Tiphys et Automedon dicar Amoris ego.\n    ille quidem ferus est et qui mihi saepe repugnet:\n\n         sed puer est, aetas mollis et apta regi.\n    Phillyrides puerum cithara perfecit Achillem,\n         atque animos placida contudit arte feros.\n    qui totiens socios, totiens exterruit hostes,\n         creditur annosum pertimuisse senem."

### Concatenate several text collections to a single one

``` r
c(ovid[1:2], ovid[3:4])
```

    ## $content
    ## $content$content
    ## [1] "    Si quis in hoc artem populo non novit amandi,\n         hoc legat et lecto carmine doctus amet.\n    arte citae veloque rates remoque moventur,\n         arte leves currus: arte regendus amor.\n\n    curribus Automedon lentisque erat aptus habenis,\n         Tiphys in Haemonia puppe magister erat:\n    me Venus artificem tenero praefecit Amori;\n         Tiphys et Automedon dicar Amoris ego.\n    ille quidem ferus est et qui mihi saepe repugnet:\n\n         sed puer est, aetas mollis et apta regi.\n    Phillyrides puerum cithara perfecit Achillem,\n         atque animos placida contudit arte feros.\n    qui totiens socios, totiens exterruit hostes,\n         creditur annosum pertimuisse senem."
    ## 
    ## $content$meta
    ##   author       : Publius Ovidius Naso
    ##   datetimestamp: 2017-11-29 00:50:06
    ##   description  : character(0)
    ##   heading      : character(0)
    ##   id           : ovid_1.txt
    ##   language     : lat
    ##   origin       : character(0)
    ## 
    ## 
    ## $meta
    ## $language
    ## [1] "lat"
    ## 
    ## attr(,"class")
    ## [1] "CorpusMeta"
    ## 
    ## $dmeta
    ## data frame with 0 columns and 2 rows
    ## 
    ## $content
    ## $content$<NA>
    ## NULL
    ## 
    ## $content$<NA>
    ## NULL
    ## 
    ## 
    ## $meta
    ## $language
    ## [1] "lat"
    ## 
    ## attr(,"class")
    ## [1] "CorpusMeta"
    ## 
    ## $dmeta
    ## data frame with 0 columns and 2 rows

### Show the number of documents in the corpus

``` r
length(ovid)
```

    ## [1] 2

### Show detailed summary of the text document collection

``` r
summary(ovid)
```

    ##         Length Class             Mode
    ## content 2      PlainTextDocument list
    ## meta    2      PlainTextDocument list

### Show predefined transformations

Can be applied to a corupus with `tm_map`

``` r
getTransformations()
```

    ## [1] "removeNumbers"     "removePunctuation" "removeWords"      
    ## [4] "stemDocument"      "stripWhitespace"

#### Remove punctuation

``` r
ovid <- tm_map(ovid, FUN = removePunctuation)
ovid[[1]][1]
```

    ## $content
    ## [1] "    Si quis in hoc artem populo non novit amandi\n         hoc legat et lecto carmine doctus amet\n    arte citae veloque rates remoque moventur\n         arte leves currus arte regendus amor\n\n    curribus Automedon lentisque erat aptus habenis\n         Tiphys in Haemonia puppe magister erat\n    me Venus artificem tenero praefecit Amori\n         Tiphys et Automedon dicar Amoris ego\n    ille quidem ferus est et qui mihi saepe repugnet\n\n         sed puer est aetas mollis et apta regi\n    Phillyrides puerum cithara perfecit Achillem\n         atque animos placida contudit arte feros\n    qui totiens socios totiens exterruit hostes\n         creditur annosum pertimuisse senem"

#### Remove numbers

``` r
ovid <- tm_map(ovid, FUN = removeNumbers)
ovid[[1]][1]
```

    ## $content
    ## [1] "    Si quis in hoc artem populo non novit amandi\n         hoc legat et lecto carmine doctus amet\n    arte citae veloque rates remoque moventur\n         arte leves currus arte regendus amor\n\n    curribus Automedon lentisque erat aptus habenis\n         Tiphys in Haemonia puppe magister erat\n    me Venus artificem tenero praefecit Amori\n         Tiphys et Automedon dicar Amoris ego\n    ille quidem ferus est et qui mihi saepe repugnet\n\n         sed puer est aetas mollis et apta regi\n    Phillyrides puerum cithara perfecit Achillem\n         atque animos placida contudit arte feros\n    qui totiens socios totiens exterruit hostes\n         creditur annosum pertimuisse senem"

#### Change to all lower case

``` r
ovid <- tm_map(ovid, FUN = content_transformer(tolower))
ovid[[1]][1]
```

    ## $content
    ## [1] "    si quis in hoc artem populo non novit amandi\n         hoc legat et lecto carmine doctus amet\n    arte citae veloque rates remoque moventur\n         arte leves currus arte regendus amor\n\n    curribus automedon lentisque erat aptus habenis\n         tiphys in haemonia puppe magister erat\n    me venus artificem tenero praefecit amori\n         tiphys et automedon dicar amoris ego\n    ille quidem ferus est et qui mihi saepe repugnet\n\n         sed puer est aetas mollis et apta regi\n    phillyrides puerum cithara perfecit achillem\n         atque animos placida contudit arte feros\n    qui totiens socios totiens exterruit hostes\n         creditur annosum pertimuisse senem"

#### Stem the corpus

``` r
# ovid <- tm_map(ovid, FUN = stemDocument)
# function is not available for language 'la'
```

#### Remove words

``` r
axe_words <- c("mater", "seu", "annis", "")
ovid <- tm_map(ovid, FUN = removeWords, axe_words)
ovid[[1]][1]
```

    ## $content
    ## [1] "    si quis in hoc artem populo non novit amandi\n         hoc legat et lecto carmine doctus amet\n    arte citae veloque rates remoque moventur\n         arte leves currus arte regendus amor\n\n    curribus automedon lentisque erat aptus habenis\n         tiphys in haemonia puppe magister erat\n    me venus artificem tenero praefecit amori\n         tiphys et automedon dicar amoris ego\n    ille quidem ferus est et qui mihi saepe repugnet\n\n         sed puer est aetas mollis et apta regi\n    phillyrides puerum cithara perfecit achillem\n         atque animos placida contudit arte feros\n    qui totiens socios totiens exterruit hostes\n         creditur annosum pertimuisse senem"

#### Remove whitespace

``` r
ovid <- tm_map(ovid, FUN = stripWhitespace)
ovid[[1]][1]
```

    ## $content
    ## [1] " si quis in hoc artem populo non novit amandi hoc legat et lecto carmine doctus amet arte citae veloque rates remoque moventur arte leves currus arte regendus amor curribus automedon lentisque erat aptus habenis tiphys in haemonia puppe magister erat me venus artificem tenero praefecit amori tiphys et automedon dicar amoris ego ille quidem ferus est et qui mihi saepe repugnet sed puer est aetas mollis et apta regi phillyrides puerum cithara perfecit achillem atque animos placida contudit arte feros qui totiens socios totiens exterruit hostes creditur annosum pertimuisse senem"

Example with `Reuters-21578`
----------------------------

``` r
reut21578 <- system.file("texts", "crude", package = "tm")
reuters <- VCorpus(DirSource(reut21578),
                   readerControl = list(reader = readReut21578XMLasPlain))
reuters[[1]][1]
```

    ## $content
    ## [1] "Diamond Shamrock Corp said that\neffective today it had cut its contract prices for crude oil by\n1.50 dlrs a barrel.\n    The reduction brings its posted price for West Texas\nIntermediate to 16.00 dlrs a barrel, the copany said.\n    \"The price reduction today was made in the light of falling\noil product prices and a weak crude oil market,\" a company\nspokeswoman said.\n    Diamond is the latest in a line of U.S. oil companies that\nhave cut its contract, or posted, prices over the last two days\nciting weak oil markets.\n Reuter"

### Eliminate Extra Whitespace

``` r
reuters <- tm_map(reuters, stripWhitespace)
```

### Convert to Lower Case

``` r
reuters <- tm_map(reuters, content_transformer(tolower))
```

### Stemming

``` r
tm_map(reuters, stemDocument)
```

    ## <<VCorpus>>
    ## Metadata:  corpus specific: 0, document level (indexed): 0
    ## Content:  documents: 20

### Remove Stop Words

``` r
reuters <- tm_map(reuters, removeWords, stopwords("english"))
reuters[[1]][1]
```

    ## $content
    ## [1] "diamond shamrock corp said  effective today   cut  contract prices  crude oil  1.50 dlrs  barrel.  reduction brings  posted price  west texas intermediate  16.00 dlrs  barrel,  copany said. \" price reduction today  made   light  falling oil product prices   weak crude oil market,\"  company spokeswoman said. diamond   latest   line  u.s. oil companies   cut  contract,  posted, prices   last two days citing weak oil markets. reuter"

### Creating Document-Term Matrices

``` r
dtm <- DocumentTermMatrix(reuters)
inspect(dtm[5:10, 740:743])
```

    ## <<DocumentTermMatrix (documents: 6, terms: 4)>>
    ## Non-/sparse entries: 8/16
    ## Sparsity           : 67%
    ## Maximal term length: 6
    ## Weighting          : term frequency (tf)
    ## Sample             :
    ##      Terms
    ## Docs  one, opec opec's opec"s
    ##   211    0    0      0      0
    ##   236    0    6      0      2
    ##   237    0    1      0      0
    ##   242    1    2      0      0
    ##   246    0    1      1      0
    ##   248    0    6      0      0

### Operations on Document-Term Matrices

``` r
findFreqTerms(dtm, 5)
```

    ##  [1] "15.8"          "abdul-aziz"    "ability"       "accord"       
    ##  [5] "agency"        "agreement"     "ali"           "also"         
    ##  [9] "analysts"      "arab"          "arabia"        "barrel."      
    ## [13] "barrels"       "billion"       "bpd"           "budget"       
    ## [17] "company"       "crude"         "daily"         "demand"       
    ## [21] "dlrs"          "economic"      "emergency"     "energy"       
    ## [25] "exchange"      "expected"      "exports"       "futures"      
    ## [29] "government"    "group"         "gulf"          "help"         
    ## [33] "hold"          "industry"      "international" "january"      
    ## [37] "kuwait"        "last"          "market"        "may"          
    ## [41] "meeting"       "minister"      "mln"           "month"        
    ## [45] "nazer"         "new"           "now"           "nymex"        
    ## [49] "official"      "oil"           "one"           "opec"         
    ## [53] "output"        "pct"           "petroleum"     "plans"        
    ## [57] "posted"        "present"       "price"         "prices"       
    ## [61] "prices,"       "prices."       "production"    "quota"        
    ## [65] "quoted"        "recent"        "report"        "research"     
    ## [69] "reserve"       "reuter"        "said"          "said."        
    ## [73] "saudi"         "sell"          "sheikh"        "sources"      
    ## [77] "study"         "traders"       "u.s."          "united"       
    ## [81] "west"          "will"          "world"

### Find Associations Between Words

Find words associated with opec, with at least 0.8 correlation

``` r
findAssocs(dtm, "opec", 0.8)
```

    ## $opec
    ##   meeting emergency       oil      15.8  analysts    buyers      said 
    ##      0.88      0.87      0.87      0.85      0.85      0.83      0.82 
    ##   ability 
    ##      0.80

### Remove Sparse Terms

``` r
inspect(removeSparseTerms(dtm, 0.4))
```

    ## <<DocumentTermMatrix (documents: 20, terms: 3)>>
    ## Non-/sparse entries: 58/2
    ## Sparsity           : 3%
    ## Maximal term length: 6
    ## Weighting          : term frequency (tf)
    ## Sample             :
    ##      Terms
    ## Docs  oil reuter said
    ##   127   5      1    1
    ##   144  11      1    9
    ##   236   7      1    6
    ##   242   3      1    3
    ##   246   4      1    4
    ##   248   9      1    5
    ##   273   5      1    5
    ##   352   5      1    1
    ##   489   4      1    2
    ##   502   4      1    2

### Dictionary: Terms to Text Mine

``` r
inspect(DocumentTermMatrix(reuters,
                           list(dictionary = c("prices", "crude", "oil"))))
```

    ## <<DocumentTermMatrix (documents: 20, terms: 3)>>
    ## Non-/sparse entries: 41/19
    ## Sparsity           : 32%
    ## Maximal term length: 6
    ## Weighting          : term frequency (tf)
    ## Sample             :
    ##      Terms
    ## Docs  crude oil prices
    ##   127     1   1      1
    ##   144     0   1      1
    ##   191     1   1      0
    ##   194     1   1      0
    ##   236     1   1      1
    ##   242     0   1      1
    ##   248     0   1      1
    ##   273     1   1      1
    ##   353     1   1      1
    ##   543     1   1      1

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
    ## [1] tm_0.7-1   NLP_0.1-11
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.13    XML_3.98-1.9    digest_0.6.12   rprojroot_1.2  
    ##  [5] SnowballC_0.5.1 slam_0.1-40     backports_1.1.1 magrittr_1.5   
    ##  [9] evaluate_0.10.1 stringi_1.1.5   rmarkdown_1.6   tools_3.4.2    
    ## [13] stringr_1.2.0   yaml_2.1.14     parallel_3.4.2  compiler_3.4.2 
    ## [17] htmltools_0.3.6 knitr_1.17
