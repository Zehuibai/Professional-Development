Quiz 1
================
Mark Blackmore
2017-10-12

``` r
blogs_file   <- "./data/final/en_US/en_US.blogs.txt"
news_file    <- "./data/final/en_US/en_US.news.txt"
twitter_file <- "./data/final/en_US/en_US.twitter.txt"
```

1. The en\_US.blogs.txt file is how many megabytes?
---------------------------------------------------

``` r
file.size(blogs_file)/(2^20)
```

    ## [1] 200.4242

2. The en\_US.twitter.txt has how many lines of text?
-----------------------------------------------------

``` r
twitter <-  readLines(twitter_file, skipNul = TRUE) 
length(twitter)
```

    ## [1] 2360148

3. What is the length of the longest line seen in any of the three en\_US data sets?
------------------------------------------------------------------------------------

``` r
# Blogs
blog_lines <- nchar(readLines(blogs_file, skipNul = TRUE))
max(blog_lines)
```

    ## [1] 40835

``` r
# News
news_lines <- nchar(readLines(news_file, skipNul = TRUE))
max(news_lines)
```

    ## [1] 5760

``` r
# Twitter
twitter_lines <- nchar(readLines(twitter_file, skipNul = TRUE))
max(twitter_lines)
```

    ## [1] 213

4. In the en\_US twitter data set, if you divide the number of lines where the word "love" (all lowercase) occurs by the number of lines the word "hate" (all lowercase) occurs, about what do you get?
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

``` r
lines_love <- grepl(".love.", readLines(twitter_file, skipNul = TRUE), 
                    ignore.case = FALSE, perl = TRUE)

n_love <- sum(lines_love)

lines_hate <- grepl(".hate.", readLines(twitter_file, skipNul= TRUE), 
                    ignore.case = FALSE, perl = TRUE)

n_hate <- sum(lines_hate)

(love_hate <- n_love/n_hate)
```

    ## [1] 3.964942

5. The one tweet in the en\_US twitter data set that matches the word "biostats" says what?
-------------------------------------------------------------------------------------------

``` r
twitter[grep("biostats", twitter)]
```

    ## [1] "i know how you feel.. i have biostats on tuesday and i have yet to study =/"

6. How many tweets have the exact characters "A computer once beat me at chess, but it was no match for me at kickboxing". (I.e. the line matches those characters exactly.)
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

``` r
sum(grepl( "A computer once beat me at chess, but it was no match for me at kickboxing", twitter))
```

    ## [1] 3
