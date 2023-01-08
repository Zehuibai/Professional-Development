#' ---
#' title: 'Task 04: Working toward a Prediction Model'
#' author: "Mark Blackmore"
#' date: "`r format(Sys.Date())`"
#' output: 
#'   github_document:
#'     toc: true
#' ---

#' This script creates the ngram files used to predict ngrams based on user input.
#' These files are used by prediction functions found in `05_Task_Script.R`

#+ setup, echo=FALSE
suppressPackageStartupMessages({
  library(tidytext)
  library(tidyverse)
  library(stringr)
  library(knitr)
  library(wordcloud)
  library(ngram)
})

#' ## Load the Data
#+ DataLoading

#' English Repository Files
blogs_file   <- "./data/final/en_US/en_US.blogs.txt"
news_file    <- "./data/final/en_US/en_US.news.txt"
twitter_file <- "./data/final/en_US/en_US.twitter.txt"  

#' Read the data files
blogs   <- readLines(blogs_file, skipNul = TRUE)
news    <- readLines(news_file,  skipNul = TRUE)
twitter <- readLines(twitter_file, skipNul = TRUE)

#' Read the data files into dataframes
blogs   <- data_frame(text = blogs)
news    <- data_frame(text = news)
twitter <- data_frame(text = twitter)

#' ## Sample the data
#+ DataSampling
set.seed(1001)
sample_pct <- 0.05

blogs_sample <- blogs %>%
  sample_n(., nrow(blogs)*sample_pct)
news_sample <- news %>%
  sample_n(., nrow(news)*sample_pct)
twitter_sample <- twitter %>%
  sample_n(., nrow(twitter)*sample_pct)

#' Create tidy repository
repo_sample <- bind_rows(mutate(blogs_sample, source = "blogs"),
                         mutate(news_sample,  source = "news"),
                         mutate(twitter_sample, source = "twitter")) 
repo_sample$source <- as.factor(repo_sample$source)

#' ## Clean the data
#' Create filters: stopwords, profanity, non-alphanumeric's, url's, repeated letters(+3x)
#+ DataCleaning
data("stop_words")
swear_words <- read_delim("./data/final/en_US/en_US.swearWords.csv", delim = "\n", col_names = FALSE)
swear_words <- unnest_tokens(swear_words, word, X1)
replace_reg <- "[^[:alpha:][:space:]]*"
replace_url <- "http[^[:space:]]*"
replace_aaa <- "\\b(?=\\w*(\\w)\\1)\\w+\\b"  

#' Clean the sample. Cleaning is separted from tidying so `unnest_tokens` function can be used for words,
#' and ngrams.
clean_sample <-  repo_sample %>%
  mutate(text = str_replace_all(text, replace_reg, "")) %>%
  mutate(text = str_replace_all(text, replace_url, "")) %>%
  mutate(text = str_replace_all(text, replace_aaa, "")) %>% 
  mutate(text = iconv(text, "ASCII//TRANSLIT"))

#' ## Create all n-grams  
#' Unigrams  
tidy_repo <- clean_sample %>%
  unnest_tokens(word, text) %>%
  anti_join(swear_words) %>%
  anti_join(stop_words)

#' Bigrams  
bigram_repo <- clean_sample  %>%
  unnest_tokens(bigram, text, token = "ngrams", n = 2)

#' Trigrams  
trigram_repo <- clean_sample  %>%
  unnest_tokens(trigram, text, token = "ngrams", n = 3)

#' Quadgrams  
quadgram_repo <- clean_sample  %>%
  unnest_tokens(quadgram, text, token = "ngrams", n = 4)

#' ## Reduce n-grams to top 50% of CDF  
#' Unigram upper half
cover_50 <- tidy_repo %>%
  count(word) %>%  
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%  
  mutate(coverage = cumsum(proportion)) %>%
  filter(coverage <= 0.5)

#' Bigram upper half
bigram_cover_50 <- bigram_repo %>%
  count(bigram) %>%  
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%  
  mutate(coverage = cumsum(proportion)) %>%
  filter(coverage <= 0.5)

#' Trigram upper half
trigram_cover_50 <- trigram_repo %>%
  count(trigram) %>%  
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%  
  mutate(coverage = cumsum(proportion)) %>%
  filter(coverage <= 0.5)

#' Quadgram upper half
quadgram_cover_50 <- quadgram_repo %>%
  count(quadgram) %>%  
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%  
  mutate(coverage = cumsum(proportion)) %>%
  filter(coverage <= 0.5)

#' ## Separate words
bi_words <- bigram_cover_50 %>%
  separate(bigram, c("word1", "word2"), sep = " ")
bi_words

tri_words <- trigram_cover_50 %>%
  separate(trigram, c("word1", "word2", "word3"), sep = " ")
tri_words

quad_words <- quadgram_cover_50 %>%
  separate(quadgram, c("word1", "word2", "word3", "word4"), sep = " ")
quad_words

#' Save separated words for prediction
saveRDS(bi_words, "./clean_repos/bi_words.rds")
saveRDS(tri_words, "./clean_repos/tri_words.rds")
saveRDS(quad_words, "./clean_repos/quad_words.rds")

#' ## Clear workspace, time load
# rm(list= ls())

go <- Sys.time()
library(tidyverse)
library(stringr)
bi_words <- readRDS("./clean_repos/bi_words.rds")
tri_words  <- readRDS("./clean_repos/tri_words.rds")
quad_words <- readRDS("./clean_repos/quad_words.rds")

stop <- Sys.time()
(how_long <- stop - go)

#' ## What does the distribution of ngrams look like?  
#' Suggests there may be a better way to subset. See `04A_Task_Script.R`
disty = data_frame(ngram = c(rep("bigrams",   nrow(bigram_cover_50)),
                             rep("trigrams",  nrow(trigram_cover_50)),
                             rep("quadgrams", nrow(quadgram_cover_50))), 
                   number = c(bigram_cover_50$n, trigram_cover_50$n, quadgram_cover_50$n))

disty$ngram <- as.factor(disty$ngram)
ggplot(data = disty, aes(y = number, x = ngram)) + geom_boxplot() + scale_y_log10()

#' -------------
#'  
#' ## Session info
sessionInfo()