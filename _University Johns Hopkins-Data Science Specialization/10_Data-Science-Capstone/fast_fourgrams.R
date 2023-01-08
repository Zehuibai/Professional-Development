
start <- Sys.time()

#+ startup, echo = FALSE 
suppressPackageStartupMessages({
  library(tidytext)
  library(tidyverse)
  library(stringr)
  library(knitr)
  library(wordcloud)
})

#' Read intermediate files for n-grams
clean_blogs <- readRDS("./data/final/en_US/clean_blogs.rds")
clean_news  <- readRDS("./data/final/en_US/clean_news.rds")
clean_twitter <- readRDS("./data/final/en_US/clean_twitter.rds")

#' Set sample percentage
pct <- 0.2

#' Create Fourgrams by source using `unnest_tokens`

blogs_fourgrams <- clean_blogs  %>%
  sample_n(., nrow(clean_blogs)*pct) %>%
  unnest_tokens(fourgram, text, token = "ngrams", n = 4)

news_fourgrams <- clean_news  %>%
  sample_n(., nrow(clean_news)*pct) %>%
  unnest_tokens(fourgram, text, token = "ngrams", n = 4)

twitter_fourgrams <- clean_twitter  %>%
  sample_n(., nrow(clean_twitter)*pct) %>%
  unnest_tokens(fourgram, text, token = "ngrams", n = 4)

#' Create tidy fourgram repository
fourgram_repo <- bind_rows(mutate(blogs_fourgrams, source = "blogs"),
                           mutate(news_fourgrams,  source = "news"),
                           mutate(twitter_fourgrams, source = "twitter"))
fourgram_repo$source <- as.factor(fourgram_repo$source)

#' Number of fourgrams to attain 90% coverage of all fourgrams in repo
fourgram_cover_90 <- fourgram_repo %>%
  count(fourgram) %>%  
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%  
  mutate(coverage = cumsum(proportion)) %>%
  filter(coverage <= 0.9)
nrow(fourgram_cover_90)

#' Fourgram distribution
fourgram_cover_90 %>%
  #count(trigram, sort = TRUE) %>%
  filter(n > 400) %>%
  mutate(fourgram = reorder(fourgram, n)) %>%
  ggplot(aes(fourgram, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip()

fourgrams_separated <- fourgram_cover_90 %>%
  separate(fourgram, c("word1", "word2", "word3", "word4"), sep = " ")
fourgrams_separated

end <- Sys.time()

(elapsed <- end - start)