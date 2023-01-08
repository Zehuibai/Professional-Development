Coursera Data Science Specialization Capstone Project
========================================================
author: Mark Blackmore
date: 2017-11-19
autosize: true

The Project
========================================================

This project involves Natural Language Processing.  The critical task is to 
take a user's input phrase (group of words) and to output a predicted next word.  

*Project deliverables:*  

- Next Word Prediction Model, as basis for an app
- Next Word Prediction App hosted at shinyapps.io
- This presentation hosted at R pubs

Next Word Prediction Model
========================================================

The next word prediction model uses the principles of "tidy data" applied to text mining in R. Key model steps: 

1. Input: raw text files for model training
2. Clean training data; separate into 2 word, 3 word, and 4 word n grams, save as tibbles
3. Sort n grams tibbles by frequency, save as repos
4. N grams function: uses a "back-off" type prediction model
  - user supplies an input phrase
  - model uses last 3, 2, or 1 words to predict the best 4th, 3rd, or 2nd match in the repos
5. Output: next word prediction

Benefits: easy to read code; uses "pipes"; fast processing of training data; able to sample up to 25% of original corpus; relatively small output repos

Next Word Prediction App
========================================================

The next word prediction app provides a simple user interface to the next word prediction model.  

*Key Features:*  

1. Text box for user input  
2. Predicted next word outputs dynamically below user input  
3. Tabs with plots of most frequent n grams in the data-set
4. Side panel with user instructions  

*Key Benefits:*  

1. Fast response  
2. Method allows for large training sets leading to better next word predictions

[Shiny App Link](https://mblackmo.shinyapps.io/ngram_match/)  

Documentation and Source Code
========================================================

Tidy Data  
"http://vita.had.co.nz/papers/tidy-data.html"

Text Mining with R: A Tidy Approach  
"http://tidytextmining.com/index.html"

Shiny App    
"https://mblackmo.shinyapps.io/ngram_match/"

Shiny App Source Code repository on Github    
"https://github.com/mark-blackmore/JHU-Data-Science-Capstone/tree/master/ngram_match"

Data Specialization Capstone repository on Github    
"https://github.com/mark-blackmore/JHU-Data-Science-Capstone"

