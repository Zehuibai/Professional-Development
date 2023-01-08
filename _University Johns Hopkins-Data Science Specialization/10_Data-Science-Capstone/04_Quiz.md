Task 04: Quiz
================
Mark Blackmore
October 30, 2017

For each of the sentence fragments below use your natural language processing algorithm to predict the next word in the sentence.

Q1. The guy in front of me just bought a pound of bacon, a bouquet, and a case of

-   cheese
-   pretzels
-   beer
-   soda

``` r
ngrams("The guy in front of me just bought a pound of bacon, a bouquet, and a case of")  
```

    ## [1] "the"

Q2. You're the reason why I smile everyday. Can you follow me please? It would mean the

-   universe
-   most
-   world
-   best

``` r
ngrams("You're the reason why I smile everyday. Can you follow me please? It would mean the")  
```

    ## [1] "world"

Q3. Hey sunshine, can you follow me and make me the

-   happiest
-   bluest
-   saddest
-   smelliest

``` r
ngrams("Hey sunshine, can you follow me and make me the")  
```

    ## [1] "most"

Q4. Very early observations on the Bills game: Offense still struggling but the

-   referees
-   players
-   crowd
-   defense

``` r
ngrams("Very early observations on the Bills game: Offense still struggling but the") 
```

    ## [1] "is"

Q5. Go on a romantic date at the

-   mall
-   movies
-   beach
-   grocery

``` r
ngrams("Go on a romantic date at the") 
```

    ## [1] "end"

Q6. Well I'm pretty sure my granny has some old bagpipes in her garage I'll dust them off and be on my

-   horse
-   phone
-   way
-   motorcycle

``` r
ngrams("Well I'm pretty sure my granny has some old bagpipes in her garage I'll dust them off and be on my") 
```

    ## [1] "way"

Q7. Ohhhhh \#PointBreak is on tomorrow. Love that film and haven't seen it in quite some

-   time
-   thing
-   weeks
-   years

``` r
ngrams("Ohhhhh #PointBreak is on tomorrow. Love that film and haven't seen it in quite some") 
```

    ## [1] "time"

Q8. After the ice bucket challenge Louis will push his long wet hair out of his eyes with his little

-   toes
-   ears
-   eyes
-   fingers

``` r
ngrams("After the ice bucket challenge Louis will push his long wet hair out of his eyes with his little") 
```

    ## [1] "?"

Q9. Be grateful for the good times and keep the faith during the

-   hard
-   bad
-   sad
-   worse

``` r
ngrams("Be grateful for the good times and keep the faith during the") 
```

    ## [1] "day"

Q10. If this isn't the cutest thing you've ever seen, then you must be

-   insane
-   asleep
-   callous
-   insensitive

``` r
ngrams("If this isn't the cutest thing you've ever seen, then you must be") 
```

    ## [1] "a"
