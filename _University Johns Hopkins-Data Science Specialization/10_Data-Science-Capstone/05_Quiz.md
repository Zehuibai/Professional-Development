Task 05: Quiz
================
Mark Blackmore
November 8, 2017

For each of the sentence fragments below use your natural language processing algorithm to predict the next word in the sentence.

Q1. When you breathe, I want to be the air for you. I'll be there for you, I'd live and I'd

-   sleep
-   give
-   die
-   eat

``` r
ngrams("When you breathe, I want to be the air for you. I'll be there for you, I'd live and I'd")
```

    ## [1] "like"

Q2. Guy at my table's wife got up to go to the bathroom and I asked about dessert and he started telling me about his

-   horticultural
-   marital
-   spiritual
-   financial

``` r
ngrams("Guy at my table's wife got up to go to the bathroom and I asked about dessert and he started telling me about his")
```

    ## [1] "own"

Q3. I'd give anything to see arctic monkeys this

-   weekend
-   decade
-   month
-   morning

``` r
ngrams("I'd give anything to see arctic monkeys this")
```

    ## [1] "is"

Q4. Talking to your mom has the same effect as a hug and helps reduce your

-   sleepiness
-   happiness
-   hunger
-   stress

``` r
ngrams("Talking to your mom has the same effect as a hug and helps reduce your")
```

    ## [1] "own"

Q5. When you were in Holland you were like 1 inch away from me but you hadn't time to take a

-   walk
-   minute
-   look
-   picture

``` r
ngrams("When you were in Holland you were like 1 inch away from me but you hadn't time to take a")
```

    ## [1] "break"

Q6. I'd just like all of these questions answered, a presentation of evidence, and a jury to settle the

-   account
-   case
-   matter
-   incident

``` r
ngrams("I'd just like all of these questions answered, a presentation of evidence, and a jury to settle the")
```

    ## [1] "of"

Q7. I can't deal with unsymetrical things. I can't even hold an uneven number of bags of groceries in each

-   hand
-   finger
-   toe
-   arm

``` r
ngrams("I can't deal with unsymetrical things. I can't even hold an uneven number of bags of groceries in each")
```

    ## [1] "of"

Q8. Every inch of you is perfect from the bottom to the

-   top
-   middle
-   center
-   side

``` r
ngrams("Every inch of you is perfect from the bottom to the")
```

    ## [1] "of"

Q9. I'm thankful my childhood was filled with imagination and bruises from playing

-   weekly
-   inside
-   outside
-   daily

``` r
ngrams("I'm thankful my childhood was filled with imagination and bruises from playing ")
```

    ## [1] "with"

Q10. I like how the same people are in almost all of Adam Sandler's

-   movies
-   stories
-   pictures
-   novels

``` r
ngrams("I like how the same people are in almost all of Adam Sandler's")
```

    ## [1] "?"
