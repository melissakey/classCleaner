
<!-- README.md is generated from README.Rmd. Please edit that file -->

# classCleaner

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The goal of classCleaner is to …

## Installation

<!-- You can install the released version of classCleaner from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("classCleaner") -->

<!-- ``` -->

`classCleaner` is currently only available on
[Github](https://github.com/), where it can be downloaded at
<!-- And the development version from [GitHub](https://github.com/) with: -->

``` r
# install.packages("devtools")
devtools::install_github("melkey/classCleaner")
```

## Example

Here’s a simple example of `classCleaner` on the Congressional Voting
Records Data Set from the UCI Machine Learning Repository. Although
`classCleaner` does not require any tidyverse package, it’s included to
clean and present the data.

The data consists of (simplified) voting records for each member of the
US House of Representatives on 16 key votes (identified by the CQA) in
1984. Each Representative’s vote is summarized as “y” (voted for), “n”
(voted against), and “?” (voted present or did not vote). We use
`classCleaner` to identify which party members tended to vote in a
manner different from their political party.

``` r
library(classCleaner)
library(tidyverse)
#> Warning: package 'tidyr' was built under R version 3.6.3
#> Warning: package 'dplyr' was built under R version 3.6.3
#> Warning: package 'forcats' was built under R version 3.6.3

## load example data set

votes <- read_csv(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/voting-records/house-votes-84.data",
  col_types = cols(),
  col_names = c('party', paste0("v", 1:16))
) %>%
  print()
#> # A tibble: 435 x 17
#>    party v1    v2    v3    v4    v5    v6    v7    v8    v9    v10   v11   v12  
#>    <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#>  1 repu~ n     y     n     y     y     y     n     n     n     y     ?     y    
#>  2 repu~ n     y     n     y     y     y     n     n     n     n     n     y    
#>  3 demo~ ?     y     y     ?     y     y     n     n     n     n     y     n    
#>  4 demo~ n     y     y     n     ?     y     n     n     n     n     y     n    
#>  5 demo~ y     y     y     n     y     y     n     n     n     n     y     ?    
#>  6 demo~ n     y     y     n     y     y     n     n     n     n     n     n    
#>  7 demo~ n     y     n     y     y     y     n     n     n     n     n     n    
#>  8 repu~ n     y     n     y     y     y     n     n     n     n     n     n    
#>  9 repu~ n     y     n     y     y     y     n     n     n     n     n     y    
#> 10 demo~ y     y     y     n     n     n     y     y     y     n     n     n    
#> # ... with 425 more rows, and 4 more variables: v13 <chr>, v14 <chr>,
#> #   v15 <chr>, v16 <chr>
```

`ClassCleaner`, requires a matrix of pairwise distances between the
Representatives. To do this, we recode the `v1` - `v16`, with `y -> 1`,
`n -> 0`, and `? -> 0.5`, and use Manhattan (“city-block”) distance to
calculate the distance between Representatives.

``` r
votes <- votes %>%
  mutate_at(vars(v1:v16), 
    ~ recode(., 
      'y' = 1,
      'n' = 0,
      '?' = .5
    )
  ) %>%
  print()
#> # A tibble: 435 x 17
#>    party    v1    v2    v3    v4    v5    v6    v7    v8    v9   v10   v11   v12
#>    <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 repu~   0       1     0   1     1       1     0     0     0     1   0.5   1  
#>  2 repu~   0       1     0   1     1       1     0     0     0     0   0     1  
#>  3 demo~   0.5     1     1   0.5   1       1     0     0     0     0   1     0  
#>  4 demo~   0       1     1   0     0.5     1     0     0     0     0   1     0  
#>  5 demo~   1       1     1   0     1       1     0     0     0     0   1     0.5
#>  6 demo~   0       1     1   0     1       1     0     0     0     0   0     0  
#>  7 demo~   0       1     0   1     1       1     0     0     0     0   0     0  
#>  8 repu~   0       1     0   1     1       1     0     0     0     0   0     0  
#>  9 repu~   0       1     0   1     1       1     0     0     0     0   0     1  
#> 10 demo~   1       1     1   0     0       0     1     1     1     0   0     0  
#> # ... with 425 more rows, and 4 more variables: v13 <dbl>, v14 <dbl>,
#> #   v15 <dbl>, v16 <dbl>

dist_mat <- votes %>% 
  as.data.frame(select(., -1)) %>%
  as.matrix() %>%
  dist(method = "manhattan") %>%
  as.matrix()
```

Running classCleaner on this produces

``` r
cc_result <- classCleaner(dist_mat, votes$party)

print(cc_result)
#> classCleaner was applied to 2 classes
#> A total of 362 instances were retained out of 435 
#> 
#>            # Instances # Retained    %
#> democrat           267        217 81.3
#> republican         168        145 86.3
```
