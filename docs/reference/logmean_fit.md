# Fit a mixture of Gaussian Distributions to log mean count of taxa.

The optimal number of components to fit is chosen using parametric
bootstrap method

## Usage

``` r
logmean_fit(logmean, sig = 0.05, max.comp = 4, max.boot = 100, verb = FALSE)
```

## Arguments

- logmean:

  vector of log mean count of taxa

- sig:

  significance level to compare against p-value to be used for
  parametric bootstrap calculation

- max.comp:

  maximum number of Gaussian components to compare sequentially

- max.boot:

  maximum number of bootstraps simulations

- verb:

  If TRUE, it prints out updates of iterations of the algorithm

## Value

A list containing the optimal number of Gaussian components fitted; and
mean and variance parameter estimates from the fit

## Examples

``` r
# \donttest{
logmean  = rnorm(100)
logmean_fit(logmean,sig=0.05,max.comp=4,max.boot=100)
#> number of iterations= 82 
#> number of iterations= 40 
#> number of iterations= 40 
#> number of iterations= 115 
#> number of iterations= 26 
#> number of iterations= 84 
#> number of iterations= 15 
#> number of iterations= 50 
#> number of iterations= 11 
#> number of iterations= 8 
#> number of iterations= 27 
#> number of iterations= 28 
#> number of iterations= 13 
#> number of iterations= 37 
#> number of iterations= 18 
#> number of iterations= 21 
#> number of iterations= 31 
#> number of iterations= 15 
#> number of iterations= 74 
#> number of iterations= 107 
#> number of iterations= 22 
#> number of iterations= 9 
#> number of iterations= 25 
#> number of iterations= 47 
#> number of iterations= 26 
#> number of iterations= 20 
#> number of iterations= 8 
#> number of iterations= 38 
#> number of iterations= 84 
#> number of iterations= 15 
#> number of iterations= 86 
#> number of iterations= 25 
#> number of iterations= 103 
#> number of iterations= 20 
#> number of iterations= 50 
#> number of iterations= 287 
#> number of iterations= 25 
#> number of iterations= 28 
#> number of iterations= 33 
#> number of iterations= 27 
#> number of iterations= 75 
#> number of iterations= 19 
#> number of iterations= 11 
#> number of iterations= 147 
#> number of iterations= 9 
#> number of iterations= 260 
#> number of iterations= 4 
#> number of iterations= 25 
#> number of iterations= 125 
#> number of iterations= 42 
#> number of iterations= 61 
#> number of iterations= 19 
#> number of iterations= 115 
#> number of iterations= 22 
#> number of iterations= 44 
#> number of iterations= 86 
#> number of iterations= 24 
#> number of iterations= 12 
#> number of iterations= 69 
#> number of iterations= 24 
#> number of iterations= 38 
#> number of iterations= 17 
#> number of iterations= 11 
#> number of iterations= 14 
#> number of iterations= 48 
#> number of iterations= 105 
#> number of iterations= 26 
#> number of iterations= 15 
#> number of iterations= 133 
#> number of iterations= 20 
#> number of iterations= 176 
#> number of iterations= 166 
#> number of iterations= 9 
#> number of iterations= 20 
#> number of iterations= 25 
#> number of iterations= 26 
#> number of iterations= 8 
#> number of iterations= 12 
#> number of iterations= 5 
#> number of iterations= 117 
#> number of iterations= 75 
#> number of iterations= 104 
#> number of iterations= 93 
#> number of iterations= 16 
#> number of iterations= 99 
#> number of iterations= 15 
#> number of iterations= 83 
#> number of iterations= 69 
#> number of iterations= 72 
#> number of iterations= 60 
#> number of iterations= 11 
#> number of iterations= 22 
#> number of iterations= 45 
#> number of iterations= 40 
#> number of iterations= 29 
#> number of iterations= 14 
#> number of iterations= 18 
#> number of iterations= 18 
#> number of iterations= 27 
#> number of iterations= 270 
#> number of iterations= 86 

#> Decision: Select 1 Component(s) 
#> $param
#>       sigma        mu
#> 1 0.9474994 0.0395493
#> 
#> $components
#> [1] 1
#> 
# }
```
