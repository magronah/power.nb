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
logmean  = rnorm(100)
logmean_fit(logmean,sig=0.05,max.comp=4,max.boot=100)
#> number of iterations= 143 
#> number of iterations= 21 
#> number of iterations= 151 
#> number of iterations= 89 
#> number of iterations= 15 
#> number of iterations= 62 
#> number of iterations= 12 
#> number of iterations= 26 
#> number of iterations= 28 
#> number of iterations= 176 
#> number of iterations= 50 
#> number of iterations= 69 
#> number of iterations= 22 
#> number of iterations= 52 
#> number of iterations= 37 
#> number of iterations= 48 
#> number of iterations= 24 
#> number of iterations= 130 
#> number of iterations= 28 
#> number of iterations= 35 
#> number of iterations= 48 
#> number of iterations= 75 
#> number of iterations= 42 
#> number of iterations= 136 
#> number of iterations= 23 
#> number of iterations= 265 
#> number of iterations= 39 
#> number of iterations= 79 
#> number of iterations= 18 
#> number of iterations= 11 
#> number of iterations= 31 
#> number of iterations= 22 
#> number of iterations= 2 
#> number of iterations= 147 
#> number of iterations= 65 
#> number of iterations= 135 
#> number of iterations= 3 
#> number of iterations= 13 
#> number of iterations= 83 
#> number of iterations= 180 
#> number of iterations= 19 
#> number of iterations= 46 
#> number of iterations= 138 
#> number of iterations= 78 
#> number of iterations= 6 
#> number of iterations= 25 
#> number of iterations= 12 
#> number of iterations= 26 
#> number of iterations= 55 
#> number of iterations= 8 
#> number of iterations= 13 
#> number of iterations= 168 
#> number of iterations= 120 
#> number of iterations= 123 
#> number of iterations= 10 
#> number of iterations= 143 
#> number of iterations= 46 
#> number of iterations= 14 
#> number of iterations= 101 
#> number of iterations= 50 
#> number of iterations= 4 
#> number of iterations= 18 
#> number of iterations= 89 
#> number of iterations= 30 
#> number of iterations= 59 
#> number of iterations= 9 
#> number of iterations= 136 
#> number of iterations= 194 
#> number of iterations= 137 
#> number of iterations= 50 
#> number of iterations= 24 
#> number of iterations= 13 
#> number of iterations= 54 
#> number of iterations= 131 
#> number of iterations= 87 
#> number of iterations= 76 
#> number of iterations= 87 
#> number of iterations= 126 
#> number of iterations= 23 
#> number of iterations= 15 
#> number of iterations= 23 
#> number of iterations= 85 
#> number of iterations= 48 
#> number of iterations= 26 
#> number of iterations= 20 
#> number of iterations= 9 
#> number of iterations= 81 
#> number of iterations= 43 
#> number of iterations= 40 
#> number of iterations= 57 
#> number of iterations= 27 
#> number of iterations= 2 
#> number of iterations= 18 
#> number of iterations= 23 
#> number of iterations= 27 
#> number of iterations= 3 
#> number of iterations= 11 
#> number of iterations= 31 
#> number of iterations= 29 
#> number of iterations= 41 
#> number of iterations= 150 

#> Decision: Select 1 Component(s) 
#> $param
#>       sigma        mu
#> 1 0.8854741 -0.096547
#> 
#> $components
#> [1] 1
#> 
```
