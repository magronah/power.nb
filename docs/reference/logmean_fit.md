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
#> number of iterations= 24 
#> number of iterations= 19 
#> number of iterations= 30 
#> number of iterations= 15 
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 52 
#> number of iterations= 29 
#> number of iterations= 13 
#> number of iterations= 17 
#> number of iterations= 15 
#> number of iterations= 9 
#> number of iterations= 6 
#> number of iterations= 63 
#> number of iterations= 21 
#> number of iterations= 46 
#> number of iterations= 131 
#> number of iterations= 42 
#> number of iterations= 13 
#> number of iterations= 9 
#> number of iterations= 35 
#> number of iterations= 32 
#> number of iterations= 102 
#> number of iterations= 117 
#> number of iterations= 74 
#> number of iterations= 7 
#> number of iterations= 163 
#> number of iterations= 97 
#> number of iterations= 47 
#> number of iterations= 25 
#> number of iterations= 21 
#> number of iterations= 17 
#> number of iterations= 18 
#> number of iterations= 40 
#> number of iterations= 114 
#> number of iterations= 19 
#> number of iterations= 37 
#> number of iterations= 8 
#> number of iterations= 47 
#> number of iterations= 16 
#> number of iterations= 153 
#> number of iterations= 7 
#> number of iterations= 10 
#> number of iterations= 8 
#> number of iterations= 130 
#> number of iterations= 36 
#> number of iterations= 41 
#> number of iterations= 11 
#> number of iterations= 10 
#> number of iterations= 13 
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 139 
#> number of iterations= 132 
#> number of iterations= 6 
#> number of iterations= 115 
#> number of iterations= 59 
#> number of iterations= 116 
#> number of iterations= 15 
#> number of iterations= 111 
#> number of iterations= 27 
#> number of iterations= 56 
#> number of iterations= 14 
#> number of iterations= 30 
#> number of iterations= 48 
#> number of iterations= 27 
#> number of iterations= 138 
#> number of iterations= 102 
#> number of iterations= 12 
#> number of iterations= 81 
#> number of iterations= 24 
#> number of iterations= 103 
#> number of iterations= 105 
#> number of iterations= 3 
#> number of iterations= 26 
#> number of iterations= 21 
#> number of iterations= 32 
#> number of iterations= 3 
#> number of iterations= 25 
#> number of iterations= 88 
#> number of iterations= 16 
#> number of iterations= 41 
#> number of iterations= 9 
#> number of iterations= 8 
#> number of iterations= 95 
#> number of iterations= 109 
#> number of iterations= 16 
#> number of iterations= 83 
#> number of iterations= 27 
#> number of iterations= 22 
#> number of iterations= 32 
#> number of iterations= 84 
#> number of iterations= 7 
#> number of iterations= 26 
#> number of iterations= 14 
#> number of iterations= 30 
#> number of iterations= 14 
#> number of iterations= 32 
#> number of iterations= 52 
#> number of iterations= 46 
#> number of iterations= 230 
#> number of iterations= 60 
#> number of iterations= 19 
#> number of iterations= 62 
#> number of iterations= 164 

#> Decision: Select 1 Component(s) 
#> $param
#>       sigma          mu
#> 1 0.8854284 -0.09361792
#> 
#> $components
#> [1] 1
#> 
```
