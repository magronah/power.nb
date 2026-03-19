# Computes the optimal number of gaussian components for log mean count

The number of gaussian components is determined using the using
parametric bootstrap

## Usage

``` r
optimal.comp(logmean, sig = 0.05, max.comp = 4, max.boot = 100)
```

## Arguments

- logmean:

  vector of log mean abundances of taxa

- sig:

  significance level to compare against p-value

- max.comp:

  maximum number of Gaussian components to compare sequentially

- max.boot:

  maximum number of bootstraps simulations

## Value

The best number of components for fitting the distribution of log mean
abundance

## Examples

``` r
logmean  = rnorm(100)
optimal.comp(logmean,sig=0.05,max.comp=4,max.boot=100)
#> number of iterations= 20 
#> number of iterations= 186 
#> number of iterations= 106 
#> number of iterations= 51 
#> number of iterations= 2 
#> number of iterations= 48 
#> number of iterations= 12 
#> number of iterations= 13 
#> number of iterations= 108 
#> number of iterations= 26 
#> number of iterations= 31 
#> number of iterations= 13 
#> number of iterations= 11 
#> number of iterations= 97 
#> number of iterations= 60 
#> number of iterations= 367 
#> number of iterations= 48 
#> number of iterations= 71 
#> number of iterations= 204 
#> number of iterations= 15 
#> number of iterations= 16 
#> number of iterations= 177 
#> number of iterations= 65 
#> number of iterations= 139 
#> number of iterations= 17 
#> number of iterations= 72 
#> number of iterations= 64 
#> number of iterations= 2 
#> number of iterations= 15 
#> number of iterations= 22 
#> number of iterations= 18 
#> number of iterations= 19 
#> number of iterations= 364 
#> number of iterations= 150 
#> number of iterations= 16 
#> number of iterations= 21 
#> number of iterations= 56 
#> number of iterations= 53 
#> number of iterations= 45 
#> number of iterations= 115 
#> number of iterations= 10 
#> number of iterations= 12 
#> number of iterations= 38 
#> number of iterations= 70 
#> number of iterations= 37 
#> number of iterations= 21 
#> number of iterations= 37 
#> number of iterations= 40 
#> number of iterations= 52 
#> number of iterations= 60 
#> number of iterations= 218 
#> number of iterations= 20 
#> number of iterations= 21 
#> number of iterations= 8 
#> number of iterations= 157 
#> number of iterations= 30 
#> number of iterations= 43 
#> number of iterations= 11 
#> number of iterations= 19 
#> number of iterations= 58 
#> number of iterations= 20 
#> number of iterations= 49 
#> number of iterations= 6 
#> number of iterations= 114 
#> number of iterations= 23 
#> number of iterations= 19 
#> number of iterations= 9 
#> number of iterations= 39 
#> number of iterations= 152 
#> number of iterations= 31 
#> number of iterations= 21 
#> number of iterations= 12 
#> number of iterations= 82 
#> number of iterations= 148 
#> number of iterations= 12 
#> number of iterations= 2 
#> number of iterations= 26 
#> number of iterations= 65 
#> number of iterations= 15 
#> number of iterations= 45 
#> number of iterations= 12 
#> number of iterations= 36 
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 19 
#> number of iterations= 160 
#> number of iterations= 23 
#> number of iterations= 38 
#> number of iterations= 10 
#> number of iterations= 14 
#> number of iterations= 26 
#> number of iterations= 62 
#> number of iterations= 13 
#> number of iterations= 10 
#> number of iterations= 19 
#> number of iterations= 20 
#> number of iterations= 19 
#> number of iterations= 25 
#> number of iterations= 10 
#> number of iterations= 29 
#> number of iterations= 141 
#> number of iterations= 42 
#> number of iterations= 21 

#> Decision: Select 1 Component(s) 
#> [1] 1
```
