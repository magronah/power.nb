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
#> number of iterations= 27 
#> number of iterations= 51 
#> number of iterations= 10 
#> number of iterations= 8 
#> number of iterations= 38 
#> number of iterations= 8 
#> number of iterations= 31 
#> number of iterations= 56 
#> number of iterations= 26 
#> number of iterations= 81 
#> number of iterations= 88 
#> number of iterations= 43 
#> number of iterations= 27 
#> number of iterations= 103 
#> number of iterations= 19 
#> number of iterations= 66 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 34 
#> number of iterations= 96 
#> number of iterations= 17 
#> number of iterations= 19 
#> number of iterations= 11 
#> number of iterations= 31 
#> number of iterations= 25 
#> number of iterations= 4 
#> number of iterations= 38 
#> number of iterations= 16 
#> number of iterations= 61 
#> number of iterations= 29 
#> number of iterations= 6 
#> number of iterations= 51 
#> number of iterations= 31 
#> number of iterations= 44 
#> number of iterations= 20 
#> number of iterations= 11 
#> number of iterations= 50 
#> number of iterations= 52 
#> number of iterations= 22 
#> number of iterations= 45 
#> number of iterations= 99 
#> number of iterations= 35 
#> number of iterations= 10 
#> number of iterations= 21 
#> number of iterations= 24 
#> number of iterations= 3 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 36 
#> number of iterations= 256 
#> number of iterations= 39 
#> number of iterations= 17 
#> number of iterations= 31 
#> number of iterations= 85 
#> number of iterations= 28 
#> number of iterations= 64 
#> number of iterations= 20 
#> number of iterations= 45 
#> number of iterations= 82 
#> number of iterations= 27 
#> number of iterations= 41 
#> number of iterations= 128 
#> number of iterations= 63 
#> number of iterations= 109 
#> number of iterations= 96 
#> number of iterations= 34 
#> number of iterations= 59 
#> number of iterations= 23 
#> number of iterations= 23 
#> number of iterations= 22 
#> number of iterations= 7 
#> number of iterations= 112 
#> number of iterations= 39 
#> number of iterations= 10 
#> number of iterations= 25 
#> One of the variances is going to zero;  trying new starting values.
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 50 
#> number of iterations= 44 
#> number of iterations= 128 
#> number of iterations= 27 
#> number of iterations= 30 
#> number of iterations= 13 
#> number of iterations= 34 
#> number of iterations= 70 
#> number of iterations= 20 
#> number of iterations= 128 
#> number of iterations= 102 
#> number of iterations= 65 
#> number of iterations= 61 
#> number of iterations= 77 
#> number of iterations= 28 
#> number of iterations= 32 
#> number of iterations= 51 
#> number of iterations= 14 
#> number of iterations= 18 
#> number of iterations= 80 
#> number of iterations= 29 
#> number of iterations= 19 
#> number of iterations= 5 
#> number of iterations= 31 
#> number of iterations= 3 
#> number of iterations= 69 
#> number of iterations= 26 
#> number of iterations= 93 
#> number of iterations= 118 
#> number of iterations= 114 

#> Decision: Select 1 Component(s) 
#> [1] 1
```
