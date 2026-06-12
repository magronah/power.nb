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
# \donttest{
logmean  = rnorm(100)
optimal.comp(logmean,sig=0.05,max.comp=4,max.boot=100)
#> number of iterations= 303 
#> number of iterations= 16 
#> number of iterations= 25 
#> One of the variances is going to zero;  trying new starting values.
#> number of iterations= 122 
#> number of iterations= 8 
#> number of iterations= 113 
#> number of iterations= 75 
#> number of iterations= 46 
#> number of iterations= 8 
#> number of iterations= 169 
#> number of iterations= 98 
#> number of iterations= 112 
#> number of iterations= 117 
#> number of iterations= 8 
#> number of iterations= 41 
#> number of iterations= 10 
#> number of iterations= 102 
#> number of iterations= 38 
#> number of iterations= 14 
#> number of iterations= 70 
#> number of iterations= 60 
#> number of iterations= 99 
#> number of iterations= 11 
#> number of iterations= 8 
#> number of iterations= 80 
#> number of iterations= 45 
#> number of iterations= 96 
#> number of iterations= 28 
#> number of iterations= 8 
#> number of iterations= 189 
#> number of iterations= 35 
#> number of iterations= 304 
#> number of iterations= 67 
#> number of iterations= 45 
#> number of iterations= 22 
#> number of iterations= 2 
#> number of iterations= 90 
#> number of iterations= 3 
#> number of iterations= 91 
#> number of iterations= 28 
#> number of iterations= 46 
#> number of iterations= 36 
#> number of iterations= 11 
#> number of iterations= 39 
#> number of iterations= 23 
#> number of iterations= 236 
#> number of iterations= 13 
#> number of iterations= 45 
#> number of iterations= 77 
#> number of iterations= 22 
#> number of iterations= 64 
#> number of iterations= 36 
#> number of iterations= 19 
#> number of iterations= 178 
#> number of iterations= 38 
#> number of iterations= 27 
#> number of iterations= 59 
#> number of iterations= 7 
#> number of iterations= 25 
#> number of iterations= 19 
#> number of iterations= 12 
#> number of iterations= 15 
#> number of iterations= 11 
#> number of iterations= 60 
#> number of iterations= 5 
#> number of iterations= 65 
#> number of iterations= 95 
#> number of iterations= 12 
#> number of iterations= 65 
#> number of iterations= 3 
#> number of iterations= 7 
#> number of iterations= 14 
#> number of iterations= 62 
#> number of iterations= 17 
#> number of iterations= 38 
#> number of iterations= 12 
#> number of iterations= 24 
#> number of iterations= 12 
#> number of iterations= 49 
#> number of iterations= 53 
#> number of iterations= 27 
#> number of iterations= 7 
#> number of iterations= 83 
#> number of iterations= 23 
#> number of iterations= 6 
#> number of iterations= 76 
#> number of iterations= 22 
#> number of iterations= 8 
#> number of iterations= 11 
#> number of iterations= 147 
#> number of iterations= 211 
#> number of iterations= 17 
#> number of iterations= 33 
#> number of iterations= 65 
#> number of iterations= 3 
#> number of iterations= 51 
#> number of iterations= 27 
#> number of iterations= 20 
#> number of iterations= 271 
#> number of iterations= 135 
#> number of iterations= 27 

#> Decision: Select 1 Component(s) 
#> [1] 1
# }
```
