# General-purpose log-likelihood function, vectorized sum(pars\*x^i)

General-purpose log-likelihood function, vectorized sum(pars\*x^i)

## Usage

``` r
polyfun(pars, x)
```

## Arguments

- pars:

  parameters

- x:

  log mean count

## Value

values representing output for the polynomial fuction (f(x))

## Examples

``` r
polyfun(pars = c(1, 2, 3), x = 1:5)
#> [1]  6 17 34 57 86
polyfun(pars = c(1, 0, 3), x = 1)
#> [1] 4
```
