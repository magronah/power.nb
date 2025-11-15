# Fit the non-linear function to dispersion estimates

Dispersion are estimated from the `DESeq2` package. The function fitted
is of the form `a + b/(mean count)` where `a` represents the asymptotic
dispersion level for high abundance taxa, and `b` captures additional
dispersion variability.

## Usage

``` r
dispersion_fit(dispersion, logmean)
```

## Arguments

- dispersion:

  dispersion estimates from deseq

- logmean:

  vector of log mean abundance

## Value

A list containing estimates for `a` and `b` and confidence intervals

## Examples

``` r
logmean    =  rnorm(100)
dispersion =  abs(rnorm(100))
dispersion_fit(dispersion,logmean)
#> Warning: singular gradient
#> Waiting for profiling to be done...
#> $param
#>   asymptDisp extraPois
#> 1   1.040575         0
#> 
#> $confint
#>            2.5 % 97.5 %
#> asymptDisp    NA     NA
#> extraPois     NA     NA
#> 
```
