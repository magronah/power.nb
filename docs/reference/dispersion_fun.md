# Calculate Dispersion for Microbiome Data

This function calculates the dispersion value for microbiome data based
on the provided parameters: mean abundance, asymptotic dispersion, and
extra Poisson dispersion.

## Usage

``` r
dispersion_fun(mean_abund, asymptDisp, extraPois)
```

## Arguments

- mean_abund:

  Numeric value representing the mean abundance of the taxa.

- asymptDisp:

  Numeric value for the asymptotic dispersion (the dispersion at high
  abundance).

- extraPois:

  Numeric value for the extra Poisson dispersion (to model
  overdispersion).

## Value

A numeric value representing the dispersion.

## Details

The dispersion is calculated using the formula: \$\$\text{dispersion} =
\text{asymptDisp} + \frac{\text{extraPois}}{\text{mean_abund}}\$\$

## Examples

``` r
mean_abund <- 10
asymptDisp <- 0.1
extraPois <- 0.05
dispersion_fun(mean_abund, asymptDisp, extraPois)
#> [1] 0.105
```
