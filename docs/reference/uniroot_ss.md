# Sample Size estimation function uisng uniroot

Sample Size estimation function uisng uniroot

## Usage

``` r
uniroot_ss(target_power, logmean, abs_lfc, model, xmin, xmax)
```

## Arguments

- target_power:

  Numeric value specifying the desired statistical power.

- logmean:

  Numeric value representing the log of the mean abundance.

- abs_lfc:

  Numeric value representing the absolute log fold change.

- model:

  A fitted GAM/SCAM model used to predict statistical power.

- xmin:

  Numeric value giving the minimum sample size considered in the search.

- xmax:

  Numeric value giving the maximum sample size considered in the search.

## Value

A numeric value corresponding to the estimated sample size required to
achieve the target power.
