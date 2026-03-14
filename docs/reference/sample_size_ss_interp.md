# Estimate sample size required to achieve a target statistical power

This function estimates the sample size required to achieve a specified
target power using predictions from a fitted GAM/SCAM power model. The
function evaluates predicted power across a grid of candidate sample
sizes and identifies the smallest sample size for which the predicted
power reaches or exceeds the target value. Linear interpolation is then
used on the log2(sample size) scale to obtain a more precise estimate.

## Usage

``` r
sample_size_ss_interp(
  target_power,
  logmean,
  abs_lfc,
  model,
  xmin = log2(5),
  xmax = log2(500),
  ngrid = 1000
)
```

## Arguments

- target_power:

  Numeric value specifying the desired statistical power.

- logmean:

  Numeric value representing the log mean abundance.

- abs_lfc:

  Numeric value representing the absolute log fold change.

- model:

  A fitted GAM/SCAM model used to predict statistical power.

- xmin:

  Numeric value giving the minimum log2(sample size) considered in the
  search. Default is `log2(5)`.

- xmax:

  Numeric value giving the maximum log2(sample size) considered in the
  search. Default is `log2(500)`.

- ngrid:

  Integer specifying the number of grid points used when searching for
  the sample size solution. Default is `1000`.

## Value

A numeric value representing the estimated sample size required to
achieve the target power. Returns `NA` if the target power is not
reached within the specified search range.

## Details

The function first constructs a grid of candidate sample sizes on the
log2 scale between `xmin` and `xmax`. Predicted power values are then
obtained from the fitted model for each grid point. The smallest sample
size at which the predicted power reaches the target value is
identified, and linear interpolation is used to refine the estimate.
