# Solve for the sample size required to achieve a target statistical power

This function estimates the sample size required to achieve a specified
target statistical power using a fitted GAM/SCAM power model. The
function first attempts to solve for the sample size using a
root-finding algorithm. If the root-finding procedure fails, a
grid-based interpolation method is used as a fallback to obtain an
approximate solution.

## Usage

``` r
ss_solver(
  target_power,
  logmean,
  abs_lfc,
  model,
  xmin = log2(5),
  xmax = log2(500)
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

## Value

A numeric value representing the estimated sample size required to
achieve the target power.

## Details

The function first attempts to compute the required sample size using a
root-finding approach implemented in
[`uniroot_ss()`](https://michaelagronah.com/power.nb/reference/uniroot_ss.md).
If this method fails (for example, due to numerical issues or if the
root cannot be bracketed within the specified interval), the function
falls back to a grid-based interpolation approach implemented in
[`sample_size_ss_interp()`](https://michaelagronah.com/power.nb/reference/sample_size_ss_interp.md).

A warning is issued if the target power is specified as 0 or 1, since
these values correspond to unrealistic design targets in most practical
applications.
