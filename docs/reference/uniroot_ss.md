# Sample Size estimation function using uniroot

Sample Size estimation function using uniroot

## Usage

``` r
uniroot_ss(
  target_power,
  logmean,
  abs_lfc,
  model,
  xmin = log2(10),
  xmax = log2(5000),
  maxiter = 10000,
  max_report = 2000
)
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

- maxiter:

  maximum number of iterations

- max_report:

  maximum group sample size to be predicted. Any predicted sample size
  that exceed max_report will be reported as "\> max_report"

## Value

A numeric value corresponding to the estimated sample size required to
achieve the target power.
