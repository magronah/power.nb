# Simulate Log Fold Change Values

This function generates simulated log fold change (LFC) values based on
the provided log mean abundance and LFC parameters. The simulation
ensures that the generated LFC values remain within a specified maximum
range by iterating until convergence or until a maximum iteration limit
is reached.

## Usage

``` r
logfoldchange_sim_fun(
  logmean_sim,
  logfoldchange_param,
  max_lfc = 15,
  max_iter = 10000,
  seed = 121
)
```

## Arguments

- logmean_sim:

  A numeric vector of simulated log mean abundances.

- logfoldchange_param:

  A list containing parameters for the log fold change simulation:

  - `par`: Optimal parameters for the log fold change fit.

  - `np`: Optimal number of components for the log fold change model.

  - `sd_ord`: Order of the polynomial used for the standard deviation
    parameter of the log fold change.

- max_lfc:

  A numeric value specifying the maximum allowable absolute log fold
  change value. Default is 15.

- max_iter:

  An integer specifying the maximum number of iterations allowed to
  ensure all simulated LFC values are within the `max_lfc` range.
  Default is 10,000.

- seed:

  random-number seed

## Value

A numeric vector of simulated log fold change values (`lfc`).

## Examples

``` r
set.seed(101)
# Define simulated log mean abundance
logmean_sim <- rnorm(100, mean = 0, sd = 1)

# Define parameters for log fold change simulation
logfoldchange_param <- list(
  par = rnorm(11),       # Example parameters
  np = 2,                # Number of components
  sd_ord = 2             # Order of polynomial for SD
)

# Simulate log fold change values
logfoldchange_sim_fun(
  logmean_sim = logmean_sim,
  logfoldchange_param = logfoldchange_param,
  max_lfc = 10,
  max_iter = 5000
)
#>   [1]  0.77536534 -0.75927908  0.57559853  0.24988592  1.29554742  3.20377118
#>   [7]  0.47951965 -0.39424034  3.08421259  3.45909061  0.07150158  1.74829991
#>  [13] -1.14778645  1.44999231  2.08759275  2.67255845  1.42411340 -0.94717905
#>  [19]  4.33031343  0.70459024  0.45920383  3.14750041  2.02422422  1.23460892
#>  [25]  2.43004015  1.54934384  2.46467222  1.99330437 -0.45609681  3.11060255
#>  [31]  0.16023216  2.57254352  3.25793035  0.60396932  1.10901345  1.37469207
#>  [37]  2.92895477  0.20506906  4.04637355  2.02622178 -0.43806886  2.50240781
#>  [43]  0.29248446  0.64951982  0.78484957  1.50003440 -0.01191195 -5.98171420
#>  [49] -2.94150866  2.62317835  1.30759331  0.71860840  0.68786991  0.61950167
#>  [55]  2.03931987  1.75635258  2.97146597  3.13447579  2.77627194  2.49630360
#>  [61]  0.09310071  1.06172800  1.06209385 -0.17932371  2.48775509  2.70721117
#>  [67]  3.06653638  0.52568456  2.00110377 -2.77202494  2.35125799  3.37302253
#>  [73]  0.28443019  1.71850784  1.06309408  1.29522830  2.30808053  0.35606973
#>  [79] -0.19398176 -0.96092257  3.57141184  0.67600711  2.93379296  2.35866360
#>  [85]  0.04805865  2.67475960  2.91065453  0.88217677  2.54776601  3.02793378
#>  [91]  3.51039604  1.25352899 -5.75334804 -1.21701458 -1.37136417  1.56562167
#>  [97] -0.64715531  0.95082237 -0.18732567  1.68143967
```
