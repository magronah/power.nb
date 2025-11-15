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
# Define simulated log mean abundance
logmean_sim <- rnorm(100, mean = 0, sd = 1)

# Define parameters for log fold change simulation
logfoldchange_param <- list(
  par = c(1, -0.5, 0.2), # Example parameters
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
#> Error in genmixpars(logmean, par, ...): number of pars (3) != expected (11) (np = 2, sd_ord = 2)
```
