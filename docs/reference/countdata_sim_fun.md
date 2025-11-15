# Simulate Count Data for Microbiome Studies

This function simulates count data for microbiome studies based on log
mean, log fold change, and dispersion parameters. It supports generating
data for multiple simulations and allows flexibility in specifying the
number of control and treatment samples or samples per group.

## Usage

``` r
countdata_sim_fun(
  logmean_param,
  logfoldchange_param,
  dispersion_param,
  nsamp_per_group = NULL,
  ncont = NULL,
  ntreat = NULL,
  notu,
  nsim = 1,
  disp_scale = 0.3,
  max_lfc = 15,
  maxlfc_iter = 1000,
  seed = NULL
)
```

## Arguments

- logmean_param:

  A list of parameters for simulating the log mean abundance.

- logfoldchange_param:

  A list of parameters for simulating log fold change, containing:

  - `par`: Optimal parameters for log fold change fitting.

  - `np`: Number of components for the log fold change model.

  - `sd_ord`: Order of the polynomial for the standard deviation
    parameter.

- dispersion_param:

  A list of dispersion parameters containing:

  - `asymptDisp`: Asymptotic dispersion parameter.

  - `extraPois`: Additional Poisson variation parameter.

- nsamp_per_group:

  Number of samples per group (control and treatment). If provided,
  `ncont` and `ntreat` must not be specified.

- ncont:

  Number of control samples. Specify along with `ntreat` when
  `nsamp_per_group` is not provided.

- ntreat:

  Number of treatment samples. Specify along with `ncont` when
  `nsamp_per_group` is not provided.

- notu:

  Number of operational taxonomic units (OTUs) to simulate.

- nsim:

  Number of simulations to run. Default is 1.

- disp_scale:

  Scale parameter for the dispersion. Default is 0.3.

- max_lfc:

  Maximum allowable log fold change. Default is 15.

- maxlfc_iter:

  Maximum number of iterations for ensuring log fold change is within
  `max_lfc`. Default is 1,000.

- seed:

  Seed value for reproducibility. Default is `NULL`.

## Value

A list containing:

- `countdata_list`: A list of count data matrices for each simulation.

- `metadata_list`: A list of metadata data frames for each simulation.

- `logmean_list`: A list of log mean vectors for each simulation.

- `logfoldchange_list`: A list of log fold change vectors for each
  simulation.

- `treat_countdata_list`: A list of treatment count data matrices for
  each simulation.

- `control_countdata_list`: A list of control count data matrices for
  each simulation.

## Examples

``` r
# Load required packages
library(foreach)
library(doParallel)
#> Loading required package: iterators
#> Loading required package: parallel
# Define parameters
logmean_param <- list(mu = 0, sigma = 1)
logfoldchange_param <- list(par = rnorm(11), np = 2, sd_ord = 2)
dispersion_param <- list(asymptDisp = 0.1, extraPois = 0.05)

# Simulate count data
result <- countdata_sim_fun(
  logmean_param = logmean_param,
  logfoldchange_param = logfoldchange_param,
  dispersion_param = dispersion_param,
  nsamp_per_group = 10,
  notu = 50,
  nsim = 2,
  seed = 123
)
#> 
#> Attaching package: ‘purrr’
#> The following objects are masked from ‘package:foreach’:
#> 
#>     accumulate, when
#> mixtools package, version 2.0.0.1, Released 2022-12-04
#> This package is based upon work supported by the National Science Foundation under Grant No. SES-0518772 and the Chan Zuckerberg Initiative: Essential Open Source Software for Science (Grant No. 2020-255193).

# Access simulation results
countdata <- result$countdata_list[[1]]
metadata <- result$metadata_list[[1]]
```
