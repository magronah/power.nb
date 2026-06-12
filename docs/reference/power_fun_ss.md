# Fit a smooth power model for sample size estimation

Fits a shape-constrained additive model (SCAM) to estimate statistical
power as a function of mean abundance, absolute log fold change, and
sample size. The function takes p-values obtained across simulation
settings, converts them to rejection indicators based on a chosen
significance level, and then models the probability of rejection using
smooth terms.

## Usage

``` r
power_fun_ss(
  pval_est_list,
  logmean_list,
  nsample_vec,
  logfoldchange_list,
  alpha_level = 0.1
)
```

## Arguments

- pval_est_list:

  A list of p-value vectors obtained from estimating the fold change
  estimates from the simulated datasets across different sample sizes.
  Each p-value corresponds to a test result for a particular combination
  of mean abundance, log fold change, and sample size.

- logmean_list:

  A list containing the simulated log mean abundance values
  corresponding to the p-values in `pval_est_list`.

- nsample_vec:

  A numeric vector of sample sizes used in the simulations.

- logfoldchange_list:

  A list containing the simulated log fold change values corresponding
  to the p-values in `pval_est_list`.

- alpha_level:

  Numeric significance level used to define rejection of the null
  hypothesis. Default is `0.1`.

## Value

A list with two components:

- `combined_data`: A tibble containing the combined simulation inputs
  and rejection indicators used for model fitting.

- `gam_mod`: The fitted `scam` model object.

## Details

The returned model can be used as a smooth approximation to the
empirical power surface, which is useful for interpolation and sample
size prediction.

The function first pools all p-values from `pval_est_list` and converts
them into binary rejection indicators: \$\$I(p \< \alpha)\$\$ where
`alpha_level` is the chosen significance threshold.

A data frame is then constructed with:

- `logmean`: log mean abundance,

- `abs_lfc`: absolute log fold change,

- `pval_reject`: binary rejection indicator,

- `sample_size`: sample size,

- `logsample_size`: base-2 logarithm of sample size.

The function attempts to fit the following SCAM model: \$\$
\mathrm{logit}\\P(\text{reject})\\ = s(\text{logmean}, \text{abs\\lfc},
\text{bs = "tedmi"}) + s(\text{abs\\lfc}, \text{logsample\\size},
\text{bs = "tedmi"}) + s(\text{logmean}, \text{logsample\\size},
\text{bs = "tedmi"}) \$\$

using a binomial family.

If the full model fails due to numerical instability, a simpler fallback
model is fitted instead: \$\$ \mathrm{logit}\\P(\text{reject})\\ =
s(\text{logmean}, \text{abs\\lfc}, \text{bs = "tedmi"}) +
s(\text{logsample\\size}, \text{bs = "mpi"}) \$\$

A warning is issued when the fallback model is used.

## See also

[`scam::scam()`](https://rdrr.io/pkg/scam/man/scam.html)

## Examples

``` r
#' @examples

# Example structure only
set.seed(101)
n = 70
pval_est_list = list(rnorm(n),rnorm(n))
logmean_list = list(rnorm(n),rnorm(n))
logfoldchange_list = list(rnorm(n),rnorm(n))
nsample_vec <- c(20, 40)
out <- power_fun_ss(
pval_est_list = pval_est_list,
logmean_list = logmean_list,
nsample_vec = nsample_vec,
logfoldchange_list = logfoldchange_list,
alpha_level = 0.1
)
#> Warning: Full SCAM model failed to converge due to numerical instability.
#> The simpler model will be used instead.
#> Original error: Model has more coefficients than data
#> Selected model: simpler model (full model did not converge).
out$combined_data
#> # A tibble: 140 × 5
#>    logmean abs_lfc pval_reject sample_size logsample_size
#>      <dbl>   <dbl>       <dbl>       <dbl>          <dbl>
#>  1  -0.909   1.17            1          20           4.32
#>  2  -0.338   2.15            0          20           4.32
#>  3  -1.41    0.342           1          20           4.32
#>  4   0.218   0.905           0          20           4.32
#>  5   0.670   1.10            0          20           4.32
#>  6  -0.288   1.47            0          20           4.32
#>  7   0.469   0.281           0          20           4.32
#>  8  -0.470   0.846           1          20           4.32
#>  9  -0.239   1.29            0          20           4.32
#> 10  -0.447   0.312           1          20           4.32
#> # ℹ 130 more rows
out$gam_mod
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> pval_reject ~ s(logmean, abs_lfc, bs = "tedmi") + s(logsample_size, 
#>     bs = "mpi")
#> <environment: 0x55959ad70c00>
#> 
#> Estimated degrees of freedom:
#> 1 1  total = 3 
#> 
#> UBRE score: 0.41649
#> rank: 57/58

```
