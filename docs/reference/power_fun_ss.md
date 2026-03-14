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
\mathrm{logit}\\P(\text{reject})\\ = s(\text{logmean},
\text{abs\\lfc}) + s(\text{abs\\lfc}, \text{logsample\\size}) +
s(\text{logmean}, \text{logsample\\size}) \$\$

using a binomial family.

If the full model fails due to numerical instability, a simpler fallback
model is fitted instead: \$\$ \mathrm{logit}\\P(\text{reject})\\ =
s(\text{logmean}, \text{abs\\lfc}) + s(\text{logsample\\size}) \$\$

A warning is issued when the fallback model is used.

## See also

[`scam::scam()`](https://rdrr.io/pkg/scam/man/scam.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example structure only
pval_est_list <- list(
  c(0.01, 0.20, 0.03),
  c(0.15, 0.04, 0.08)
)

logmean_list <- list(
  c(2, 2, 2),
  c(3, 3, 3)
)

logfoldchange_list <- list(
  c(1, 1, 1),
  c(2, 2, 2)
)

nsample_vec <- c(20, 40, 80)

out <- power_fun_ss(
  pval_est_list = pval_est_list,
  logmean_list = logmean_list,
  nsample_vec = nsample_vec,
  logfoldchange_list = logfoldchange_list,
  alpha_level = 0.1
)

out$combined_data
out$gam_mod
} # }
```
