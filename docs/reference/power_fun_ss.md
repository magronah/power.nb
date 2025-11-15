# Power function

Power function

## Usage

``` r
power_fun_ss(
  deseq_est_list,
  true_logfoldchange,
  true_logmean,
  sample_vec,
  alpha_level = 0.1,
  notu
)
```

## Arguments

- deseq_est_list:

  a list containing fold change, pvalues and other estimates from
  `DESeq2`

- true_logfoldchange:

  list containing simulated log fold change used for simulating the
  count data

- true_logmean:

  list containing simulated log mean counte used for simulating the
  count data

- sample_vec:

  vector of sample sizes

- alpha_level:

  sign containing simulated log fold change used for simulating the
  count data

- notu:

  number of OTUs

  fit GAM with covariates as tensor product (ie,interaction between log
  mean abundance and absolute log fold changes and then a spline for the
  sample sizes log mean abundance and log fold changes are related
  directly,hence the interaction but sample size is not quite related to
  log mean abundance and log fold changes directly
