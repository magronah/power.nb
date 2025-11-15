# Title

Title

## Usage

``` r
gam_fit(
  deseq_est_list,
  true_lfoldchange_list,
  true_lmean_list,
  grid_len = 50,
  alpha_level = 0.1
)
```

## Arguments

- deseq_est_list:

  a list containing fold change, pvalues and other estimates from
  `DESeq2`

- true_lfoldchange_list:

  list containing simulated log fold change used for simulating the
  count data

- true_lmean_list:

  list containing simulated log mean count used for simulating the count
  data

- grid_len:

  number of grids for

- alpha_level:

  significance level for power calculations

## Value

A list fit_2d is the fitted scam object

power_estimate predicted power estimated using fit_2d

combined_data tibble containing pvlaues and used for GAM fit
