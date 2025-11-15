# Generate Parameter Names for Mixture Model

This function generates parameter names for a Gaussian mixture model
based on the number of components (`np`) and the order of the polynomial
function (`sd_ord`) used to model the standard deviation parameters.

## Usage

``` r
gen_parnames(np, sd_ord)
```

## Arguments

- np:

  Integer. The number of Gaussian components in the mixture model.

- sd_ord:

  Integer. The order of the polynomial function used to model the
  standard deviation parameters. Possible values are:

  - `1`: Linear function.

  - `2`: Quadratic function.

## Value

A character vector of parameter names, including:

- Logit-transformed probabilities (`logitprob_1`, ...,
  `logitprob_(np-1)`).

- Mean parameters (`mu_int_1`, `mu_slope_1`, ..., for each component).

- Log-transformed standard deviations (`logsd_.1_1`, `logsd_.L_1`, ...,
  depending on `sd_ord` and `np`).

## Examples

``` r
# Generate parameter names for a 3-component mixture with linear standard deviation function
gen_parnames(np = 3, sd_ord = 1)
#>  [1] "logitprob_1" "logitprob_2" "mu_int_1"    "mu_int_2"    "mu_int_3"   
#>  [6] "mu_slope_1"  "mu_slope_2"  "mu_slope_3"  "logsd_.1_1"  "logsd_.1_2" 
#> [11] "logsd_.1_3"  "logsd_.L_1"  "logsd_.L_2"  "logsd_.L_3" 

# Generate parameter names for a 4-component mixture with quadratic standard deviation function
gen_parnames(np = 4, sd_ord = 2)
#>  [1] "logitprob_1" "logitprob_2" "logitprob_3" "mu_int_1"    "mu_int_2"   
#>  [6] "mu_int_3"    "mu_int_4"    "mu_slope_1"  "mu_slope_2"  "mu_slope_3" 
#> [11] "mu_slope_4"  "logsd_.1_1"  "logsd_.1_2"  "logsd_.1_3"  "logsd_.1_4" 
#> [16] "logsd_.L_1"  "logsd_.L_2"  "logsd_.L_3"  "logsd_.L_4"  "logsd_.Q_1" 
#> [21] "logsd_.Q_2"  "logsd_.Q_3"  "logsd_.Q_4" 
```
