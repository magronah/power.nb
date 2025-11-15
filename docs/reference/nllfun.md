# Objective function

Objective function

## Usage

``` r
nllfun(par, vals, logmean, np, sd_ord)
```

## Arguments

- par:

  parameters of the mixture of Gaussian

- vals:

  values of log fold change

- logmean:

  log mean count

- np:

  number of Gaussian components

- sd_ord:

  order of polynomial function to model standard deviation parameters
  (1 - linear function and 2- quad)

## Value

value for the objective function
