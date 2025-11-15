# generate normal mixture parameters (prob vector, mean vector, sd vector for a specified set of 'x' values (logmean)

generate normal mixture parameters (prob vector, mean vector, sd vector
for a specified set of 'x' values (logmean)

## Usage

``` r
genmixpars(x, pars, np = 2, sd_ord = 2)
```

## Arguments

- x:

  independent variable

- pars:

  parameter vector: first logit-probs (np-1), then mean parameters (2
  per component: intercepts, slopes), then var parameters (varord + 1
  per component: intercepts, slopes, quad coeffs, etc.)

- np:

  number of components in mixture

- sd_ord:

  order of logsd model (2 = quadratic)

## Value

A list

     probs: mixture proportions ()

     muvals: mean values

     sdvals: standard deviation values
