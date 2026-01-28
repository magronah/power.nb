# Fit a mixture of Gaussian distributions to log fold change

The standard deviation parameters are modeled either by linear or
quadratic functions of log mean count and the mean parameter is modeled
by linear functions of log mean count

## Usage

``` r
logfoldchange_fit(
  logmean,
  logfoldchange,
  ncore = 2,
  max_sd_ord = 2,
  max_np = 5,
  minval = -5,
  maxval = 5,
  itermax = 100,
  NP = 800,
  seed = 100
)
```

## Arguments

- logmean:

  vector of log mean abundance

- logfoldchange:

  vector of log fold change

- ncore:

  number of cores to use

- max_sd_ord:

  the maximum order of polynomial function to fit to standard deviation
  parameter. This must be either 1 (linear) or 2(quadratic)

- max_np:

  maximum number of Gaussian components to check for

- minval:

  minimum value for DEoptim search

- maxval:

  maximum value for DEoptim search

- itermax:

  maximum number of iterations

- NP:

  the number of population members for DEoptim

- seed:

  seed value

## Value

A list.

       par is a vector of the estimates of the mixture proportion,
       the mean and standard deviation parameters,

       np is the number of gaussian components fitted

       sd_ord is the order for the function for the standard deviation

       aic is the aic of the best fit

## Examples

``` r
logmean        =  rnorm(100)
logfoldchange  =  rnorm(100)
if (FALSE) { # \dontrun{
logfoldchange_fit(logmean,logfoldchange)
} # }
```
