# Density of a Normal Mixture Model

This function calculates the density of a normal mixture model for a
given vector of parameters on the unconstrained scale (softmax(prob),
mean, log(sd)).

## Usage

``` r
dnormmix(x, par, logmean, ..., log = FALSE)
```

## Arguments

- x:

  Numeric vector of values at which to evaluate the density.

- par:

  A vector of parameters on the unconstrained scale, including:

  - `softmax(prob)`: Mixture probabilities (on the unconstrained scale,
    transformed via softmax).

  - `mean`: Means of the normal components.

  - `log(sd)`: Logarithms of standard deviations of the normal
    components.

- logmean:

  Numeric value representing the log of the mean parameter.

- ...:

  Additional arguments passed to the `genmixpars` function. Defaults:
  two components (`np = 2`) and a quadratic model for standard deviation
  parameters (`sd_ord = 2`).

- log:

  Logical. If `TRUE`, the logarithm of the density is returned. Default
  is `FALSE`.

## Value

A numeric vector of density values (or log-density values if
`log = TRUE`) for the mixture model.

## Examples

``` r
# Example parameters
x <- seq(-3, 3, length.out = 100)
## par <- c(-0.5, 0.5, log(0.8), log(1.2))  # Example: softmax probabilities, mean, log(sd)
set.seed(101); par <- rnorm(11)
logmean <- rep(0.1, length(x))  ## constant log mean

# Calculate density
density <- dnormmix(x, par, logmean)

# Calculate log-density
log_density <- dnormmix(x, par, logmean, log = TRUE)
```
