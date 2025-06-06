% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{dispersion_fun}
\alias{dispersion_fun}
\title{Calculate Dispersion for Microbiome Data}
\usage{
dispersion_fun(mean_abund, asymptDisp, extraPois)
}
\arguments{
\item{mean_abund}{Numeric value representing the mean abundance of the taxa.}

\item{asymptDisp}{Numeric value for the asymptotic dispersion (the dispersion at high abundance).}

\item{extraPois}{Numeric value for the extra Poisson dispersion (to model overdispersion).}
}
\value{
A numeric value representing the dispersion.
}
\description{
This function calculates the dispersion value for microbiome data based on the
provided parameters: mean abundance, asymptotic dispersion, and extra Poisson dispersion.
}
\details{
The dispersion is calculated using the formula:
\deqn{\text{dispersion} = \text{asymptDisp} + \frac{\text{extraPois}}{\text{mean_abund}}}
}
\examples{
mean_abund <- 10
asymptDisp <- 0.1
extraPois <- 0.05
dispersion_fun(mean_abund, asymptDisp, extraPois)

}
