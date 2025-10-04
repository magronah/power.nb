#######################################################
#' Simulate Log Means for OTUs
#'
#' This function generates log means for a specified number of OTUs (Operational Taxonomic Units) based on the provided parameters.
#' If a single mean is specified, the log means are drawn from a normal distribution.
#' If multiple means and corresponding weights are specified, the log means are drawn from a mixture of normal distributions.
#'
#' @param logmean_param A list containing the parameters for the distribution:
#' \itemize{
#'   \item \code{mu}: A single value or a vector of mean(s) for the normal or mixture distribution.
#'   \item \code{sigma}: The standard deviation(s) for the normal or mixture distribution.
#'   \item \code{lambda}: (Optional) A vector of weights for the components of the mixture distribution. Required only if \code{mu} has more than one value.
#' }
#' @param notu An integer specifying the number of OTUs to simulate.
#'
#' @return A numeric vector of simulated log means for the specified number of OTUs.
#'
#' @export
#'
#' @examples
#' # Example 1: Single normal distribution
#' params_single <- list(mu = 0, sigma = 1)
#' logmean_sim_fun(logmean_param = params_single, notu = 100)
#'
#' # Example 2: Mixture of normal distributions
#' params_mixture <- list(
#'   mu = c(-1, 1),
#'   sigma = c(0.5, 0.5),
#'   lambda = c(0.4, 0.6)
#' )
#' logmean_sim_fun(logmean_param = params_mixture, notu = 100)
logmean_sim_fun = function(logmean_param,notu){

  l=length(logmean_param$mu)
  if(l == 1){
    logmean = rnorm(n = notu,
                    mean = logmean_param$mu,
                    sd   = logmean_param$sigma)
  }else{
    logmean = mixtools::rnormmix(n=notu,
                                 lambda = logmean_param$lambda,
                                 mu = logmean_param$mu,
                                 sigma = logmean_param$sigma)
  }
  logmean
}

##############################################################
#' Simulate Log Fold Change Values
#'
#' This function generates simulated log fold change (LFC) values based on the provided log mean abundance and LFC parameters.
#' The simulation ensures that the generated LFC values remain within a specified maximum range by iterating until convergence or until a maximum iteration limit is reached.
#'
#' @param logmean_sim A numeric vector of simulated log mean abundances.
#' @param logfoldchange_param A list containing parameters for the log fold change simulation:
#' \itemize{
#'   \item \code{par}: Optimal parameters for the log fold change fit.
#'   \item \code{np}: Optimal number of components for the log fold change model.
#'   \item \code{sd_ord}: Order of the polynomial used for the standard deviation parameter of the log fold change.
#' }
#'
#' @param max_lfc A numeric value specifying the maximum allowable absolute log fold change value. Default is 15.
#' @param max_iter An integer specifying the maximum number of iterations allowed to ensure all simulated LFC values are within the \code{max_lfc} range. Default is 10,000.
#'
#' @return A numeric vector of simulated log fold change values (\code{lfc}).
#'
#' @export
#'
#' @examples
#' # Define simulated log mean abundance
#' logmean_sim <- rnorm(100, mean = 0, sd = 1)
#'
#' # Define parameters for log fold change simulation
#' logfoldchange_param <- list(
#'   par = c(1, -0.5, 0.2), # Example parameters
#'   np = 2,                # Number of components
#'   sd_ord = 2             # Order of polynomial for SD
#' )
#'
#' # Simulate log fold change values
#' logfoldchange_sim_fun(
#'   logmean_sim = logmean_sim,
#'   logfoldchange_param = logfoldchange_param,
#'   max_lfc = 10,
#'   max_iter = 5000
#' )
logfoldchange_sim_fun <- function(logmean_sim, logfoldchange_param,
                                  max_lfc = 15, max_iter = 10000) {

  par   =   logfoldchange_param$par
  np    =   logfoldchange_param$np
  sd_ord  =   logfoldchange_param$sd_ord

  lfc <- myrnormmix(par, logmean_sim, np = np, sd_ord = sd_ord)
  r <- range(lfc);   iteration_count <- 0

  while (max(abs(r)) > max_lfc && iteration_count < max_iter) {
    lfc <- myrnormmix(par, logmean_sim, np = np, sd_ord = sd_ord)
    r <- range(lfc)
    iteration_count <- iteration_count + 1
  }

  if (iteration_count == max_iter) {
    warning("Maximum number of iterations reached without convergence.")
  }

  return(lfc)
}

##############################################################
#' Simulate Count Data for Microbiome Studies
#'
#' This function simulates count data for microbiome studies based on log mean, log fold change,
#' and dispersion parameters. It supports generating data for multiple simulations and allows
#' flexibility in specifying the number of control and treatment samples or samples per group.
#'
#' @param logmean_param A list of parameters for simulating the log mean abundance.
#' @param logfoldchange_param A list of parameters for simulating log fold change, containing:
#'   \itemize{
#'     \item \code{par}: Optimal parameters for log fold change fitting.
#'     \item \code{np}: Number of components for the log fold change model.
#'     \item \code{sd_ord}: Order of the polynomial for the standard deviation parameter.
#'   }
#' @param dispersion_param A list of dispersion parameters containing:
#'   \itemize{
#'     \item \code{asymptDisp}: Asymptotic dispersion parameter.
#'     \item \code{extraPois}: Additional Poisson variation parameter.
#'   }
#' @param nsamp_per_group Number of samples per group (control and treatment). If provided,
#'   \code{ncont} and \code{ntreat} must not be specified.
#' @param ncont Number of control samples. Specify along with \code{ntreat} when \code{nsamp_per_group} is not provided.
#' @param ntreat Number of treatment samples. Specify along with \code{ncont} when \code{nsamp_per_group} is not provided.
#' @param notu Number of operational taxonomic units (OTUs) to simulate.
#' @param nsim Number of simulations to run. Default is 1.
#' @param disp_scale Scale parameter for the dispersion. Default is 0.3.
#' @param max_lfc Maximum allowable log fold change. Default is 15.
#' @param maxlfc_iter Maximum number of iterations for ensuring log fold change is within \code{max_lfc}. Default is 1,000.
#' @param seed Seed value for reproducibility. Default is \code{NULL}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{countdata_list}: A list of count data matrices for each simulation.
#'     \item \code{metadata_list}: A list of metadata data frames for each simulation.
#'     \item \code{logmean_list}: A list of log mean vectors for each simulation.
#'     \item \code{logfoldchange_list}: A list of log fold change vectors for each simulation.
#'     \item \code{treat_countdata_list}: A list of treatment count data matrices for each simulation.
#'     \item \code{control_countdata_list}: A list of control count data matrices for each simulation.
#'   }
#'
#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach %do% %dopar% foreach
#' @examples
#' # Load required packages
#' library(foreach)
#' library(doParallel)
#' # Define parameters
#' logmean_param <- list(mu = 0, sigma = 1)
#' logfoldchange_param <- list(par = rnorm(11), np = 2, sd_ord = 2)
#' dispersion_param <- list(asymptDisp = 0.1, extraPois = 0.05)
#'
#' # Simulate count data
#' result <- countdata_sim_fun(
#'   logmean_param = logmean_param,
#'   logfoldchange_param = logfoldchange_param,
#'   dispersion_param = dispersion_param,
#'   nsamp_per_group = 10,
#'   notu = 50,
#'   nsim = 2,
#'   seed = 123
#' )
#'
#' # Access simulation results
#' countdata <- result$countdata_list[[1]]
#' metadata <- result$metadata_list[[1]]

countdata_sim_fun <- function(logmean_param, logfoldchange_param, dispersion_param,
                              nsamp_per_group = NULL, ncont = NULL,ntreat = NULL,
                              notu, nsim = 1,   disp_scale=0.3, max_lfc = 15,
                              maxlfc_iter = 1000, seed = NULL){

  # Check if both nsamp_per_group and (ncont or ntreat) are provided
  if (!is.null(nsamp_per_group) && (!is.null(ncont) || !is.null(ntreat))) {
    stop("Please specify either 'nsamp_per_group' or 'ncont' and 'ntreat', but not both.")
  }

  if (is.null(nsamp_per_group)) {
    if (is.null(ncont) || is.null(ntreat)) {
      stop("Please specify both 'ncont' and 'ntreat' when 'nsamp_per_group' is not provided.")
    }
  }

  par     =   logfoldchange_param$par
  np      =   logfoldchange_param$np
  sd_ord  =   logfoldchange_param$sd_ord

  set.seed(seed)

  res = foreach(j = 1:nsim,
                .export = c("logfoldchange_sim_fun","logmean_sim_fun",
                            "myrnormmix", "rnormmix0","genmixpars",
                            "dispersion_fun"),
                .packages = c("stats","purrr", "mixtools"),.verbose = F) %do% {

                  ####Simulate log mean and log fold change
                  logmean        =   logmean_sim_fun(logmean_param,notu)
                  logfoldchange  =   logfoldchange_sim_fun(logmean,logfoldchange_param,
                                                           max_lfc = max_lfc,
                                                           max_iter = maxlfc_iter)

                  ####Calculate control and treatment mean abundances
                  mean_abund =  2^logmean; foldchange = 2^logfoldchange
                  control    =  (2*mean_abund)/(1+ foldchange)

                  ####predict dispersion
                  dispersion  =   mean_abund |>
                    purrr::map(function(x) dispersion_fun(x, dispersion_param$asymptDisp,
                                                   dispersion_param$extraPois)) |>  unlist()

                  dd = data.frame(control = control,
                                  treatment = control*foldchange,
                                  dispersion = dispersion)

                  ####Simulate count data
                  if(!is.null(nsamp_per_group)){
                    countdata = data.frame(apply(dd, 1, function(x) {
                      rnbinom(n=2*nsamp_per_group,
                              mu = rep(c(x["control"],x["treatment"]),
                                       each= nsamp_per_group),
                              size = 1/(disp_scale*x["dispersion"]))
                    }))

                    metadata =  data.frame(Samples=paste0("sample_",1:(2*nsamp_per_group)),
                                           Groups=factor(rep(c("control", "treatment"),
                                                             each=(nsamp_per_group))))

                    rownames(countdata) =   paste0("sample_",1:(2*nsamp_per_group))
                    colnames(countdata) =   paste0("otu_",1:notu)

                    control_count       =   countdata[1:nsamp_per_group, ]
                    treat_count         =   countdata[(nsamp_per_group+1):(2*nsamp_per_group), ]

                  }else{
                    n  =  ncont + ntreat
                    countdata = data.frame(apply(dd, 1, function(x) {
                      stats::rnbinom(n = n,
                              mu = rep(c(x["control"],x["treatment"]),
                                       c(ncont,ntreat)),
                              size = 1/(disp_scale*x["dispersion"]))
                    }))

                    metadata =  data.frame(Samples=paste0("sample_",1:n),
                                           Groups=factor(rep(c("control", "treatment"),
                                                             c(ncont,ntreat))))

                    rownames(countdata) = paste0("sample_",1:n)
                    colnames(countdata) = paste0("otu_",1:notu)

                    control_count = countdata[1:ncont,]
                    treat_count   = countdata[(ncont+1):n,]
                  }
                  names(logmean)  =   paste0("otu_",1:notu)
                  names(logfoldchange)  = paste0("otu_",1:notu)


                  list(countdata     =   t(countdata),
                       control_count =   t(control_count),
                       treat_count   =   t(treat_count),
                       metadata      =   metadata,
                       logmean       =   logmean,
                       logfoldchange =   logfoldchange)
                }


  names(res)      =   paste0("sim_",1:nsim)

  ## extract individual results into a list
  list(  countdata_list  =   read_data(res,"countdata"),
         metadata_list   =   read_data(res,"metadata"),
         logmean_list    =   read_data(res,"logmean"),
         logfoldchange_list     =  read_data(res,"logfoldchange"),
         treat_countdata_list   =  read_data(res,"treat_count"),
         control_countdata_list =  read_data(res,"control_count")
  )

}

#######################################################
#' Calculate Dispersion for Microbiome Data
#'
#' This function calculates the dispersion value for microbiome data based on the
#' provided parameters: mean abundance, asymptotic dispersion, and extra Poisson dispersion.
#'
#' The dispersion is calculated using the formula:
#' \deqn{\text{dispersion} = \text{asymptDisp} + \frac{\text{extraPois}}{\text{mean_abund}}}
#'
#' @param mean_abund Numeric value representing the mean abundance of the taxa.
#' @param asymptDisp Numeric value for the asymptotic dispersion (the dispersion at high abundance).
#' @param extraPois Numeric value for the extra Poisson dispersion (to model overdispersion).
#'
#' @return A numeric value representing the dispersion.
#'
#' @examples
#' mean_abund <- 10
#' asymptDisp <- 0.1
#' extraPois <- 0.05
#' dispersion_fun(mean_abund, asymptDisp, extraPois)
#'
#' @export
dispersion_fun <- function(mean_abund, asymptDisp, extraPois) {
  asymptDisp + extraPois / mean_abund
}


