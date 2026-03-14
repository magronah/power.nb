#
# #######################################################################
# inverse_fun <- function(target,lmb,abs_lfc, model,xmin, xmax) {
#
#   if (length(target) != length(lmb) || length(target) != length(abs_lfc)) {
#     stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
#   }
#
#   # Initialize an empty list to store the roots
#   roots <- list()
#
#   # Loop through each element and solve for the root
#   for (i in 1:length(target)) {
#     roots[[i]] <- uniroot(function(s) {
#       predict(model,
#               type = "response",
#               newdata = data.frame(sample_size = s,
#                                    lmean_abund = lmb[i],
#                                    abs_lfc = abs_lfc[i])
#       ) - target[i]
#     },
#     interval = c(xmin, xmax),
#     extendInt = "yes")$root
#   }
#
#   unlist(roots)
# }
#
#
# uniroot_lmb =  function(target_power,sample_size,abs_lfc,model,xmin,xmax){
#
#   root <- uniroot(function(lmb) {
#     predict(model,
#             type = "response",
#             newdata = data.frame(sample_size = sample_size,
#                                  lmean_abund = lmb,
#                                  abs_lfc = abs_lfc)
#     ) - target_power
#   },
#   interval = c(xmin, xmax),
#   extendInt = "yes")$root
#   root
# }
#
#
#
# inverse_lfc <- function(target,lmb,samp_size, model,xmin, xmax) {
#
#   if (length(target) != length(lmb) || length(target) != length(samp_size)) {
#     stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
#   }
#
#   # Initialize an empty list to store the roots
#   roots <- list()
#
#   # Loop through each element and solve for the root
#   for (i in 1:length(target)) {
#     roots[[i]] <- uniroot(function(s) {
#       predict(model,
#               type = "response",
#               newdata = data.frame(sample_size =  samp_size[i],
#                                    lmean_abund =  lmb[i],
#                                    abs_lfc     =  s)
#       ) - target[i]
#     },
#     interval = c(xmin, xmax),
#     extendInt = "yes")$root
#   }
#
#   unlist(roots)
# }
#
#
#
# inverse_lmb <- function(target,abs_lfc,samp_size, model,xmin, xmax) {
#
#   if (length(target) != length(abs_lfc) || length(target) != length(samp_size)) {
#     stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
#   }
#
#   # Initialize an empty list to store the roots
#   roots <- list()
#
#   # Loop through each element and solve for the root
#   for (i in 1:length(target)) {
#     roots[[i]] <- uniroot(function(lmb) {
#       predict(model,
#               type = "response",
#               newdata = data.frame(sample_size =  samp_size[i],
#                                    lmean_abund =  lmb,
#                                    abs_lfc     =  abs_lfc[i])
#       ) - target[i]
#     },
#     interval = c(xmin, xmax),
#     extendInt = "yes")$root
#   }
#
#   unlist(roots)
# }
#
#
#
#
#



# # Function to compute the missing input
# power.nb <- function(power, sample_size, logfoldchange, logmean_abundance, gam_mod) {
#
#   if (is.numeric(power) && is.numeric(sample_size) && is.numeric(logfoldchange) && is.numeric(logmean_abundance)) {
#     if (!missing(power) && !missing(sample_size) && !missing(logfoldchange)) {
#       if (missing(logmean_abundance)) {
#         # Compute the missing input
#
#       }
#     } else if (!missing(power) && !missing(sample_size) && !missing(logmean_abundance)) {
#       # Compute the missing input
#       logfoldchange <- logmean_abundance - power - sample_size
#     } else if (!missing(power) && !missing(logfoldchange) && !missing(logmean_abundance)) {
#       # Compute the missing input
#       sample_size <- logmean_abundance - power - logfoldchange
#     } else if (!missing(sample_size) && !missing(logfoldchange) && !missing(logmean_abundance)) {
#       # Compute the missing input
#       power <- logmean_abundance - sample_size - logfoldchange
#     } else {
#       stop("Please specify at least 3 inputs.")
#     }
#
#     # Return the computed inputs
#     return(list(power = power, sample_size = sample_size, logfoldchange = logfoldchange, logmean_abundance = logmean_abundance))
#   } else {
#     stop("Inputs must be numeric values.")
#   }
# }



#' Fit a smooth power model for sample size estimation
#'
#' Fits a shape-constrained additive model (SCAM) to estimate statistical power
#' as a function of mean abundance, absolute log fold change, and sample size.
#' The function takes p-values obtained across simulation settings, converts them
#' to rejection indicators based on a chosen significance level, and then models
#' the probability of rejection using smooth terms.
#'
#' The returned model can be used as a smooth approximation to the empirical power
#' surface, which is useful for interpolation and sample size prediction.
#'
#' @param pval_est_list A list of p-value vectors obtained from
#'   estimating the fold change estimates from the simulated datasets
#'   across different sample sizes.
#'   Each p-value corresponds to a test result for a particular combination of
#'   mean abundance, log fold change, and sample size.
#' @param logmean_list A list containing the simulated log mean abundance values
#' corresponding to the p-values in `pval_est_list`.
#' @param nsample_vec A numeric vector of sample sizes used in the simulations.
#' @param logfoldchange_list A list containing the simulated log fold change values
#'   corresponding to the p-values in `pval_est_list`.
#' @param alpha_level Numeric significance level used to define rejection of the
#'   null hypothesis. Default is `0.1`.
#'
#' @details
#' The function first pools all p-values from `pval_est_list` and converts them
#' into binary rejection indicators:
#' \deqn{I(p < \alpha)}
#' where `alpha_level` is the chosen significance threshold.
#'
#' A data frame is then constructed with:
#' \itemize{
#'   \item `logmean`: log mean abundance,
#'   \item `abs_lfc`: absolute log fold change,
#'   \item `pval_reject`: binary rejection indicator,
#'   \item `sample_size`: sample size,
#'   \item `logsample_size`: base-2 logarithm of sample size.
#' }
#'
#' The function attempts to fit the following SCAM model:
#' \deqn{
#' \mathrm{logit}\{P(\text{reject})\} =
#' s(\text{logmean}, \text{abs\_lfc}) +
#' s(\text{abs\_lfc}, \text{logsample\_size}) +
#' s(\text{logmean}, \text{logsample\_size})
#' }
#'
#' using a binomial family.
#'
#' If the full model fails due to numerical instability, a simpler fallback model
#' is fitted instead:
#' \deqn{
#' \mathrm{logit}\{P(\text{reject})\} =
#' s(\text{logmean}, \text{abs\_lfc}) +
#' s(\text{logsample\_size})
#' }
#'
#' A warning is issued when the fallback model is used.
#'
#' @return A list with two components:
#' \itemize{
#'   \item `combined_data`: A tibble containing the combined simulation inputs
#'   and rejection indicators used for model fitting.
#'   \item `gam_mod`: The fitted `scam` model object.
#' }
#'
#' @examples
#' \dontrun{
#' # Example structure only
#' pval_est_list <- list(
#'   c(0.01, 0.20, 0.03),
#'   c(0.15, 0.04, 0.08)
#' )
#'
#' logmean_list <- list(
#'   c(2, 2, 2),
#'   c(3, 3, 3)
#' )
#'
#' logfoldchange_list <- list(
#'   c(1, 1, 1),
#'   c(2, 2, 2)
#' )
#'
#' nsample_vec <- c(20, 40, 80)
#'
#' out <- power_fun_ss(
#'   pval_est_list = pval_est_list,
#'   logmean_list = logmean_list,
#'   nsample_vec = nsample_vec,
#'   logfoldchange_list = logfoldchange_list,
#'   alpha_level = 0.1
#' )
#'
#' out$combined_data
#' out$gam_mod
#' }
#'
#' @seealso [scam::scam()]
#' @export


power_fun_ss <- function(pval_est_list,
                         logmean_list,
                         nsample_vec,
                         logfoldchange_list,
                         alpha_level=0.1){

  # concatenate all p-values from all the sample size
        p_val    =   unname(unlist(pval_est_list))
  pval_reject    =   (!is.na(p_val) & p_val < alpha_level)
  sample_size    =   rep(nsample_vec,
                         times = sapply(lapply(logfoldchange_list, unlist),
                                        length))

  # create a table with all the information
  comb   = tibble::tibble(logmean      =   unname(unlist(logmean_list)),
                          abs_lfc      =   abs(unname(unlist(logfoldchange_list))),
                          pval_reject  =   as.numeric(pval_reject),
                          sample_size  =   sample_size)

  comb$logsample_size = log2(comb$sample_size)

  fit_2d <- tryCatch(
    {
      scam::scam(
        pval_reject ~
          s(logmean, abs_lfc, bs = "tedmi") +
          s(abs_lfc, logsample_size, bs = "tedmi") +
          s(logmean, logsample_size, bs = "tedmi"),
        data = comb,
        family = binomial
      )
    },
    error = function(e) {
      warning(
        "Full SCAM model failed to converge due to numerical instability. ",
        "A simpler fallback model was fitted instead (removing the interactions ",
        "between sample size and log mean abundance and between sample size and ",
        "effect size).\n",
        "To attempt the full model, consider refitting with a reduced set of sample sizes \n",
        "Original error: ", e$message,
        call. = FALSE)
      scam::scam(
        pval_reject ~
          s(logmean, abs_lfc, bs = "tedmi") +
          s(logsample_size, bs = "mpi"),
        data = comb,
        family = binomial
      )
    }
  )
  list(combined_data=comb, gam_mod = fit_2d)
}


###############################################################
#' Sample Size estimation function uisng uniroot
#'
#' @param target_power Numeric value specifying the desired statistical power.
#' @param logmean Numeric value representing the log of the mean abundance.
#' @param abs_lfc Numeric value representing the absolute log fold change.
#' @param model A fitted GAM/SCAM model used to predict statistical power.
#' @param xmin Numeric value giving the minimum sample size considered in the search.
#' @param xmax Numeric value giving the maximum sample size considered in the search.
#' @return A numeric value corresponding to the estimated sample size required to achieve the target power.
#' @export
uniroot_ss =  function(target_power,logmean, abs_lfc,model,xmin,xmax){

  root <- stats::uniroot(function(ss) {

    data = data.frame(logsample_size = ss,
                      logmean = logmean,
                      abs_lfc = abs_lfc)

    ##' predict.scam seem to have a bug and
    ##' would not predict row data with only one row
    pred = predict(model,
                   type = "response",
                   newdata = data[rep(1,2),])

    pred[[1]] - target_power
  },
  interval = c(xmin, xmax),
  extendInt = "no")$root
  2^root
}
###############################################################
#' Estimate sample size required to achieve a target statistical power
#'
#' This function estimates the sample size required to achieve a specified
#' target power using predictions from a fitted GAM/SCAM power model.
#' The function evaluates predicted power across a grid of candidate sample
#' sizes and identifies the smallest sample size for which the predicted
#' power reaches or exceeds the target value. Linear interpolation is then
#' used on the log2(sample size) scale to obtain a more precise estimate.
#'
#' @param target_power Numeric value specifying the desired statistical power.
#' @param logmean Numeric value representing the log mean abundance.
#' @param abs_lfc Numeric value representing the absolute log fold change.
#' @param model A fitted GAM/SCAM model used to predict statistical power.
#' @param xmin Numeric value giving the minimum log2(sample size) considered
#'   in the search. Default is \code{log2(5)}.
#' @param xmax Numeric value giving the maximum log2(sample size) considered
#'   in the search. Default is \code{log2(500)}.
#' @param ngrid Integer specifying the number of grid points used when
#'   searching for the sample size solution. Default is \code{1000}.
#'
#' @return A numeric value representing the estimated sample size required
#'   to achieve the target power. Returns \code{NA} if the target power is
#'   not reached within the specified search range.
#'
#' @details
#' The function first constructs a grid of candidate sample sizes on the
#' log2 scale between \code{xmin} and \code{xmax}. Predicted power values are
#' then obtained from the fitted model for each grid point. The smallest
#' sample size at which the predicted power reaches the target value is
#' identified, and linear interpolation is used to refine the estimate.
#'
#' @export
#'
sample_size_ss_interp <- function(target_power, logmean, abs_lfc, model,
                                  xmin = log2(5), xmax = log2(500),
                                  ngrid = 1000) {

  ss_grid <- seq(xmin, xmax, length.out = ngrid)

  nd <- data.frame(
    logsample_size = ss_grid,
    logmean = logmean,
    abs_lfc = abs_lfc
  )

  pred <- predict(model, type = "response", newdata = nd)

  idx <- which(pred >= target_power)[1]

  if (is.na(idx)) return(NA_real_)     # target never reached
  if (idx == 1) return(2^ss_grid[1])   # already achieved at minimum

  x0 <- ss_grid[idx - 1]
  x1 <- ss_grid[idx]
  y0 <- pred[idx - 1]
  y1 <- pred[idx]

  # linear interpolation on log2(sample size) scale
  x_star <- x0 + (target_power - y0) * (x1 - x0) / (y1 - y0)

  2^x_star
}
##########################################################
#' Solve for the sample size required to achieve a target statistical power
#'
#' This function estimates the sample size required to achieve a specified
#' target statistical power using a fitted GAM/SCAM power model. The function
#' first attempts to solve for the sample size using a root-finding algorithm.
#' If the root-finding procedure fails, a grid-based interpolation method is
#' used as a fallback to obtain an approximate solution.
#'
#' @param target_power Numeric value specifying the desired statistical power.
#' @param logmean Numeric value representing the log mean abundance.
#' @param abs_lfc Numeric value representing the absolute log fold change.
#' @param model A fitted GAM/SCAM model used to predict statistical power.
#' @param xmin Numeric value giving the minimum log2(sample size) considered
#'   in the search. Default is \code{log2(5)}.
#' @param xmax Numeric value giving the maximum log2(sample size) considered
#'   in the search. Default is \code{log2(500)}.
#'
#' @return A numeric value representing the estimated sample size required
#'   to achieve the target power.
#'
#' @details
#' The function first attempts to compute the required sample size using a
#' root-finding approach implemented in \code{uniroot_ss()}. If this method
#' fails (for example, due to numerical issues or if the root cannot be
#' bracketed within the specified interval), the function falls back to
#' a grid-based interpolation approach implemented in
#' \code{sample_size_ss_interp()}.
#'
#' A warning is issued if the target power is specified as 0 or 1, since
#' these values correspond to unrealistic design targets in most practical
#' applications.
#'
#' @export
ss_solver <- function(target_power, logmean, abs_lfc, model,
                      xmin = log2(5), xmax = log2(500)) {

  if(target_power == 0 || target_power == 1){
    warning("statistical power for  0% or 100% is  unlikely")
  }
  out <- tryCatch({

    uniroot_ss(
      target_power = target_power,
      logmean = logmean,
      abs_lfc = abs_lfc,
      model = model,
      xmin = xmin,
      xmax = xmax
    )

  }, error = function(e) {

    sample_size_ss_interp(
      target_power = target_power,
      logmean = logmean,
      abs_lfc = abs_lfc,
      model = model,
      xmin = xmin,
      xmax = xmax
    )

  })

  out
}
