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
#' s(\text{logmean}, \text{abs\_lfc}, \text{bs = "tedmi"}) +
#' s(\text{abs\_lfc}, \text{logsample\_size}, \text{bs = "tedmi"}) +
#' s(\text{logmean}, \text{logsample\_size}, \text{bs = "tedmi"})
#' }
#'
#' using a binomial family.
#'
#' If the full model fails due to numerical instability, a simpler fallback model
#' is fitted instead:
#' \deqn{
#' \mathrm{logit}\{P(\text{reject})\} =
#' s(\text{logmean}, \text{abs\_lfc}, \text{bs = "tedmi"}) +
#' s(\text{logsample\_size}, \text{bs = "mpi"})
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
#' @export
#'
#' @examples
#' #' @examples
#'
#' # Example structure only
#' set.seed(101)
#' n = 70
#' pval_est_list = list(rnorm(n),rnorm(n))
#' logmean_list = list(rnorm(n),rnorm(n))
#' logfoldchange_list = list(rnorm(n),rnorm(n))
#' nsample_vec <- c(20, 40)
#' out <- power_fun_ss(
#' pval_est_list = pval_est_list,
#' logmean_list = logmean_list,
#' nsample_vec = nsample_vec,
#' logfoldchange_list = logfoldchange_list,
#' alpha_level = 0.1
#' )
#' out$combined_data
#' out$gam_mod
#'
#'
#' @seealso [scam::scam()]

power_fun_ss <- function(pval_est_list,
                         logmean_list,
                         nsample_vec,
                         logfoldchange_list,
                         alpha_level = 0.1) {

  p_val <- unname(unlist(pval_est_list))
  pval_reject <- (!is.na(p_val) & p_val < alpha_level)

  sample_size <- rep(
    nsample_vec,
    times = sapply(lapply(logfoldchange_list, unlist), length)
  )

  comb <- tibble::tibble(
    logmean = unname(unlist(logmean_list)),
    abs_lfc = abs(unname(unlist(logfoldchange_list))),
    pval_reject = as.numeric(pval_reject),
    sample_size = sample_size
  )

  comb$logsample_size <- log2(comb$sample_size)

  # Fit simpler model (always)
  simple_mod <- tryCatch(
    {
      scam::scam(
        pval_reject ~
          s(logmean, abs_lfc, bs = "tedmi") +
          s(logsample_size, bs = "mpi"),
        data = comb,
        family = binomial
      )
    },
    error = function(e) {
      stop(
        "The simpler SCAM model failed to converge.\n",
        "Original error: ", e$message,
        call. = FALSE
      )
    }
  )

  # Try full model
  full_mod <- tryCatch(
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
        "Full SCAM model failed to converge due to numerical instability.\n",
        "The simpler model will be used instead.\n",
        "Original error: ", e$message,
        call. = FALSE
      )
      NULL
    }
  )

  # Model selection
  if (is.null(full_mod)) {

    aic_table <- data.frame(
      model = "simple",
      AIC = stats::AIC(simple_mod)
    )

    selected_model <- "simple"
    best_mod <- simple_mod

    message("Selected model: simpler model (full model did not converge).")

  } else {

    aic_simple <- stats::AIC(simple_mod)
    aic_full   <- stats::AIC(full_mod)

    aic_table <- data.frame(
      model = c("simple", "full"),
      AIC = c(aic_simple, aic_full)
    )

    if (aic_simple <= aic_full) {

      selected_model <- "simple"
      best_mod <- simple_mod

      message(
        "Selected model: simpler model.\n",
        "Reason: lower (or equal) AIC compared to the full model.\n",
        sprintf("AIC(simple) = %.2f, AIC(full) = %.2f.", aic_simple, aic_full)
      )

    } else {

      selected_model <- "full"
      best_mod <- full_mod

      message(
        "Selected model: full model.\n",
        sprintf("AIC(full) = %.2f is lower than AIC(simple) = %.2f.", aic_full, aic_simple)
      )
    }
  }

  list(
    combined_data = comb,
    simple_mod = simple_mod,
    full_mod = full_mod,
    aic_table = aic_table,
    selected_model = selected_model,
    gam_mod = best_mod
  )
}
###############################################################
#' Sample Size estimation function using uniroot
#'
#' @param target_power Numeric value specifying the desired statistical power.
#' @param logmean Numeric value representing the log of the mean abundance.
#' @param abs_lfc Numeric value representing the absolute log fold change.
#' @param model A fitted GAM/SCAM model used to predict statistical power.
#' @param xmin Numeric value giving the minimum sample size considered in the search.
#' @param xmax Numeric value giving the maximum sample size considered in the search.
#' @param maxiter maximum number of iterations
#' @param max_report maximum group sample size to be predicted.
#'                   Any predicted sample size that exceed max_report will be
#'                   reported as "> max_report"
#' @return A numeric value corresponding to the estimated sample size required
#'        to achieve the target power.
#' @export

uniroot_ss <- function(target_power, logmean, abs_lfc, model,
                        xmin = log2(10), xmax = log2(5000),
                        maxiter = 10000,
                        max_report = 2000) {

  f <- function(ss) {

    data <- data.frame(
      logsample_size = ss,
      logmean = logmean,
      abs_lfc = abs_lfc
    )

    pred <- predict(
      model,
      type = "response",
      newdata = data[rep(1, 2), ]
    )

    pred[[1]] - target_power
  }

  f_min <- f(xmin)
  f_max <- f(xmax)

  pred_min <- f_min + target_power
  pred_max <- f_max + target_power

  if (is.na(f_min) || is.na(f_max)) {
    return(list(
      sample_size_per_group = NA,
      conclusion = "Sample size could not be estimated because the model returned NA predictions.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }

  if (f_min >= 0) {
    ss_est <- 2^xmin

    return(list(
      sample_size_per_group = ss_est,
      conclusion = "Target power is already achieved at the minimum sample size.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }

  if (f_max < 0) {
    return(list(
      sample_size_per_group = NA,
      conclusion = paste0(
        "No sign change was found. The predicted power does not reach the target power ",
        "within the specified sample-size interval."
      ),
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }

  root <- tryCatch(
    stats::uniroot(
      f,
      interval = c(xmin, xmax),
      maxiter = maxiter
    )$root,
    error = function(e) NA_real_
  )

  if (is.na(root)) {
    return(list(
      sample_size_per_group = NA,
      conclusion = "Root finding failed even though the endpoint checks suggested a solution.",
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }

  ss_est <- 2^root

  if (ss_est > max_report) {
    return(list(
      sample_size_per_group = paste0("> ", max_report),
      conclusion = paste0(
        "The estimated sample size exceeds ",
        max_report,
        "."
      ),
      predicted_power_min_n = pred_min,
      predicted_power_max_n = pred_max
    ))
  }

  list(
    sample_size_per_group = ss_est,
    conclusion = "Success.",
    predicted_power_min_n = pred_min,
    predicted_power_max_n = pred_max
  )
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


