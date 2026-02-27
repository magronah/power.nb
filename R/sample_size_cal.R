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
uniroot_ss =  function(target_power,logmean, abs_lfc,model,xmin,xmax){

  root <- uniroot(function(ss) {

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
  #2^root
}


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
#
#
#
#
#' Power function
#'
#' @param deseq_est_list a list containing fold change, pvalues and other estimates from `DESeq2`
#' @param true_logfoldchange  list containing simulated log fold change used  for simulating the count data
#' @param true_logmean list containing simulated log mean counte used  for simulating the count data
#' @param nsample_vec vector of sample sizes
#' @param alpha_level sign containing simulated log fold change used  for simulating the count data
#' @param notu number of OTUs
#'
#'
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

  #' fit GAM with covariates as tensor product (ie,interaction between
  #' log mean abundance and absolute log fold changes
  #' and then a spline for the sample sizes
  #' log mean abundance and log fold changes are related directly,hence the
  #' interaction but sample size is not quite related to log mean abundance
  #' and log fold changes directly

  comb$logsample_size = log2(comb$sample_size)

  fit_3d <- scam::scam(pval_reject ~ s(logmean, abs_lfc,bs="tedmi") +
                                     s(logsample_size,bs="mpi"),
                                     data = comb, family = binomial)

  list(combined_data=comb, gam_mod = fit_3d)
}


