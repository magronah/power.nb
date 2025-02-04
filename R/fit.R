#' Fit a mixture of Gaussian Distributions to log mean count of taxa.
#'
#' The optimal number of components to fit is chosen using parametric bootstrap method
#'
#' @param logmean vector of log mean count of taxa
#' @param sig significance level to compare against p-value to be used for
#'             parametric bootstrap calculation
#' @param max.comp maximum number of Gaussian components to compare sequentially
#' @param max.boot maximum number of bootstraps simulations
#' @return A list containing the optimal number of Gaussian components fitted; and mean and
#'          variance parameter estimates from the fit


#' @export
#'
#' @examples
#' logmean  = rnorm(100)
#' logmean_fit(logmean,sig=0.05,max.comp=4,max.boot=100)
#'
logmean_fit <- function(logmean,sig=0.05,max.comp=4,max.boot=100){

  ncomp_opt = optimal.comp(logmean,sig,max.comp,max.boot)
  if(length(ncomp_opt) == 0){stop("zero number of component")}

  if(ncomp_opt == 1){
    mixmdl = fitdistrplus::fitdist(logmean,"norm")

    param  = data.frame(sigma =  coef(mixmdl)[["sd"]],
                        mu  = coef(mixmdl)[["mean"]])

    sim_mean = rnorm(n = length(logmean),
                     mean = param$mu,
                     sd   = param$sigma)
  }else{
    mixmdl = normalmixEM(logmean, k=ncomp_opt, maxrestarts=1e3)

    param  = data.frame(lambda =  mixmdl$lambda,
                        sigma =  mixmdl$sigma,
                        mu   =  mixmdl$mu)

    sim_mean=mixtools::rnormmix(n=length(logmean),lambda=mixmdl$lambda,
                                mu =mixmdl$mu,
                                sigma=mixmdl$sigma)
  }

  list(logmean_param=param,components=ncomp_opt)
}


##################################################################
#' Fit the non-linear function to dispersion estimates
#'
#' Dispersion are estimated from the `DESeq2` package. The function fitted is of the
#' form  `a + b/(mean count)` where `a` represents the asymptotic dispersion level
#' for high abundance taxa, and `b` captures additional dispersion variability.
#'
#' @param dispersion dispersion estimates from deseq
#' @param logmean  vector of log mean abundance
#'
#' @return A list containing estimates for `a` and `b` and confidence intervals
#' @export
#'
#' @examples
#' logmean    =  rnorm(100)
#' dispersion =  abs(rnorm(100))
#' dispersion_fit(dispersion,logmean)

dispersion_fit <- function(dispersion,logmean){

  dat = data.frame(dispersion=dispersion,mean_abund = 2^logmean)

  param= nls(dispersion~ asymptDisp + extraPois/mean_abund, data = dat,
             start = list(asymptDisp = 0.1, extraPois = 0.1))

  dd=list(param=data.frame(asymptDisp= coef(param)[[1]],
                           extraPois = coef(param)[[2]]),
          confint=confint(param))
  dd
}

##################################################################
#' Fit a mixture of Gaussian distributions to log fold change
#'
#' The standard deviation parameters are modeled either by linear or quadratic functions of log mean count and the mean parameter is modeled
#' by linear functions of log mean count
#'
#' @param logmean vector of log mean abundance
#' @param logfoldchange vector of log fold change
#' @param ncore number of cores to use
#' @param max_sd_ord the maximum order of polynomial function to fit to
#'                    standard deviation parameter.
#'                    This must be either 1 (linear) or 2(quadratic)
#' @param max_np maximum number of Gaussian components to check for
#' @param minval minimum value for DEoptim search
#' @param maxval maximum value for DEoptim search
#' @param NP the number of population members for DEoptim
#' @param itermax maximum number of iterations
#' @param seed  seed value
#'
#' @return A list.
#'
#'        par is a vector of the estimates of the mixture proportion,
#'        the mean and standard deviation parameters,
#'
#'        np is the number of gaussian components fitted
#'
#'        sd_ord is the order for the function for the standard deviation
#'
#'        aic is the aic of the best fit
#'
#' @export
#'
#' @examples
#' logmean        =  rnorm(100)
#' logfoldchange  =  rnorm(100)
#' logfoldchange_fit(logmean,logfoldchange)
#'
logfoldchange_fit = function(logmean,logfoldchange,ncore = 2,
                             max_sd_ord = 2, max_np=5,
                             minval = -5, maxval = 5,
                             itermax = 100,NP=800, seed  = 100){

  res  =  list(); cnt = 0

  for(sd_ord in 1:max_sd_ord){

    for(np in 2:max_np){

      l   =  (np-1)+2*np+(sd_ord+1)*np
      set.seed(seed)

      cl <- makeCluster(ncore)
      fun_names <- c("polyfun", "dnormmix", "dnormmix0", "genmixpars")
      clusterExport(cl, c("logmean", "logfoldchange", fun_names))

      opt  <- DEoptim(fn = nllfun, lower = rep(minval, l), upper = rep(maxval, l),
                      vals = logfoldchange, logmean = logmean, np = np, sd_ord = sd_ord,
                      control = DEoptim.control(NP = NP, itermax = itermax,
                                                cluster = cl))

      pn <- gen_parnames(np = np, sd_ord = sd_ord)
      par <- stats::setNames(opt$optim$bestmem, pn)
      aic = 2*length(par) + 2*opt$optim$bestval


      cnt = cnt + 1
      res[[cnt]]  =  list(par = par, np = np, sd_ord= sd_ord,aic = aic)
    }
  }
  aic_values <- sapply(res, function(x) x$aic)
  index_of_min_aic <- which.min(aic_values)

  res[[index_of_min_aic]]
}

##################################################################
#' Title
#'
#' @param deseq_est_list a list containing fold change, pvalues and other estimates from `DESeq2`
#' @param true_lfoldchange_list list containing simulated log fold change used  for simulating the count data
#' @param true_lmean_list   list containing simulated log mean count used  for simulating the count data
#' @param grid_len     number of grids for
#' @param alpha_level significance level for power calculations
#'
#' @return A list
#'
#' fit_2d is the fitted scam object
#'
#' power_estimate predicted power estimated using fit_2d
#'
#' combined_data tibble containing pvlaues and used for GAM fit

#' @export
#'
gam_fit <- function(deseq_est_list,
                    true_lfoldchange_list,
                    true_lmean_list,
                    grid_len = 50,
                    alpha_level=0.1){

  p_val   = foreach::foreach(k = 1:length(deseq_est_list),.combine = "c") %do%{
    deseq_est_list[[k]]$padj
  }

  pval_reject       =   (!is.na(p_val) & p_val < 0.1)
  true_lfoldchange    =    unlist(true_lfoldchange_list)
  true_lmean_abund    =    unlist(true_lmean_list)

  comb     =   tibble(lmean_abund  =  true_lmean_abund,
                      abs_lfc   =  abs(true_lfoldchange),
                      pval_reject  =  as.numeric(pval_reject))


  fit_2d       =    scam::scam(pval_reject ~ s(lmean_abund, abs_lfc, bs="tedmi"),
                               data = comb, family = binomial)

  pp      =   with(comb,
                   expand.grid(lmean_abund = seq(min(lmean_abund),
                                                 max(lmean_abund),
                                                 length  = grid_len),
                               abs_lfc   =  seq(min(abs_lfc),
                                                max(abs_lfc),
                                                length  =  grid_len)))
  pp$power <- predict(fit_2d, newdata = pp,type = "response")

  p=list(combined_data = comb, power_estimate = pp, fit_2d=fit_2d)
  p
}



