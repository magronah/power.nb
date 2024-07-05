#' Computes the optimal number of gaussian components for log mean count
#'
#' The number of gaussian components is determined using the
#' using parametric bootstrap
#'
#' @param logmean vector of log mean abundances of taxa
#' @param sig: significance level to compare against p-value
#' @param max.comp: maximum number of Gaussian components to compare sequentially
#' @param max.boot: maximum number of bootstraps simulations
#'
#' @return The best number of components for fitting the distribution of log mean abundance
#' @export
#'
#' @examples
#' logmean  = rnorm(100)
#' optimal.comp(logmean,sig=0.05,max.comp=4,max.boot=100)

optimal.comp <- function(logmean,sig=0.05,max.comp=4,max.boot=100){
  a <- mixtools::boot.comp(y = logmean, max.comp = max.comp, B = max.boot,
                           mix.type = "normalmix",epsilon = 1e-3)

  pvals=a$p.values; l=length(pvals)
  ncomp <- if (pvals[length(pvals)]<sig) length(pvals)+1 else length(pvals)
  return(ncomp)
}

###############################################################################
#' Filter to remove low abundant taxa
#'
#'
#' Filter to retain only taxa with at least `abund_thresh` counts in at least `sample_thresh` samples
#'
#' @param countdata  otu table
#' @param metadata   dataframe with 2 rows sample names and group names
#' @param abund_thresh  minimum number of taxa abundance threshold
#' @param sample_thresh    minimum number of sample threshold
#' @return filtered otu count data
#'
#' @export
#'

filter_low_count <- function(countdata, metadata,abund_thresh=5, sample_thresh=3){

  ## sanity check
  if(all((metadata$Samples)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countdata,metadata, ~Groups)
  keep <- rowSums(counts(dds) >= abund_thresh) >= sample_thresh

  dds    =   dds[keep,]
  data.frame(counts(dds))
}

########################################################
#' Estimate log fold changes using `DESeq2`.
#'
#' It is stongly recommended to keep defaults for `minReplicatesForReplace=Inf`,
#' `cooksCutoff=TRUE`, `independentFiltering=TRUE`, especially when
#' estimating fold change in order to fit the mixture of Gaussian distributions.
#'
#'
#' @param countdata otu count data
#' @param metadata dataframe with 2 rows sample names and group names
#' @param alpha_level significance level
#' @param ref_name: reference for fold change calculation
#' @param minReplicatesForReplace: DESeq2's parameter to control minimum number of
#'     replicates needed for the replacement of outliers during dispersion estimation.
#' @param cooksCutoff: DESeq2's outlier removal or shrinkage.
#' @param independentFiltering: DESeq2's independent filtering.
#' @param shrinkage_method:  DESeq2's shrinkage method
#'
#' @return  A list
#'
#' logfoldchange log fold change estimates
#'
#' logmean  is the log mean count for taxa
#' (arithmetic mean for taxa across all subjects)
#'
#' dispersion: dispersion estimates for each taxa
#'
#' deseq_estimate is a  dataframe containing results from deseq
#'         baseMean,log2FoldChange, lfcSE, pvalue, padj
#'
#' normalised_count is the normalised count data
#' @export
#'
#' @examples
#' countdata = matrix(100, ncol = ..)
#' deseqfun(countdata,metadata)
deseqfun <- function(countdata,metadata,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf,
                     cooksCutoff = TRUE,
                     independentFiltering = TRUE,
                     shrinkage_method="normal"){

  #check otu table is in otu by samples format
  if(all((metadata$Samples)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }

  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  metadata= metadata[keep, ]

  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,metadata, ~Groups)
  dds$Groups <- relevel(dds$Groups, ref = ref_name)

  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace)

  res <- results(dds, cooksCutoff=cooksCutoff,
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)

  reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)

  deseq_est = data.frame(reslt)
  disp = dispersions(dds)

  logfoldchange = deseq_est$log2FoldChange
  names(logfoldchange) = rownames(deseq_est)

  normalised_count     =  counts(dds, normalized=TRUE)
  logmean = log2(rowMeans(normalised_count))

  list(logfoldchange = logfoldchange,
       logmean = logmean,
       dispersion = disp,
       deseq_estimate=deseq_est,
       normalised_count = normalised_count)
}



#' General-purpose log-likelihood function, vectorized sum(pars*x^i)
#'
#'
#' @param pars parameters
#' @param x log mean count
#'
#' @return values representing output for the polynomial fuction (f(x))
#'
#'

polyfun <- function(pars, x) {
  ord <- length(pars)
  xmat <- sapply(0:(ord-1), function(i) x^i)
  xmat |> sweep(MARGIN = 2, FUN = "*", pars) |> rowSums()
}


#' generate normal mixture parameters (prob vector, mean vector, sd vector
#' for a specified set of 'x' values (logmean)
#'
#' @param x independent variable
#' @param pars parameter vector: first logit-probs (np-1),
#'             then mean parameters (2 per component: intercepts, slopes),
#'             then var parameters (varord + 1 per component: intercepts,
#'             slopes, quad coeffs, etc.)
#' @param np number of components in mixture
#' @param sd_ord order of logsd model (2 = quadratic)
#'
#'@return   A list
#'
#'      probs: mixture proportions ()
#'
#'      muvals: mean values
#'
#'      sdvals: standard deviation values


genmixpars <- function(x, pars, np = 2, sd_ord = 2){
  ## complain if pars is wrong length
  expected_pars <- (np-1)+2*np+(sd_ord+1)*np
  if (length(pars) != expected_pars) {
    stop(sprintf("number of pars (%d) != expected (%d) (np = %d, sd_ord = %d)",
                 length(pars), expected_pars, np, sd_ord))
  }
  ## softmax function
  cp <- 1  ## parameter index
  probs <- exp(c(0,pars[1:(np-1)])) # first is zero
  probs <- probs/sum(probs)
  cp <- cp + (np-1)
  ## mu: always linear (2 params per component)
  ## intercepts first, then slopes
  mupars <- pars[cp:(cp + 2*np - 1)]
  mupars <- matrix(mupars, ncol = 2)
  muvals <- apply(mupars, 1, polyfun, x = x)
  cp <- cp + 2*np
  ## logsd: similar
  logsdpars <- pars[cp:(cp + (sd_ord+1)*np - 1)]
  logsdpars <- matrix(logsdpars, ncol = (sd_ord+1))
  logsdvals <- apply(logsdpars, 1, polyfun, x = x)
  list(probs = probs, muvals = muvals, sdvals = exp(logsdvals))
}

#' general-purpose normal-mixture deviate generator: takes _matrices_
#' of probabilities, means, sds
#'
#' @param n number of observations
#' @param probs mixture proportions
#' @param muvals values for the mean
#' @param sdvals values for the standard deviation
#'
#' @return simulations
#' @export
#'
rnormmix0 <- function(n, probs, muvals, sdvals) {
  np <- length(probs)
  component <- sample(np, size = n, prob = probs, replace = TRUE)
  inds <- cbind(seq(n), component)
  rnorm(n, mean = muvals[inds], sd = sdvals[inds])
}


#' Simulating from a mixture of Gaussian
#'
#' @param par parameters (mean, standard deviation and mixture proportion of the mixture of Gaussian)
#' @param logmean log mean count of taxa
#' @param ... other parameters taken by `genmixpars`
#'
#' @return random values from a  a mixture of Gaussian
#' @export
#'
myrnormmix <- function(par, logmean, ...) {
  g0 <- genmixpars(logmean, par, ...)
  do.call(rnormmix0, c(list(n = length(logmean)), g0))
}

#' Density function for the mixture of Gaussian distributions
#'
#' takes pars as three vectors (prob, mean, sd), on constrained scale
## i.e. prob (0,1); mean; sdvals (0, Inf)
#'
#'
#' @param x vector
#' @param probs mixture proportions
#' @param muvals values of the mean
#' @param sdvals values for standard deviataions
#' @param log log scale
#'
#' @return likelihood
#' @export
#'

dnormmix0 <- function(x, probs, muvals, sdvals, log = FALSE) {
  np <- length(probs)
  pmat <- matrix(NA, nrow = length(x), ncol = np)
  for (i in 1:np) {
    pmat[,i] <- dnorm(x, mean = muvals[,i], sd = sdvals[,i])
  }
  lik <- pmat |> sweep(MARGIN = 2, FUN = "*", probs) |> rowSums()
  if (log) log(lik) else lik
}

#' Takes a single vector of parameters on the unconstrained scale
#  (i.e. softmax(prob), mean, log(sd))
#'
#'
#' @param x
#' @param par   softmax(prob), mean, log(sd) on unconstrained scale
#' @param logmean
#' @param ... other input calues taken by `genmixpars`
#' @param log log scale
#'
#' @return
#'
dnormmix <- function(x, par, logmean, ..., log = FALSE) {
  g0 <- genmixpars(logmean, par, ...)
  do.call(dnormmix0, c(list(x), g0, list(log = log)))
}

#' Objective function
#'
#' @param par parameters of the mixture of Gaussian
#' @param vals values of log fold change
#' @param logmean log mean count
#' @param np  number of Gaussian components
#' @param sd_ord order of polynomial function to model standard deviation parameters
#'               (1 - linear function and 2- quad)
#'
#' @return value for the objective function
#'
nllfun <- function(par, vals, logmean, np, sd_ord) {
  -sum(dnormmix(x = vals, par, logmean, np = np, sd_ord = sd_ord, log = TRUE))
}


#' mm
#'
#' @param np  number of Gaussian components
#' @param sd_ord order of polynomial function to model standard deviation parameters
#'               (1 - linear function and 2- quadratic function)
#'
#' @return
#'
#'
#'
gen_parnames <- function(np, sd_ord) {
  ## taken from contr.poly:
  sdlabs <- c(".1", ".L", ".Q", ".C", paste0(".", 4:10))[1:(sd_ord+1)]
  c(paste0("logitprob_", 1:(np-1)),
    c(outer(1:np, c("int", "slope"), function(x,y) sprintf("mu_%s_%d", y, x))),
    c(outer(1:np, sdlabs, function(x,y) sprintf("logsd_%s_%d", y, x))))
}



#############################################################################
#' Fold change and p-value estimations for a many simulations
#'
#' @param metadata_list : list of metadata
#' @param countdata_list : list of otu count data
#' @param num_cores : number of cores
#' @param ref_name reference for fold change calculation
#'
#' @return  A list
#'  logfoldchange log fold change estimates
#'
#' logmean  is the log mean count for taxa
#' (arithmetic mean for taxa across all subjects)
#'
#' dispersion: dispersion estimates for each taxa
#'
#' deseq_estimate is a  dataframe containing results from deseq
#'         baseMean,log2FoldChange, lfcSE, pvalue, padj
#'
#' normalised_count is the normalised count data
#' @export
#'
#'
deseq_fun_est <-function(metadata_list,  countdata_list,
                         num_cores=2, ref_name= "control"){

  registerDoParallel(cores = num_cores)
  l = length(countdata_list)

  dds  =  foreach(i= 1:l, .packages = "DESeq2", .export = "deseqfun") %dopar%{
    countdata =  countdata_list[[i]]
    metadata  =  metadata_list[[i]]
    stopifnot(!is.na(sum(countdata)))
    stopifnot(colnames(countdata) == metadata$Samples)
    deseqfun(countdata,metadata, ref_name=ref_name,
             minReplicatesForReplace = Inf,
             cooksCutoff = FALSE,
             independentFiltering = FALSE)
  }

  stopImplicitCluster()
  unregister_dopar()
  names(dds)  =  names(countdata_list)
  dds
}

#############################################################
#' Contour plot for showing predicted power
#'
#' @param combined_data data used for fitting gam
#' @param power_estimate predicted power
#' @param cont_breaks breaks for contour plot
#'
#' @return ggplot object
#' @export
#'
#'
contour_plot_fun <- function(combined_data,
                             power_estimate,
                             cont_breaks){

  combined_data$pvalue_reject <- factor(combined_data$pval_reject)

  gg_2dimc <- (ggplot(combined_data)
               + aes(lmean_abund, abs_lfc)
               + rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)"))
               + ylab(TeX("|$\\log_2$(fold change)|"))
               + scale_colour_manual(values = c("black", "red"))
               + geom_contour(data = power_estimate,
                              aes(z=power),lwd=1,
                              breaks = cont_breaks)
               + geom_label_contour(data = power_estimate,
                                    aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               )

  )
  gg_2dimc

}




