#' Computes the optimal number of gaussian components for log mean count
#'
#' The number of gaussian components is determined using the
#' using parametric bootstrap
#'
#' @param logmean vector of log mean abundances of taxa
#' @param sig significance level to compare against p-value
#' @param max.comp maximum number of Gaussian components to compare sequentially
#' @param max.boot maximum number of bootstraps simulations
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
#' @param sample_colname  column names of the samples
#' @param group_colname   column names of the groups or conditions
#' @return filtered otu count data
#'
#' @export
#'

filter_low_count <- function (countdata, metadata, abund_thresh = 5, sample_thresh = 3,
                              sample_colname, group_colname)
{
  stopifnot(is.data.frame(metadata) || is.matrix(metadata))
  stopifnot(is.data.frame(countdata) || is.matrix(countdata))
  if (is.null(sample_colname) || !(sample_colname %in% names(metadata))) {
    stop("Could not find a sample ID column in metadata.\n
         Pass `sample_colname=` explicitly.")
  }
  if (is.null(group_colname) || !(group_colname %in% names(metadata))) {
    stop("Could not find a group/condition column in metadata.\n
         Pass `group_colname=` explicitly.")
  }
  if (all((metadata[[sample_colname]]) == colnames(countdata))) {
    countdata = t(countdata)
    }

  metadata[[group_colname]]  = as.factor(metadata[[group_colname]])
  design_formula <- stats::as.formula(paste("~", group_colname))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata, design = design_formula)
  keep <- rowSums(DESeq2::counts(dds) >= abund_thresh) >= sample_thresh
  dds = dds[keep, ]
  as.data.frame(DESeq2::counts(dds))
}

########################################################
#' Estimate log fold changes using `DESeq2`.
#'
#' This function estimates log fold changes (LFC) for microbiome count data using the DESeq2 package.
#' It is strongly recommended to keep the following defaults:
#' - `minReplicatesForReplace=Inf`
#' - `cooksCutoff=TRUE`
#' - `independentFiltering=TRUE`
#' These options are particularly useful when estimating fold changes to fit the mixture of Gaussian distributions.
#' @param countdata A matrix of OTU count data where rows represent taxa and columns represent samples.
#' @param metadata A dataframe containing sample information with two rows: one for sample names and one for group names.
#' @param alpha_level The significance level for determining differential expression. Default is 0.1.
#' @param sample_colname  column names of the samples
#' @param group_colname   column names of the groups or conditions
#' @param ref_name The reference group for calculating fold changes.
#' @param minReplicatesForReplace DESeq2's parameter to control the minimum number of replicates required for replacing outliers during dispersion estimation. Default is `Inf` (no replacement).
#' @param cooksCutoff DESeq2's parameter for removing outliers based on Cook's distance. Default is `TRUE` (outlier removal enabled).
#' @param independentFiltering DESeq2's parameter for independent filtering. Default is `TRUE`.
#' @param shrinkage_method DESeq2's shrinkage method for fold changes. Default is `"normal"`. Other options include `"apeglm"` or `"ashr"`.
#' @importFrom DESeq2 counts
#' @return A list containing the following elements:
#' - `logfoldchange`: A vector of log fold change estimates for each taxa.
#' - `dispersion`: A vector of dispersion estimates for each taxa.
#' - `deseq_estimate`: A dataframe containing DESeq2 results, including `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, and `padj`.
#' - `normalised_count`: A matrix of normalized count data.
#' - `dds`: DESeq object
#'
#' @export
#'
#' @examples
#' # Example usage
#' set.seed(101)
#' nr = 10; nc =50
#' countdata <- matrix(rpois(500, 3), ncol = nc, nrow = nr)
#' # Simulated OTU count data with 50 taxa and 10 samples
#' countdata  <- as.data.frame(countdata)
#' rownames(countdata)  <- paste("Sample", 1:nr, sep = "_")
#' colnames(countdata) <- paste("otu", 1:nc, sep = "_")
#' metadata <- data.frame(Samples = paste("Sample", 1:nr, sep = "_"),
#'              Groups = rep(c("Control", "Treatment"), each = 5))
#' sample_colname = "Samples"
#' group_colname  = "Groups"
#'
#' result <- deseqfun(countdata, metadata, ref_name = "Control",
#'                     minReplicatesForReplace = Inf,
#'                     cooksCutoff = TRUE,
#'                     sample_colname = "Samples",
#'                     group_colname  = "Groups",
#'                     independentFiltering = TRUE,
#'                     shrinkage_method="normal")
#'
#' # Examine the results
#' result$logfoldchange  # Log fold changes
#' result$logmean  # Log mean counts
#' result$dispersion  # Dispersion estimates
#' result$deseq_estimate  # DESeq2 results
#' result$normalised_count  # Normalized count data
#'

deseqfun <- function (countdata, metadata, alpha_level = 0.1,
                      ref_name,
                      group_colname,
                      sample_colname,
                      minReplicatesForReplace = Inf,
                      cooksCutoff = TRUE,
                      independentFiltering = TRUE,
                      shrinkage_method = "normal"){

  stopifnot(is.data.frame(metadata) || is.matrix(metadata))
  stopifnot(is.data.frame(countdata) || is.matrix(countdata))

  if (is.null(sample_colname) || !(sample_colname %in% names(metadata))) {
    stop("Could not find a sample ID column in metadata.\n
         Pass `sample_colname=` explicitly.")
  }

  if (is.null(group_colname) || !(group_colname %in% names(metadata))) {
    stop("Could not find a group/condition column in metadata.\n
         Pass `group_colname=` explicitly.")
  }

  if (!all(metadata[[sample_colname]] %in% colnames(countdata))) {
    countdata <- t(countdata)
  }

  keep <- (colSums(countdata) > 0)

  if(sum(keep) < length(keep)){
    countdata = countdata[, keep]
    #metadata  = metadata[metadata[[sample_colname]] %in% rownames(countdata), ]
  }

  metadata[[group_colname]]  = as.factor(metadata[[group_colname]])
  design_formula <- stats::as.formula(paste("~", group_colname))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = design_formula)

  dds[[group_colname]] <- relevel(dds[[group_colname]], ref = ref_name)

  dds <- DESeq2::DESeq(dds, sfType = "poscounts",
                       minReplicatesForReplace = minReplicatesForReplace)

  res <- DESeq2::results(dds, cooksCutoff = cooksCutoff,
                         independentFiltering = independentFiltering,
                         alpha = alpha_level)

  reslt <- DESeq2::lfcShrink(dds, res = res, coef = 2, type = shrinkage_method)
  deseq_est = data.frame(reslt)
  disp = DESeq2::dispersions(dds)
  logfoldchange = deseq_est$log2FoldChange
  names(logfoldchange) = rownames(deseq_est)
  normalised_count = DESeq2::counts(dds, normalized = TRUE)
  #logmean = log2(rowMeans(normalised_count))

  list(logfoldchange  =  logfoldchange,
       #logmean      =  logmean,
       dispersion     =  disp,
       deseq_estimate =  deseq_est,
       normalised_count =  normalised_count,
       deseq_object  =  dds)
}


#' General-purpose log-likelihood function, vectorized sum(pars*x^i)
#'
#'
#' @param pars parameters
#' @param x log mean count
#'
#' @return values representing output for the polynomial fuction (f(x))
#' @export
#' @examples
#' polyfun(pars = c(1, 2, 3), x = 1:5)
#' polyfun(pars = c(1, 0, 3), x = 1)

polyfun <- function(pars, x) {
  ord <- length(pars)
  xmat <- sapply(0:(ord-1), function(i) x^i)
  if (!is.matrix(xmat)) xmat <- matrix(xmat, nrow = 1)
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
#'@export


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
  stats::rnorm(n, mean = muvals[inds], sd = sdvals[inds])
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


#' Density of a Normal Mixture Model
#'
#' This function calculates the density of a normal mixture model for a given vector of parameters
#' on the unconstrained scale (softmax(prob), mean, log(sd)).
#'
#' @param x Numeric vector of values at which to evaluate the density.
#' @param par A vector of parameters on the unconstrained scale, including:
#'   \itemize{
#'     \item \code{softmax(prob)}: Mixture probabilities (on the unconstrained scale, transformed via softmax).
#'     \item \code{mean}: Means of the normal components.
#'     \item \code{log(sd)}: Logarithms of standard deviations of the normal components.
#'   }
#' @param logmean Numeric value representing the log of the mean parameter.
#' @param ... Additional arguments passed to the \code{genmixpars} function. Defaults: two components (\code{np = 2}) and a quadratic model for standard deviation parameters (\code{sd_ord = 2}).
#' @param log Logical. If \code{TRUE}, the logarithm of the density is returned. Default is \code{FALSE}.
#'
#' @return A numeric vector of density values (or log-density values if \code{log = TRUE}) for the mixture model.
#'
#' @export
#'
#' @examples
#' # Example parameters
#' x <- seq(-3, 3, length.out = 100)
#' ## par <- c(-0.5, 0.5, log(0.8), log(1.2))  # Example: softmax probabilities, mean, log(sd)
#' set.seed(101); par <- rnorm(11)
#' logmean <- rep(0.1, length(x))  ## constant log mean
#'
#' # Calculate density
#' density <- dnormmix(x, par, logmean)
#'
#' # Calculate log-density
#' log_density <- dnormmix(x, par, logmean, log = TRUE)

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


#' Generate Parameter Names for Mixture Model
#'
#' This function generates parameter names for a Gaussian mixture model based on the number
#' of components (\code{np}) and the order of the polynomial function (\code{sd_ord})
#' used to model the standard deviation parameters.
#'
#' @param np Integer. The number of Gaussian components in the mixture model.
#' @param sd_ord Integer. The order of the polynomial function used to model the
#'   standard deviation parameters. Possible values are:
#'   \itemize{
#'     \item \code{1}: Linear function.
#'     \item \code{2}: Quadratic function.
#'   }
#'
#' @return A character vector of parameter names, including:
#'   \itemize{
#'     \item Logit-transformed probabilities (\code{logitprob_1}, ..., \code{logitprob_(np-1)}).
#'     \item Mean parameters (\code{mu_int_1}, \code{mu_slope_1}, ..., for each component).
#'     \item Log-transformed standard deviations (\code{logsd_.1_1}, \code{logsd_.L_1}, ...,
#'           depending on \code{sd_ord} and \code{np}).
#'   }
#'
#' @export
#'
#' @examples
#' # Generate parameter names for a 3-component mixture with linear standard deviation function
#' gen_parnames(np = 3, sd_ord = 1)
#'
#' # Generate parameter names for a 4-component mixture with quadratic standard deviation function
#' gen_parnames(np = 4, sd_ord = 2)

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
#' @param sample_colname  column names of the samples
#' @param group_colname   column names of the groups or conditions
#' @param alpha_level The significance level for determining differential expression. Default is 0.1.
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
                         alpha_level = 0.1,
                         group_colname,
                         sample_colname,
                         num_cores=2, ref_name= "control"){

  registerDoParallel(cores = num_cores)
  l = length(countdata_list)
  #i  <- NULL
  dds  =  foreach(i= 1:l, .packages = "DESeq2", .export = "deseqfun") %dopar%{
    countdata =  countdata_list[[i]]
    metadata  =  metadata_list[[i]]
    stopifnot(!is.na(sum(countdata)))
    stopifnot(all(colnames(countdata) %in% metadata$Samples))

    deseqfun(countdata,metadata,alpha_level = alpha_level,
             group_colname = group_colname,
             sample_colname = sample_colname,
             ref_name=ref_name,
             minReplicatesForReplace = Inf,
             cooksCutoff = FALSE,
             independentFiltering = FALSE)
  }

  stopImplicitCluster()
  # unregister_dopar()
  names(dds)  =  names(countdata_list)
  dds
}

#############################################################
#' Contour plot for showing predicted power
#'
#' @param combined_data data used for fitting gam
#' @param power_estimate predicted power
#' @param cont_breaks breaks for contour plot
#' @importFrom ggplot2 ggplot aes geom_contour geom_point xlab ylab scale_colour_manual
#' @importFrom ggplot2 after_stat
#' @return ggplot2 object
#' @importFrom latex2exp TeX
#' @export
#'
#'
contour_plot_fun <- function(combined_data,
                             power_estimate,
                             cont_breaks){

    ## utils::globalVariables(c("lmean_abund", "abs_lfc"))
    ## deal with code checking: not sure why 'globalVariables' not working
    ##lmean_abund <- abs_lfc <- power <- power_estimate <-
    ##    pvalue_reject <- level <- NULL

  combined_data$pvalue_reject <- factor(combined_data$pval_reject)

  gg_2dimc <- (ggplot(combined_data)
               + aes(lmean_abund, abs_lfc)
               + ggrastr::rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)"))
               + ylab(TeX("|$\\log_2$(fold change)|"))
               + scale_colour_manual(values = c("black", "red"))
               + geom_contour(data = power_estimate,
                              aes(z=power),lwd=1,
                              breaks = cont_breaks)
               + metR::geom_label_contour(data = power_estimate,
                                    aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               )

  )
  gg_2dimc

}



#' Extract specified data from a list of datasets
#'
#' This function extracts a specific component (data) from a list of datasets.
#' The component to extract is specified by the `extract_name` parameter, and
#' the function returns a list containing the extracted data from each dataset.
#'
#' @param dataset_list A list of datasets from which data will be extracted.
#' Each element of the list is assumed to be a dataset (typically a list or dataframe).
#' @param extract_name A string representing the name of the component or column
#' to be extracted from each dataset in the `dataset_list`. The function looks
#' for this name within each dataset.
#'
#' @return A list containing the extracted data. Each element corresponds to
#' the extracted component from the datasets in `dataset_list`. The names of
#' the list elements are taken from the names of `dataset_list`.
#'
#' @examples
#' # Example dataset list
#' dataset1 <- list(countdata = matrix(1:9, nrow = 3), metadata = data.frame(id = 1:3))
#' dataset2 <- list(countdata = matrix(10:18, nrow = 3), metadata = data.frame(id = 4:6))
#' dataset_list <- list(dataset1 = dataset1, dataset2 = dataset2)
#'
#' # Extract 'countdata' from each dataset in the list
#' result <- read_data(dataset_list, "countdata")
#' print(result)
#'
#' @export

read_data <- function(dataset_list, extract_name){
  extract_data_list= list()
  for(n in 1:length(dataset_list)){
    dataset <- dataset_list[[n]]
    extract_data_list[[n]] <- dataset[[extract_name]]
  }
  names(extract_data_list) <- names(dataset_list)
  extract_data_list
}
