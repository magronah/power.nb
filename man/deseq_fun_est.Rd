% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{deseq_fun_est}
\alias{deseq_fun_est}
\title{Fold change and p-value estimations for a many simulations}
\usage{
deseq_fun_est(
  metadata_list,
  countdata_list,
  num_cores = 2,
  ref_name = "control"
)
}
\arguments{
\item{metadata_list}{: list of metadata}

\item{countdata_list}{: list of otu count data}

\item{num_cores}{: number of cores}

\item{ref_name}{reference for fold change calculation}
}
\value{
A list
logfoldchange log fold change estimates

logmean  is the log mean count for taxa
(arithmetic mean for taxa across all subjects)

dispersion: dispersion estimates for each taxa

deseq_estimate is a  dataframe containing results from deseq
baseMean,log2FoldChange, lfcSE, pvalue, padj

normalised_count is the normalised count data
}
\description{
Fold change and p-value estimations for a many simulations
}
