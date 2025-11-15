# Fold change and p-value estimations for a many simulations

Fold change and p-value estimations for a many simulations

## Usage

``` r
deseq_fun_est(
  metadata_list,
  countdata_list,
  alpha_level = 0.1,
  group_colname,
  sample_colname,
  num_cores = 2,
  ref_name = NULL
)
```

## Arguments

- metadata_list:

  : list of metadata

- countdata_list:

  : list of otu count data

- alpha_level:

  The significance level for determining differential expression.
  Default is 0.1.

- group_colname:

  column names of the groups or conditions

- sample_colname:

  column names of the samples

- num_cores:

  : number of cores

- ref_name:

  reference level for fold change calculation. If NULL, the reference
  level is determined automatically — by default, the factor level that
  comes first is used as the reference.

## Value

A list logfoldchange log fold change estimates

logmean is the log mean count for taxa (arithmetic mean for taxa across
all subjects)

dispersion: dispersion estimates for each taxa

deseq_estimate is a dataframe containing results from deseq
baseMean,log2FoldChange, lfcSE, pvalue, padj

normalised_count is the normalised count data
