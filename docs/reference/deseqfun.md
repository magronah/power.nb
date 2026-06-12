# Estimate log fold changes using `DESeq2`.

This function estimates log fold changes (LFC) for microbiome count data
using the DESeq2 package. It is strongly recommended to keep the
following defaults:

- `minReplicatesForReplace=Inf`

- `cooksCutoff=TRUE`

- `independentFiltering=TRUE` These options are particularly useful when
  estimating fold changes to fit the mixture of Gaussian distributions.

## Usage

``` r
deseqfun(
  countdata,
  metadata,
  alpha_level = 0.1,
  ref_name = NULL,
  group_colname,
  sample_colname,
  minReplicatesForReplace = Inf,
  cooksCutoff = TRUE,
  independentFiltering = TRUE,
  shrinkage_method = "normal"
)
```

## Arguments

- countdata:

  A matrix of OTU count data where rows represent taxa and columns
  represent samples.

- metadata:

  A dataframe containing sample information with two rows: one for
  sample names and one for group names.

- alpha_level:

  The significance level for determining differential expression.
  Default is 0.1.

- ref_name:

  reference level for fold change calculation. If NULL, the reference
  level is determined automatically — by default, the factor level that
  comes first is used as the reference.

- group_colname:

  column names of the groups or conditions

- sample_colname:

  column names of the samples

- minReplicatesForReplace:

  DESeq2's parameter to control the minimum number of replicates
  required for replacing outliers during dispersion estimation. Default
  is `Inf` (no replacement).

- cooksCutoff:

  DESeq2's parameter for removing outliers based on Cook's distance.
  Default is `TRUE` (outlier removal enabled).

- independentFiltering:

  DESeq2's parameter for independent filtering. Default is `TRUE`.

- shrinkage_method:

  DESeq2's shrinkage method for fold changes. Default is `"normal"`.
  Other options include `"apeglm"` or `"ashr"`.

## Value

A list containing the following elements:

- `logfoldchange`: A vector of log fold change estimates for each taxa.

- `dispersion`: A vector of dispersion estimates for each taxa.

- `deseq_estimate`: A dataframe containing DESeq2 results, including
  `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, and `padj`.

- `normalised_count`: A matrix of normalized count data.

- `dds`: DESeq object

## Examples

``` r
# \donttest{
# Example usage
set.seed(101)
nr = 10; nc = 35
countdata <- matrix(rpois(350, 3), ncol = nc, nrow = nr)
# Simulated OTU count data with 50 taxa and 10 samples
countdata  <- as.data.frame(countdata)
rownames(countdata)  <- paste("Sample", 1:nr, sep = "_")
colnames(countdata) <- paste("otu", 1:nc, sep = "_")
metadata <- data.frame(Samples = paste("Sample", 1:nr, sep = "_"),
             Groups = rep(c("Control", "Treatment"), each = 5))
sample_colname = "Samples"
group_colname  = "Groups"

result <- deseqfun(countdata, metadata, ref_name = "Control",
                    minReplicatesForReplace = Inf,
                    cooksCutoff = TRUE,
                    sample_colname = "Samples",
                    group_colname  = "Groups",
                    independentFiltering = TRUE,
                    shrinkage_method="normal")
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
#> using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
#> 
#> Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
#> See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
#> Reference: https://doi.org/10.1093/bioinformatics/bty895

# Examine the results
result$logfoldchange  # Log fold changes
#>        otu_1        otu_2        otu_3        otu_4        otu_5        otu_6 
#>  0.006539958 -0.503562487 -0.456073305  0.191041559 -0.077182700  0.087067055 
#>        otu_7        otu_8        otu_9       otu_10       otu_11       otu_12 
#> -0.330786588 -0.362875098 -0.249232240 -0.211065497  0.043266255  0.218668027 
#>       otu_13       otu_14       otu_15       otu_16       otu_17       otu_18 
#>  0.223783495  0.225181598 -0.511200242  0.083506136 -0.069297689  0.004134123 
#>       otu_19       otu_20       otu_21       otu_22       otu_23       otu_24 
#> -0.366580938  0.132955536 -0.140957759  0.003734321  0.042751885  0.167098245 
#>       otu_25       otu_26       otu_27       otu_28       otu_29       otu_30 
#>  0.062900036  0.047768583  0.175193461 -0.117555247  0.530107061 -0.360696897 
#>       otu_31       otu_32       otu_33       otu_34       otu_35 
#>  0.055798248 -0.217883070  0.037330312 -0.250857707  0.077417639 
result$logmean  # Log mean counts
#> NULL
result$dispersion  # Dispersion estimates
#>  [1] 2.559110e-01 1.023829e-02 5.455305e-03 4.622220e-01 2.344478e-01
#>  [6] 1.680165e-02 3.090915e-02 1.082847e-02 3.542332e-01 3.095074e-01
#> [11] 8.680201e-01 1.564384e-02 1.925445e-01 2.516616e-02 1.885939e-02
#> [16] 1.589767e-01 2.708736e-01 1.000000e+01 2.455034e-01 4.030049e-01
#> [21] 2.265075e-02 4.350329e-01 2.168889e-01 1.721007e-01 1.379342e-01
#> [26] 2.120916e-01 2.271774e-02 3.313791e-01 4.038254e-02 1.216709e-02
#> [31] 1.843717e-01 6.521776e-02 3.429549e-01 2.328457e-01 6.359506e-04
result$deseq_estimate  # DESeq2 results
#>        baseMean log2FoldChange      lfcSE        stat     pvalue      padj
#> otu_1  2.508851    0.006539958 0.25499845  0.02567473 0.97951678 0.9879252
#> otu_2  3.467602   -0.503562487 0.26782243 -1.84450569 0.06510948 0.7596106
#> otu_3  3.532148   -0.456073305 0.26739635 -1.67993479 0.09297000 0.8134875
#> otu_4  2.729537    0.191041559 0.24401716  0.78618176 0.43176101 0.9444772
#> otu_5  2.673495   -0.077182700 0.25792249 -0.29876445 0.76511978 0.9650802
#> otu_6  3.267586    0.087067055 0.26802918  0.32445509 0.74559352 0.9650802
#> otu_7  4.082019   -0.330786588 0.26657283 -1.23314774 0.21752065 0.9444772
#> otu_8  3.467089   -0.362875098 0.26751575 -1.34401546 0.17894339 0.9079409
#> otu_9  3.088081   -0.249232240 0.25345437 -0.98124514 0.32647187 0.9444772
#> otu_10 2.709810   -0.211065497 0.25340686 -0.83171370 0.40557057 0.9444772
#> otu_11 2.686769    0.043266255 0.22141279  0.19555800 0.84495614 0.9650802
#> otu_12 3.422602    0.218668027 0.26788432  0.81243624 0.41654136 0.9444772
#> otu_13 2.381876    0.223783495 0.25733714  0.86609366 0.38643880 0.9444772
#> otu_14 3.223632    0.225181598 0.26820390  0.83545049 0.40346411 0.9444772
#> otu_15 3.326514   -0.511200242 0.26819128 -1.86788052 0.06177873 0.7596106
#> otu_16 3.061197    0.083506136 0.26444236  0.31524023 0.75257927 0.9650802
#> otu_17 4.129279   -0.069297689 0.26299717 -0.26380024 0.79193385 0.9650802
#> otu_18 4.953225    0.004134123 0.09549467  0.04328012 0.96547824 0.9879252
#> otu_19 2.439132   -0.366580938 0.25420878 -1.43014377 0.15267576 0.9079409
#> otu_20 2.758927    0.132955536 0.24789042  0.53910418 0.58981497 0.9650802
#> otu_21 3.272299   -0.140957759 0.26817542 -0.52543345 0.59928190 0.9650802
#> otu_22 2.608128    0.003734321 0.24540411  0.01513408 0.98792522 0.9879252
#> otu_23 3.021698    0.042751885 0.26111519  0.16395294 0.86976819 0.9650802
#> otu_24 2.652647    0.167098245 0.26121817  0.63743674 0.52384040 0.9650802
#> otu_25 2.136434    0.062900036 0.25848262  0.24250249 0.80839083 0.9650802
#> otu_26 3.035776    0.047768583 0.26141202  0.18309084 0.85472674 0.9650802
#> otu_27 3.230226    0.175193461 0.26812698  0.65088965 0.51511772 0.9650802
#> otu_28 3.344766   -0.117555247 0.25636633 -0.45957426 0.64582184 0.9650802
#> otu_29 1.964100    0.530107061 0.26020283  1.94016453 0.05235970 0.7596106
#> otu_30 3.456470   -0.360696897 0.26755941 -1.33588101 0.18158817 0.9079409
#> otu_31 2.413772    0.055798248 0.25822438  0.21628368 0.82876663 0.9650802
#> otu_32 3.184141   -0.217883070 0.26792220 -0.81057166 0.41761170 0.9444772
#> otu_33 2.785745    0.037330312 0.25221625  0.14797936 0.88235906 0.9650802
#> otu_34 2.152825   -0.250857707 0.25173542 -0.99522773 0.31962552 0.9444772
#> otu_35 3.704353    0.077417639 0.26628982  0.29047742 0.77145102 0.9650802
result$normalised_count  # Normalized count data
#>        Sample_1  Sample_2  Sample_3 Sample_4 Sample_5  Sample_6 Sample_7
#> otu_1    2.2906 0.0000000 3.9414682 4.009005 2.216844 1.8679453 3.239763
#> otu_2    5.7265 3.7699446 3.9414682 6.013508 3.325267 2.8019179 5.399606
#> otu_3    4.5812 5.6549168 1.9707341 4.009005 6.650533 3.7358906 1.079921
#> otu_4    3.4359 1.8849723 0.9853671 1.002251 3.325267 4.6698632 2.159842
#> otu_5    3.4359 2.8274584 3.9414682 2.004503 2.216844 1.8679453 1.079921
#> otu_6    1.1453 2.8274584 3.9414682 4.009005 3.325267 3.7358906 4.319685
#> otu_7    4.5812 3.7699446 2.9561012 8.018011 5.542111 3.7358906 0.000000
#> otu_8    2.2906 5.6549168 2.9561012 5.011257 5.542111 1.8679453 0.000000
#> otu_9    4.5812 3.7699446 3.9414682 7.015760 0.000000 5.6038358 2.159842
#> otu_10   4.5812 4.7124307 3.9414682 1.002251 2.216844 4.6698632 2.159842
#> otu_11   1.1453 0.0000000 1.9707341 5.011257 4.433689 1.8679453 1.079921
#> otu_12   1.1453 0.0000000 2.9561012 5.011257 5.542111 6.5378085 5.399606
#> otu_13   2.2906 0.9424861 0.9853671 3.006754 2.216844 3.7358906 3.239763
#> otu_14   4.5812 1.8849723 1.9707341 2.004503 3.325267 2.8019179 3.239763
#> otu_15   6.8718 1.8849723 5.9122023 2.004503 5.542111 0.9339726 3.239763
#> otu_16   1.1453 2.8274584 2.9561012 6.013508 1.108422 1.8679453 3.239763
#> otu_17   8.0171 3.7699446 4.9268353 3.006754 2.216844 9.3397264 2.159842
#> otu_18   4.5812 6.5974030 2.9561012 4.009005 5.542111 7.4717811 4.319685
#> otu_19   5.7265 3.7699446 2.9561012 2.004503 2.216844 0.0000000 2.159842
#> otu_20   2.2906 0.0000000 1.9707341 2.004503 5.542111 1.8679453 4.319685
#> otu_21   3.4359 2.8274584 1.9707341 2.004503 7.758955 0.9339726 5.399606
#> otu_22   2.2906 4.7124307 2.9561012 3.006754 0.000000 2.8019179 2.159842
#> otu_23   2.2906 2.8274584 2.9561012 3.006754 3.325267 0.9339726 6.479527
#> otu_24   2.2906 1.8849723 1.9707341 3.006754 2.216844 3.7358906 2.159842
#> otu_25   2.2906 4.7124307 1.9707341 1.002251 0.000000 2.8019179 2.159842
#> otu_26   2.2906 2.8274584 1.9707341 5.011257 2.216844 2.8019179 6.479527
#> otu_27   2.2906 3.7699446 2.9561012 2.004503 3.325267 4.6698632 1.079921
#> otu_28   0.0000 2.8274584 5.9122023 2.004503 7.758955 2.8019179 6.479527
#> otu_29   1.1453 0.9424861 1.9707341 0.000000 1.108422 1.8679453 3.239763
#> otu_30   2.2906 7.5398891 3.9414682 5.011257 2.216844 2.8019179 2.159842
#> otu_31   2.2906 1.8849723 2.9561012 2.004503 2.216844 2.8019179 4.319685
#> otu_32   0.0000 4.7124307 3.9414682 2.004503 7.758955 3.7358906 2.159842
#> otu_33   3.4359 1.8849723 3.9414682 3.006754 1.108422 4.6698632 1.079921
#> otu_34   4.5812 3.7699446 1.9707341 2.004503 1.108422 0.9339726 3.239763
#> otu_35   4.5812 5.6549168 1.9707341 1.002251 4.433689 5.6038358 5.399606
#>         Sample_8  Sample_9 Sample_10
#> otu_1  1.9598548 2.6069276 2.9561012
#> otu_2  1.9598548 1.7379517 0.0000000
#> otu_3  1.9598548 1.7379517 3.9414682
#> otu_4  3.9197095 0.0000000 5.9122023
#> otu_5  4.8996369 3.4759034 0.9853671
#> otu_6  2.9397822 3.4759034 2.9561012
#> otu_7  1.9598548 4.3448793 5.9122023
#> otu_8  1.9598548 3.4759034 5.9122023
#> otu_9  2.9397822 0.8689759 0.0000000
#> otu_10 1.9598548 0.8689759 0.9853671
#> otu_11 0.0000000 3.4759034 7.8829364
#> otu_12 2.9397822 1.7379517 2.9561012
#> otu_13 2.9397822 3.4759034 0.9853671
#> otu_14 5.8795643 2.6069276 3.9414682
#> otu_15 3.9197095 0.0000000 2.9561012
#> otu_16 3.9197095 2.6069276 4.9268353
#> otu_17 4.8996369 0.0000000 2.9561012
#> otu_18 4.8996369 5.2138551 3.9414682
#> otu_19 0.9799274 2.6069276 1.9707341
#> otu_20 4.8996369 1.7379517 2.9561012
#> otu_21 1.9598548 3.4759034 2.9561012
#> otu_22 2.9397822 5.2138551 0.0000000
#> otu_23 0.9799274 3.4759034 3.9414682
#> otu_24 1.9598548 4.3448793 2.9561012
#> otu_25 0.9799274 3.4759034 1.9707341
#> otu_26 3.9197095 0.8689759 1.9707341
#> otu_27 3.9197095 4.3448793 3.9414682
#> otu_28 2.9397822 1.7379517 0.9853671
#> otu_29 3.9197095 3.4759034 1.9707341
#> otu_30 5.8795643 1.7379517 0.9853671
#> otu_31 2.9397822 1.7379517 0.9853671
#> otu_32 0.9799274 2.6069276 3.9414682
#> otu_33 3.9197095 0.8689759 3.9414682
#> otu_34 3.9197095 0.0000000 0.0000000
#> otu_35 0.9799274 3.4759034 3.9414682
# }
```
