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
# Example usage
set.seed(101)
nr = 10; nc =50
countdata <- matrix(rpois(500, 3), ncol = nc, nrow = nr)
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
#> -0.007959025 -0.445880472 -0.394760811  0.233294799 -0.091829484  0.058651776 
#>        otu_7        otu_8        otu_9       otu_10       otu_11       otu_12 
#> -0.308246138 -0.308896891 -0.251029922 -0.252289726  0.039529520  0.160924136 
#>       otu_13       otu_14       otu_15       otu_16       otu_17       otu_18 
#>  0.194703220  0.192716051 -0.464400160  0.068904985 -0.097712013  0.053908758 
#>       otu_19       otu_20       otu_21       otu_22       otu_23       otu_24 
#> -0.363959631  0.156101373 -0.144220223 -0.013117419  0.028005740  0.150500957 
#>       otu_25       otu_26       otu_27       otu_28       otu_29       otu_30 
#>  0.042978792  0.027651694  0.141706746 -0.177098228  0.456759743 -0.309509409 
#>       otu_31       otu_32       otu_33       otu_34       otu_35       otu_36 
#>  0.035671304 -0.228363278  0.030497385 -0.256186889  0.044228160 -0.083204924 
#>       otu_37       otu_38       otu_39       otu_40       otu_41       otu_42 
#>  0.251044749  0.271340396 -0.148152551 -0.417094381  0.003548958  0.101448484 
#>       otu_43       otu_44       otu_45       otu_46       otu_47       otu_48 
#> -0.023449549 -0.062952952 -0.208034179 -0.104467686 -0.019733425  0.030671416 
#>       otu_49       otu_50 
#>  0.273876722  0.149311832 
result$logmean  # Log mean counts
#> NULL
result$dispersion  # Dispersion estimates
#>  [1] 0.23726334 0.11818891 0.12779467 0.20192426 0.17595129 0.07919315
#>  [7] 0.11134903 0.14821130 0.42777852 0.19140022 0.85818539 0.12657805
#> [13] 0.24094148 0.07205231 0.09815225 0.07600581 0.14666306 0.04598340
#> [19] 0.25463325 0.16538272 0.08420186 0.23473744 0.06856281 0.17004895
#> [25] 0.22618425 0.06933439 0.07672390 0.12162285 0.14324368 0.14440827
#> [31] 0.23368195 0.08659100 0.11722574 0.24416527 0.13947440 0.28604738
#> [37] 0.08057311 0.11586475 0.07294106 0.07079641 0.13274784 0.07402237
#> [43] 0.14003320 0.11678770 0.13505417 0.06700914 0.06674076 0.10383276
#> [49] 0.21379294 0.17696335
result$deseq_estimate  # DESeq2 results
#>        baseMean log2FoldChange     lfcSE        stat     pvalue      padj
#> otu_1  2.523916   -0.007959025 0.2519496 -0.03164913 0.97475186 0.9892808
#> otu_2  3.436181   -0.445880472 0.2632859 -1.67359558 0.09421011 0.9892808
#> otu_3  3.542517   -0.394760811 0.2633018 -1.48732907 0.13692792 0.9892808
#> otu_4  2.702925    0.233294799 0.2553894  0.91069198 0.36245769 0.9892808
#> otu_5  2.694906   -0.091829484 0.2572253 -0.35677492 0.72126031 0.9892808
#> otu_6  3.287608    0.058651776 0.2642521  0.22185078 0.82443004 0.9892808
#> otu_7  4.150850   -0.308246138 0.2647080 -1.16052986 0.24583314 0.9892808
#> otu_8  3.555650   -0.308896891 0.2626516 -1.17228830 0.24108133 0.9892808
#> otu_9  3.090175   -0.251029922 0.2446089 -1.02344368 0.30609811 0.9892808
#> otu_10 2.675081   -0.252289726 0.2561637 -0.97892159 0.32761872 0.9892808
#> otu_11 2.749900    0.039529520 0.2178270  0.18227420 0.85536754 0.9892808
#> otu_12 3.390569    0.160924136 0.2630611  0.61090560 0.54126208 0.9892808
#> otu_13 2.356214    0.194703220 0.2500851  0.77559202 0.43798991 0.9892808
#> otu_14 3.205661    0.192716051 0.2641421  0.72725190 0.46707166 0.9892808
#> otu_15 3.301207   -0.464400160 0.2633714 -1.74236099 0.08144530 0.9892808
#> otu_16 3.130311    0.068904985 0.2638552  0.26125914 0.79389267 0.9892808
#> otu_17 4.059487   -0.097712013 0.2640602 -0.36948165 0.71176875 0.9892808
#> otu_18 4.936842    0.053908758 0.2616827  0.20581280 0.83693715 0.9892808
#> otu_19 2.433494   -0.363959631 0.2494307 -1.44767679 0.14770746 0.9892808
#> otu_20 2.736376    0.156101373 0.2577149  0.60598515 0.54452462 0.9892808
#> otu_21 3.243073   -0.144220223 0.2639479 -0.54652933 0.58470213 0.9892808
#> otu_22 2.633467   -0.013117419 0.2535946 -0.05154981 0.95888741 0.9892808
#> otu_23 3.018337    0.028005740 0.2637688  0.10623129 0.91539885 0.9892808
#> otu_24 2.658367    0.150500957 0.2572649  0.58326522 0.55971477 0.9892808
#> otu_25 2.140747    0.042978792 0.2488796  0.17205591 0.86339357 0.9892808
#> otu_26 3.019768    0.027651694 0.2637470  0.10490246 0.91645321 0.9892808
#> otu_27 3.255246    0.141706746 0.2642600  0.53474221 0.59282813 0.9892808
#> otu_28 3.335342   -0.177098228 0.2630580 -0.67289343 0.50101509 0.9892808
#> otu_29 1.956007    0.456759743 0.2498016  1.76909120 0.07687866 0.9892808
#> otu_30 3.517628   -0.309509409 0.2627584 -1.17253684 0.24098159 0.9892808
#> otu_31 2.388890    0.035671304 0.2509830  0.14209553 0.88700455 0.9892808
#> otu_32 3.218515   -0.228363278 0.2638465 -0.86312957 0.38806623 0.9892808
#> otu_33 2.797655    0.030497385 0.2609264  0.11688845 0.90694845 0.9892808
#> otu_34 2.127118   -0.256186889 0.2466555 -1.03492011 0.30070622 0.9892808
#> otu_35 3.643876    0.044228160 0.2634726  0.16759310 0.86690340 0.9892808
#> otu_36 2.505223   -0.083204924 0.2488969 -0.33428091 0.73816759 0.9892808
#> otu_37 2.873749    0.251044749 0.2627610  0.94886216 0.34269072 0.9892808
#> otu_38 3.371951    0.271340396 0.2635050  1.02123215 0.30714446 0.9892808
#> otu_39 3.037822   -0.148152551 0.2636697 -0.56179228 0.57425755 0.9892808
#> otu_40 3.034731   -0.417094381 0.2635909 -1.56179566 0.11833613 0.9892808
#> otu_41 3.759194    0.003548958 0.2639096  0.01343493 0.98928080 0.9892808
#> otu_42 3.121752    0.101448484 0.2641010  0.38304575 0.70168584 0.9892808
#> otu_43 3.885608   -0.023449549 0.2639198 -0.08883389 0.92921392 0.9892808
#> otu_44 3.571744   -0.062952952 0.2639640 -0.23829804 0.81164994 0.9892808
#> otu_45 3.560160   -0.208034179 0.2633832 -0.78692916 0.43132333 0.9892808
#> otu_46 2.933184   -0.104467686 0.2634991 -0.39666457 0.69161483 0.9892808
#> otu_47 3.075949   -0.019733425 0.2640999 -0.07469519 0.94045724 0.9892808
#> otu_48 2.803747    0.030671416 0.2615479  0.11730384 0.90661928 0.9892808
#> otu_49 2.045793    0.273876722 0.2476936  1.09437629 0.27378998 0.9892808
#> otu_50 2.686524    0.149311832 0.2569471  0.58024839 0.56174712 0.9892808
result$normalised_count  # Normalized count data
#> NULL
```
