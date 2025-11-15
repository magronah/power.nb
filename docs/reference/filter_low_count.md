# Filter to remove low abundant taxa

Filter to retain only taxa with at least `abund_thresh` counts in at
least `sample_thresh` samples

## Usage

``` r
filter_low_count(
  countdata,
  metadata,
  abund_thresh = 5,
  sample_thresh = 3,
  sample_colname,
  group_colname
)
```

## Arguments

- countdata:

  otu table

- metadata:

  dataframe with 2 rows sample names and group names

- abund_thresh:

  minimum number of taxa abundance threshold

- sample_thresh:

  minimum number of sample threshold

- sample_colname:

  column names of the samples

- group_colname:

  column names of the groups or conditions

## Value

filtered otu count data
