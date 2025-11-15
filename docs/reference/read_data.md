# Extract specified data from a list of datasets

This function extracts a specific component (data) from a list of
datasets. The component to extract is specified by the `extract_name`
parameter, and the function returns a list containing the extracted data
from each dataset.

## Usage

``` r
read_data(dataset_list, extract_name)
```

## Arguments

- dataset_list:

  A list of datasets from which data will be extracted. Each element of
  the list is assumed to be a dataset (typically a list or dataframe).

- extract_name:

  A string representing the name of the component or column to be
  extracted from each dataset in the `dataset_list`. The function looks
  for this name within each dataset.

## Value

A list containing the extracted data. Each element corresponds to the
extracted component from the datasets in `dataset_list`. The names of
the list elements are taken from the names of `dataset_list`.

## Examples

``` r
# Example dataset list
dataset1 <- list(countdata = matrix(1:9, nrow = 3), metadata = data.frame(id = 1:3))
dataset2 <- list(countdata = matrix(10:18, nrow = 3), metadata = data.frame(id = 4:6))
dataset_list <- list(dataset1 = dataset1, dataset2 = dataset2)

# Extract 'countdata' from each dataset in the list
result <- read_data(dataset_list, "countdata")
print(result)
#> $dataset1
#>      [,1] [,2] [,3]
#> [1,]    1    4    7
#> [2,]    2    5    8
#> [3,]    3    6    9
#> 
#> $dataset2
#>      [,1] [,2] [,3]
#> [1,]   10   13   16
#> [2,]   11   14   17
#> [3,]   12   15   18
#> 
```
