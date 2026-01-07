# Preprocess Features for Machine Learning

This function preprocesses a dataset by removing features with near-zero
variance, highly correlated features, and features that are constant
within any class of the target variable. It ensures that the resulting
feature set is more suitable for machine learning models.

## Usage

``` r
preprocess_features(
  data,
  target_col = NULL,
  time_var = NULL,
  event_var = NULL,
  cor_thresh = 0.9
)
```

## Arguments

- data:

  A data frame containing predictor features and the target variable.

- target_col:

  A character string specifying the name of the target column in `data`.
  Default is `"target"`.

- time_var:

  A character string specifying the name of the time-to-event column in
  `data`. Used only for survival analysis.

- event_var:

  A character string specifying the name of the event indicator column
  in `data` (e.g., 0/1). Used only for survival analysis.

- cor_thresh:

  A numeric value between 0 and 1 specifying the correlation threshold
  for removing highly correlated features. Default is `0.9`.

## Value

A data frame with the preprocessed features and the target column.

## Details

The preprocessing steps include:

1.  Removing near-zero variance features (using
    [`caret::nearZeroVar`](https://rdrr.io/pkg/caret/man/nearZeroVar.html)).

2.  Removing highly correlated features above the specified threshold
    (using
    [`caret::findCorrelation`](https://rdrr.io/pkg/caret/man/findCorrelation.html)).

3.  Removing features that are constant within any class of the target
    variable (i.e., provide no discriminatory power across classes).

After preprocessing, the target column is re-attached to the dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
library(caret)
library(dplyr)

set.seed(123)
df <- data.frame(
  feature1 = c(1, 1, 1, 1, 1),             # constant
  feature2 = c(1, 2, 3, 4, 5),             # numeric
  feature3 = c(1, 2, 3, 4, 5) * 2,         # highly correlated with feature2
  target   = c("A", "A", "B", "B", "B")
)

clean_df <- preprocess_features(df, target_col = "target", cor_thresh = 0.9)
} # }
```
