# Preprocess Features for Machine Learning

This function preprocesses a set of features by removing features with
near-zero variance and highly correlated features. It ensures that the
resulting feature set is more suitable for machine learning models. It
receives the features only: the callers set the outcome columns aside
and add them back, so the outcome is not used.

## Usage

``` r
preprocess_features(data, cor_thresh = 0.9)
```

## Arguments

- data:

  A data frame with the numeric features only (no outcome column).

- cor_thresh:

  A numeric value between 0 and 1 specifying the correlation threshold
  for removing highly correlated features. Default is `0.9`.

## Value

A data frame with the features that are kept (same rows).

## Details

The preprocessing steps include:

1.  Removing near-zero variance features (using
    [`caret::nearZeroVar`](https://rdrr.io/pkg/caret/man/nearZeroVar.html)).

2.  Removing highly correlated features: for each pair with an absolute
    correlation above the threshold, one of the two is removed (using
    [`caret::findCorrelation`](https://rdrr.io/pkg/caret/man/findCorrelation.html)).

The function stops if a feature is not numeric, or if no feature is left
after preprocessing.

## Examples

``` r
if (FALSE) { # \dontrun{
library(caret)
library(dplyr)

set.seed(123)
df <- data.frame(
  feature1 = c(1, 1, 1, 1, 1),             # constant
  feature2 = c(1, 2, 3, 4, 5),             # numeric
  feature3 = c(1, 2, 3, 4, 5) * 2          # highly correlated with feature2
)

clean_df <- preprocess_features(df, cor_thresh = 0.9)
} # }
```
