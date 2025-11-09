# Extract ML models from a directory based on specific AUC score

This function searches a directory for machine learning models and
filters them based on a specified AUC threshold for either the ROC or
Precision-Recall curves. It returns a list of model file names that meet
the specified AUC criteria.

## Usage

``` r
find.ML.models(folder_path, metric, AUC)
```

## Arguments

- folder_path:

  A character string specifying the directory path where the machine
  learning models are stored.

- metric:

  A character string indicating which AUC metric to use. Choose either
  "ROC" or "PRC".

- AUC:

  A numeric value representing the minimum acceptable AUC score for the
  models.

## Value

A character vector with the file paths of the ML models that meet the
AUC criteria.
