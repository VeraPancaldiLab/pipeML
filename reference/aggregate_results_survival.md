# Aggregate Cross-Validation Results for Survival Models

This function aggregates and summarizes cross-validation results across
folds (and parameter configurations, if present) for multiple survival
models. It computes summary statistics such as the **median** and
**median absolute deviation (MAD)** of the Concordance Index (C-index)
per hyperparameter configuration and identifies the best-performing
configuration for each model.

## Usage

``` r
aggregate_results_survival(all_loaded)
```

## Arguments

- all_loaded:

  A nested list of model evaluation results from cross-validation:

  - If parameter tuning was used: `all_loaded[[fold]][[param]][[model]]`
    contains results for each model.

  - If no parameter tuning was used: `all_loaded[[fold]][[model]]`
    contains the results directly.

  Each result object must be a data frame containing columns:
  `"Resample"`, `"model"`, `"c_index"`, and any hyperparameter columns.

## Value

A list (one element per model), where each element contains:

- `Prediction_folds`:

  Data frame of all per-fold predictions and metrics.

- `Results_folds`:

  Summary table of median and MAD of C-index per configuration.

- `Besttune`:

  Data frame with the best hyperparameter configuration.

- `Resample_matrix`:

  Per-fold C-index values for the best configuration.

## Details

The function supports two types of data structures:

1.  Folds → Parameters → Models (when multiple parameter combinations
    are tested)

2.  Folds → Models (when there are no parameter combinations)

For each survival model, the function:

1.  Concatenates results across all folds (and parameter combinations if
    available).

2.  Groups by hyperparameter configuration to compute:

    - Median C-index (`c_index_median`)

    - MAD of C-index (`c_index_mad`)

3.  Selects the best hyperparameter configuration (highest median
    C-index).

4.  Extracts the per-fold results for the best configuration.

The output provides both a detailed and summarized view of model
performance, enabling comparison of tuned configurations and stability
across folds.

## See also

[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
[`compute_cv_CINDEX()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_CINDEX.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: aggregate cross-validation results after nested CV
aggregated <- aggregate_results_survival(all_loaded)

# Access the best hyperparameters for model 1
aggregated[[1]]$Besttune

# View per-fold C-index for the best configuration
aggregated[[1]]$Resample_matrix
} # }
```
