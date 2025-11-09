# Train model with optimized hyperparameters for classification tasks

This function wraps cross-validation, hyperparameter optimization, and
final training into a single workflow. It identifies the best
hyperparameters using a custom cross-validation function, reconstructs
the training set, preprocesses features, and retrains the model on the
complete training data with the selected hyperparameters.

## Usage

``` r
wrapper_train_best_hyperparams_classification(
  train_data,
  optimized,
  ml_method,
  fold_construction_fun,
  fold_construction_args_fixed
)
```

## Arguments

- train_data:

  A data frame containing the full training dataset, including
  predictors and the target variable.

- ml_method:

  A character string specifying the machine learning method to be passed
  to [`caret::train`](https://rdrr.io/pkg/caret/man/train.html).

- fold_construction_fun:

  A function used to (re)construct training data partitions given the
  best hyperparameters.

- fold_construction_args_fixed:

  A named list of additional fixed arguments to pass to
  `fold_construction_fun`.

- fold_data:

  A list or object containing pre-constructed folds for
  cross-validation, typically created by `fold_construction_fun`.

- tuneGrid:

  (optional) A data frame of hyperparameter values to evaluate. If
  `NULL`, defaults are used.

- ncores:

  (optional) Integer specifying the number of cores for parallel
  processing during cross-validation. If `NULL`, defaults to serial
  execution.

## Value

A list with the following components:

- `Model`:

  A trained `caret` model object with results, predictions, and
  resampling info attached.

- `training_set`:

  The final preprocessed training dataset.

- `custom_output`:

  Additional output from `fold_construction_fun`.

## Details

The workflow proceeds in the following steps:

1.  Runs cross-validation using `compute_custom_k_fold_CV` to identify
    the best hyperparameter set.

2.  Reconstructs the training set using `fold_construction_fun` and the
    selected hyperparameters.

3.  Preprocesses the training features by removing near-zero variance,
    highly correlated, and constant-within-class features (via
    [`preprocess_features`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md)).

4.  Retrains the model on the complete training data with the optimized
    hyperparameters.

The returned object mimics a `caret` model object but includes
additional elements derived from the custom cross-validation.

## See also

[`compute_custom_k_fold_CV`](https://verapancaldilab.github.io/pipeML/reference/compute_custom_k_fold_CV.md),
[`preprocess_features`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(caret)

# Example placeholders
train_data <- your_training_data
fold_data  <- your_prepared_folds

result <- wrapper_train_best_hyperparams_classification(
  train_data = train_data,
  fold_data  = fold_data,
  ml_method  = "rf",
  fold_construction_fun = your_fold_fun,
  fold_construction_args_fixed = list(arg1 = "value"),
  tuneGrid = expand.grid(mtry = 2:4),
  ncores = 4
)

result$Model
} # }
```
