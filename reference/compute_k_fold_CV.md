# Perform repeated stratified k-fold cross-validation for model training and tuning

Internal function that performs repeated stratified k-fold
cross-validation to train and tune hyperparameters across multiple
machine learning models. Optionally, it can perform model stacking and
Boruta-based feature selection. Model performance is evaluated using
user-specified metrics such as Accuracy, AUROC, or AUPRC.

## Usage

``` r
compute_k_fold_CV(
  train_data,
  k_folds,
  n_rep,
  stacking = FALSE,
  metric = "Accuracy",
  file_name = NULL,
  LODO = FALSE,
  ncores = NULL,
  return = FALSE,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL
)
```

## Arguments

- train_data:

  A data frame containing predictor variables and a column named
  `target` corresponding to the response variable.

- k_folds:

  Integer. Number of folds used for k-fold cross-validation. Default is
  5.

- n_rep:

  Integer. Number of repetitions of the k-fold cross-validation. Default
  is 100.

- stacking:

  Logical. Whether to perform model stacking. Default is `FALSE`.

- metric:

  Character. Performance metric used for hyperparameter tuning and model
  evaluation. Supported values include `"Accuracy"`, `"AUROC"`, and
  `"AUPRC"`.

- file_name:

  Character. File name used when saving output plots in the `Results/`
  directory.

- LODO:

  Logical. If `TRUE`, performs Leave-One-Dataset-Out (LODO)
  cross-validation by stratifying folds based on cohort membership.

- ncores:

  Integer. Number of cores used for parallel computation. If `NULL`,
  `parallel::detectCores() - 1` will be used.

- return:

  Logical. Whether to return the results and generated plots.

- fold_construction_fun:

  Function used to construct cross-validation folds. The function must
  accept a `bestune` argument, which is used internally to inject
  optimized parameters after hyperparameter tuning. If `bestune = NULL`,
  the function explores a parameter grid across folds (parallelized with
  `foreach`). If `bestune` is provided, the optimized parameters are
  applied to rebuild features on the full training data.

- fold_construction_args_fixed:

  List of arguments passed to `fold_construction_fun` that remain fixed
  during both cross-validation and final training.

- fold_construction_args_tunable:

  List of arguments passed to `fold_construction_fun` that define
  hyperparameters to be tuned during cross-validation. Each element
  should contain candidate values.

## Value

A list containing:

- Features used during training

- The selected machine learning model

- All trained machine learning models

If `stacking = TRUE`, the list will also include:

- Base models

- Meta-learner

- Matrix of weighted feature importance (see
  [`calculate_feature_importance_stacking()`](https://verapancaldilab.github.io/pipeML/reference/calculate_feature_importance_stacking.md))
