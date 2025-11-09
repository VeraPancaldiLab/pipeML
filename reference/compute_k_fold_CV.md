# Perform repeated stratified k-fold cross-validation for model training and tuning

This function performs repeated stratified k-fold cross-validation on a
dataset to train and tune hyperparameters for 13 machine learning
methods. Optionally, it can also perform model stacking and Boruta-based
feature selection. Performance is evaluated using user-specified metrics
such as Accuracy, AUROC, or AUPRC.

## Usage

``` r
compute_k_fold_CV(
  model,
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

- model:

  A data frame containing features and a target column named 'target'
  corresponding to the response variable to predict.

- k_folds:

  Integer. Number of folds for k-fold cross-validation. Default is 5.

- n_rep:

  Integer. Number of repetitions of the k-fold cross-validation. Default
  is 100.

- stacking:

  Logical. Whether to perform model stacking. Default is FALSE.

- metric:

  Character. Metric used for hyperparameter tuning and model evaluation.
  Supported values are "Accuracy", "AUROC", and "AUPRC".

- file_name:

  Character. File name used for saving output plots in the `Results/`
  directory.

- LODO:

  Logical. If TRUE, performs Leave-One-Dataset-Out (LODO)
  cross-validation by stratifying folds based on cohort membership.

- ncores:

  Integer. Number of cores to use for parallelization. If not given,
  detectCores() - 1 will be used.

- return:

  Logical. Whether to return the results and generated plots.

- fold_construction_fun:

  Function. A custom function used to construct the cross-validation
  folds. This function must accept a `bestune` argument, which is used
  internally to inject optimized parameters after hyperparameter tuning.
  If `bestune = NULL`, the function will explore a parameter grid across
  folds (parallelized with `foreach`); if `bestune` is provided, the
  optimized parameters will be applied to rebuild the features on the
  full training data.

- fold_construction_args_fixed:

  List. A list of arguments passed to `fold_construction_fun` that
  remain fixed during both cross-validation and final training.

- fold_construction_args_tunable:

  List. A list of arguments passed to `fold_construction_fun` that
  define the hyperparameters to be tuned during cross-validation. Each
  element should contain candidate values for tuning.

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
