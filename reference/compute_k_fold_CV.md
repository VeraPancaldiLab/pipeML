# Perform repeated stratified k-fold cross-validation for model training and tuning

Internal function that performs repeated stratified k-fold
cross-validation to train and tune hyperparameters across multiple
machine learning models. Hyperparameters are tuned and the best model is
selected with the user-specified metric (AUROC or AUPRC).

## Usage

``` r
compute_k_fold_CV(
  train_data,
  k_folds,
  n_rep,
  metric = "AUROC",
  file_name = NULL,
  LODO = FALSE,
  ncores = NULL,
  return = FALSE,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL,
  seed = 123,
  preprocess = TRUE
)
```

## Arguments

- train_data:

  A data frame containing predictor variables and a column named
  `target` corresponding to the response variable.

- k_folds:

  Integer. Number of folds used for k-fold cross-validation.

- n_rep:

  Integer. Number of repetitions of the k-fold cross-validation.

- metric:

  Character. Performance metric used for hyperparameter tuning and model
  selection: `"AUROC"` (default) or `"AUPRC"`.

- file_name:

  Character. File name used when saving output plots in the `Results/`
  directory.

- LODO:

  Logical. If `TRUE`, the folds are stratified by cohort and by `target`
  (for Leave-One-Dataset-Out analyses). `train_data` must then contain a
  column named `dataset` with the cohort of each sample; it is removed
  before training.

- ncores:

  Integer. Number of cores used for parallel computation. If `NULL`
  (default), the computation is sequential. Not used when
  `fold_construction_fun` is provided (the models are trained
  sequentially).

- return:

  Logical. Whether to save the cross-validation performance plots in the
  `Results/` directory.

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

- seed:

  Integer. Random seed set before the folds are drawn, so fold
  assignment and model fitting are reproducible. `NULL` leaves the
  random number generator untouched.

- preprocess:

  Logical. If `TRUE` (default), near-zero variance and highly correlated
  features are removed with
  [`preprocess_features()`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md):
  once on all training samples before the cross-validation, or, with
  `fold_construction_fun`, on the training part of each fold and on the
  final training set.

## Value

A list containing:

- `Model`: the selected machine learning model, trained on all training
  samples

- `ML_Models`: all trained machine learning models (models that predict
  the same value for all training samples are excluded)

- `AUROC_median`, `AUPRC_median`: median and MAD of the cross-validation
  AUROC and AUPRC of each model

- `Custom_output`: output of `fold_construction_fun` for the selected
  model (only with a custom function)
