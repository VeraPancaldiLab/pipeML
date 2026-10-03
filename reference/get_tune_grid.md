# Get Hyperparameter Grid for a Classification Method

Internal helper that returns the hyperparameter grid evaluated for each
classification method during custom-fold cross-validation:
[`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
passes it to
[`compute_custom_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_custom_k_fold_CV.md),
which trains every row in every fold. The grids are deliberately small
so that the full grid can be evaluated in every fold. The standard path
(no fold construction function) does not use these grids: it uses
caret's default grids, except for lasso and ridge, which use the same
`lambda` values as here.

## Usage

``` r
get_tune_grid(method, train_data)
```

## Arguments

- method:

  Character. One of `"glmnet"`, `"lasso"`, `"ridge"`, `"rf"`,
  `"svmRadial"`, `"treebag"`, `"C5.0"`, `"knn"`, `"rpart"`,
  `"svmLinear"`, or `"xgbTree"`.

- train_data:

  Data frame with the training features and the `target` column (the
  only non-feature column). Only used by `"rf"`, to scale `mtry` to the
  number of features (`ncol - 1`), and by `"svmRadial"`, to estimate
  `sigma`.
  [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
  passes the fold training data with the fewest features after
  preprocessing, the same in all folds, so that the grid is the same in
  all folds.

## Value

A data frame with one row per hyperparameter combination, suitable as
the `tuneGrid` argument of
[`caret::train()`](https://rdrr.io/pkg/caret/man/train.html):

- glmnet:

  `alpha` of 0 or 1, times 20 `lambda` values from 0.001 to 1, evenly
  spaced on the log scale (40 rows).

- lasso:

  `alpha = 1`, with the same 20 `lambda` values.

- ridge:

  `alpha = 0`, with the same 20 `lambda` values.

- rf:

  Up to 3 `mtry` values spanning 20-90 percent of the number of features
  (with 1 or 2 features, the values include 0, which randomForest resets
  to 1).

- svmRadial:

  Up to 3 `sigma` values estimated from the features (the 90, 50 and 10
  percent quantiles of \\1/\\x - x'\\^2\\ over all pairs of samples, on
  scaled features, as
  [`kernlab::sigest()`](https://rdrr.io/pkg/kernlab/man/sigest.html)
  does on a random subset), times `C` of 0.5, 1 or 2.

- treebag:

  A single placeholder row (`parameter = "none"`); no tunable
  parameters.

- C5.0:

  `trials` of 1, 5 or 10, times `winnow` TRUE or FALSE, with
  `model = "tree"`.

- knn:

  `k` of 3, 5, 7, 9 or 11 (odd values avoid ties).

- rpart:

  10 `cp` values from 0.001 to 0.1.

- svmLinear:

  `C` of 0.25, 0.5, 1, 2 or 4.

- xgbTree:

  `nrounds` of 100, 300 or 500, times `max_depth` of 3, 6 or 9, times
  `eta` of 0.01, 0.1 or 0.3, with `gamma = 0`, `colsample_bytree = 0.8`,
  `min_child_weight = 1` and `subsample = 0.8` fixed (27 rows).

## Details

The grids are deterministic (no random numbers are used): the same
`train_data` always gives the same grid. An error is raised for
unsupported methods.

## See also

[`compute_custom_k_fold_CV`](https://verapancaldilab.github.io/pipeML/reference/compute_custom_k_fold_CV.md),
[`get_default_hyperparams`](https://verapancaldilab.github.io/pipeML/reference/get_default_hyperparams.md)
for the survival-model equivalent.
