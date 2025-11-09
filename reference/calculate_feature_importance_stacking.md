# Compute Weighted Feature Importance from Base Models and Meta-Learner for Stacking Models

This function computes the feature importance by weighing the feature
importances from multiple base models in a stacking ensemble, combined
with the meta-learner's model importance. The final importance score for
each feature is calculated by multiplying the base model's feature
importance with the meta-learner's weight for each base model.

## Usage

``` r
calculate_feature_importance_stacking(
  base_importance,
  base_models,
  meta_learner
)
```

## Arguments

- base_importance:

  A list where each element corresponds to a base model and contains a
  data frame with feature importances. Each data frame should have a
  column called `importance` (either for the positive class or overall,
  depending on the type of model).

- base_models:

  A character vector with the names of the base models whose feature
  importances are provided in `base_importance`.

- meta_learner:

  A `caret` object representing the trained meta-learner model. This
  model is used to obtain weights for each base model, based on their
  performance in the ensemble.

## Value

A data frame with two columns:

- features:

  The feature names.

- final_importance:

  The final weighted importance score for each feature, calculated by
  summing the weighted importances across all base models. Features are
  sorted in descending order of their final importance score.

## Details

The function extracts feature importance values from the base models,
then computes the weighted importance for each feature based on the
meta-learner's performance. The meta-learner's feature importance is
normalized to ensure the sum of the importances across all models is 1,
and it is used as the weight for each base model. The feature
importances from all base models are then aggregated and weighted by
their respective meta-learner importance scores.

## See also

[`varImp`](https://rdrr.io/pkg/caret/man/varImp.html)
