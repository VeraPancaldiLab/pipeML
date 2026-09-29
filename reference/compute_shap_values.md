# Compute SHAP Values for Machine Learning Models

This function calculates SHAP (SHapley Additive exPlanations) values to
assess feature importance for the final model selected by
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md),
i.e. the model trained on all training samples with the tuned
hyperparameters. SHAP values are computed for every training sample, on
the features of the final model, so each feature has the same definition
for all samples (including features built by a custom fold construction
function, which are computed once on all training samples). It supports
both classification and survival tasks.

## Usage

``` r
compute_shap_values(model_trained, task_type = "classification", seed = 123)
```

## Arguments

- model_trained:

  The trained model returned as `$Model` by
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (e.g. `res$Model`). Everything else is taken from it:

  - Classification: the caret `train` object; features of the training
    samples from `$trainingData`; SHAP values explain the predicted
    probability of the positive class (`"yes"`).

  - Survival: the final fitted model in `$Model_object`; features of the
    training samples from `$trainingData` (the `time` and `event`
    columns are not used as features); SHAP values explain the predicted
    risk score.

- task_type:

  Character. Either `"classification"` (default) or `"survival"`.

- seed:

  Integer. Random seed for reproducible SHAP values (they are Monte
  Carlo estimates). Default: `123`. Use `NULL` to leave the random
  number generator untouched.

## Value

A data frame of SHAP values, with rows corresponding to training samples
(sample IDs as rownames) and columns to the features of the final model.
Values are in the units of the model output (probability of the positive
class for classification, risk score for survival). The average
prediction of the model on the training samples is stored in the
attribute `"baseline"`: for each sample, the baseline plus the sum of
its SHAP values equals its prediction. If the model predicts the same
value for all samples, a warning is issued and `NULL` is returned.

## Details

SHAP values are estimated with
[`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
(100 Monte Carlo simulations per feature), using the training samples
both as the samples to explain and as background data.

Because the final model was trained on the samples it explains, SHAP
values describe how the model uses each feature on its training data.
With flexible models and small datasets, a model that overfits can
assign importance to features that help fit the training samples but do
not generalize: compare the training performance with the
cross-validation performance before interpreting them.

## Examples

``` r
if (FALSE) { # \dontrun{
data(data_example_classification)
X <- data_example_classification[, setdiff(colnames(data_example_classification), "target")]
res <- compute_features.training.ML(features_train = X,
                                    target_var = data_example_classification$target,
                                    task_type = "classification",
                                    trait.positive = "1",
                                    k_folds = 5,
                                    n_rep = 2)

shap <- compute_shap_values(res$Model, task_type = "classification")
head(shap)

# Global importance: mean absolute SHAP value of each feature
sort(colMeans(abs(shap)), decreasing = TRUE)
} # }
```
