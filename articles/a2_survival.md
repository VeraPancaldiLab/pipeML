# Survival analysis

``` r

library(pipeML)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## **Survival outcomes**

`pipeML` also trains models that predict time-to-event outcomes. The
outcome has two components:

- **time**: follow-up or survival time
- **event**: whether the event occurred (1) or the observation was
  censored (0)

Survival models are built with `parsnip` and its `censored` extension.
`censored` must be installed, but you don’t need to load it: `pipeML`
loads it when a survival task runs.

``` r

install.packages("censored")  # only needed once
```

## **Data**

`data_example_survival` is the lung cancer dataset of the `survival`
package: survival `time` in days, event indicator `status` (1 = death, 0
= censored) and 8 covariates.

``` r

data <- pipeML::data_example_survival
X <- data %>% dplyr::select(-time, -status)
time <- data$time
event <- data$status
```

Split the samples into a training and a test set, stratified by event:

``` r

set.seed(123)
train_idx <- caret::createDataPartition(event, p = 0.7, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
time_train <- time[train_idx]
time_test  <- time[-train_idx]
event_train <- event[train_idx]
event_test  <- event[-train_idx]
```

## **Train models**

With `task_type = "survival"`,
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
trains and tunes six survival models (Cox proportional hazards,
elastic-net Cox, parametric accelerated failure time model, conditional
inference tree, bagged CART and oblique random survival forest) with
repeated k-fold cross-validation stratified by event. `time_var` and
`event_var` are the time and event of each training sample. The best
model is selected by the concordance index (C-index).

Survival models are slower to train than classification models: here we
use 3 folds and 1 repetition to keep the example fast. Use more folds
and repetitions (e.g. `k_folds = 5`, `n_rep = 5`) for real analyses.

``` r

res_survival <- compute_features.training.ML(features_train = X_train,
                                             task_type = "survival",
                                             time_var = time_train,
                                             event_var = event_train,
                                             k_folds = 3,
                                             n_rep = 1,
                                             ncores = 2,
                                             seed = 123,
                                             file_name = "Example_survival",
                                             return = TRUE)
```

All trained models and the name of the selected one:

``` r

names(res_survival$ML_Models)
res_survival$Model$model
```

The selected model trained on all training samples with the tuned
hyperparameters (a fitted `workflows` object), and those
hyperparameters:

``` r

res_survival$Model$Model_object
res_survival$Model$bestTune
```

Cross-validation performance: median C-index of the selected model, and
its C-index per resample:

``` r

res_survival$C_index_median
head(res_survival$Model$Resample_matrix)
```

With `return = TRUE`, the C-index of all models across resamples is
saved in `Results/` (named with `file_name`).

![Figure 1. Cross-validation C-index of the trained
models.](figures/cindex_survival.png)

Figure 1. Cross-validation C-index of the trained models.

## **Predict on test data**

``` r

pred_survival <- compute_prediction(model = res_survival$Model,
                                    test_data = X_test,
                                    task_type = "survival",
                                    time_var = time_test,
                                    event_var = event_test,
                                    file.name = "Example_survival",
                                    return = TRUE)
```

C-index on the test set, with its 95% confidence interval:

``` r

pred_survival$c_index
c(pred_survival$c_index_lower, pred_survival$c_index_upper)
```

Predicted risk scores of the test samples:

``` r

head(pred_survival$preds)
```

Depending on the model, survival models predict a risk score, a survival
time or a survival probability.
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
converts them into a risk score so that all models are read the same
way:

- risk score (**linear_pred**), e.g. Cox models, is used as is;
- predicted survival time (**time**), e.g. parametric models, is
  reversed: a longer survival means a lower risk;
- survival probability (**survival**) is also reversed: a higher
  probability of survival means a lower risk.

**Higher prediction values always correspond to higher predicted risk.**

## **Kaplan-Meier curves by risk group**

With `return = TRUE`,
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
splits the test samples into two groups at the median predicted risk and
saves their Kaplan-Meier curves in `Results/`, with the C-index and the
log-rank test p-value.

![Figure 2. Kaplan-Meier curves of the test samples by predicted risk
group.](figures/KM.png)

Figure 2. Kaplan-Meier curves of the test samples by predicted risk
group.

To use another number of groups, call
[`plot_survival_performance()`](https://verapancaldilab.github.io/pipeML/reference/plot_survival_performance.md)
with the observed outcome of the test samples and the output of
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md).
Groups are defined by quantiles of the predicted risk:

``` r

km <- plot_survival_performance(df_test = data.frame(time = time_test, event = event_test),
                                prediction = pred_survival,
                                n_groups = 3,
                                file_name = "Example_survival_3groups")
```

## **SHAP values**

[`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md)
works the same way as for classification, with `task_type = "survival"`.
SHAP values are in risk-score units: positive values push the prediction
towards a higher risk.

``` r

shap_survival <- compute_shap_values(model_trained = res_survival$Model,
                                     task_type = "survival",
                                     seed = 123)
head(shap_survival)
```

## **Training and prediction in one step**

[`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
trains and predicts in one step. For survival, `time_var` and
`event_var` are the names of the time and event columns of `coldata`,
whose row names must be the sample names of the feature tables:

``` r

res_onestep_survival <- compute_features.ML(features_train = X_train,
                                            features_test = X_test,
                                            coldata = data,
                                            task_type = "survival",
                                            time_var = "time",
                                            event_var = "status",
                                            k_folds = 3,
                                            n_rep = 1,
                                            ncores = 2,
                                            file_name = "Example_onestep_survival",
                                            return = FALSE)
```

`res_onestep_survival$Model` is the output of
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md),
`C_index` the C-index on the test set and `Prediction` the predicted
risk scores:

``` r

res_onestep_survival$C_index
head(res_onestep_survival$Prediction)
```
