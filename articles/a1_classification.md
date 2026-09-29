# Classification

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

## **Data**

`data_example_classification` is the Breast Cancer Wisconsin dataset: 9
cell features per sample and a `target` column (`1` = malignant, `0` =
benign).

``` r

data <- pipeML::data_example_classification
X <- data %>% dplyr::select(-target)
y <- data$target
```

Split the samples into a training and a test set:

``` r

set.seed(123)
train_idx <- caret::createDataPartition(y, p = 0.7, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]
```

## **Train models**

[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
trains and tunes 11 classification algorithms (bagged trees, random
forest, C5.0, elastic net, lasso, ridge, k-nearest neighbours, CART,
linear and radial SVM, and XGBoost) with repeated stratified k-fold
cross-validation, and selects the best one.

Key parameters:

- **trait.positive**: value of `target_var` considered as the positive
  class
- **metric**: metric used to tune and select the models, `"AUROC"`
  (default) or `"AUPRC"`
- **k_folds**, **n_rep**: number of folds and repetitions of the
  cross-validation
- **ncores**: number of cores used to run the folds in parallel (`NULL`
  runs sequentially)
- **seed**: random seed for the folds and the model fitting (default
  `123`), so results are reproducible
- **return**: whether to save the cross-validation performance plots in
  `Results/` (named with `file_name`)

``` r

res <- compute_features.training.ML(features_train = X_train,
                                    target_var = y_train,
                                    task_type = "classification",
                                    trait.positive = "1",
                                    metric = "AUROC",
                                    k_folds = 5,
                                    n_rep = 2,
                                    ncores = 2,
                                    seed = 123,
                                    file_name = "Example_classification",
                                    return = TRUE)
```

The selected model is a `caret` `train` object, trained on all training
samples with the tuned hyperparameters:

``` r

res$Model
res$Model$bestTune
```

All trained and tuned models:

``` r

names(res$ML_Models)
```

`res$AUROC_median` and `res$AUPRC_median` compare the cross-validation
performance of all the algorithms (median and MAD across resamples). The
performance of the selected model per resample is in
`res$Model$resample`:

``` r

res$AUROC_median
head(res$Model$resample)
```

![Figure 1. Cross-validation performance of the trained
models.](figures/AUROC_classification.png)

Figure 1. Cross-validation performance of the trained models.

## **Predict on test data**

[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
applies the selected model to the test set. The test data must contain
the features used for training:

``` r

pred <- compute_prediction(model = res$Model,
                           test_data = X_test,
                           target_var = y_test,
                           task_type = "classification",
                           trait.positive = "1",
                           file.name = "Example_classification",
                           return = TRUE)
```

AUROC and AUPRC on the test set. Each is a list with the `estimate` on
the full test set and the `lower` and `upper` bounds of its 95%
confidence interval, from 1000 bootstrap resamples of the test samples:

``` r

pred$AUC
```

Predicted probabilities of each class (`yes` is the positive class):

``` r

head(pred$Predictions)
```

Accuracy, sensitivity, specificity, precision, recall, F1 score and MCC
at each probability threshold:

``` r

head(pred$Metrics)
```

With `return = TRUE`, the ROC and precision-recall curves are saved in
`Results/` (named with `file.name`). The shaded areas are pointwise 95%
bootstrap confidence bands, also returned in `pred$Curve_bands`.

![Figure 2. ROC curve on the test set.](figures/ROC.png)

Figure 2. ROC curve on the test set.

![Figure 3. Precision-recall curve on the test set.](figures/PR.png)

Figure 3. Precision-recall curve on the test set.

The curves can also be drawn directly with
[`get_curves()`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md):

``` r

get_curves(data = pred$Metrics,
           color = "model",
           auc_roc = pred$AUC$AUROC,
           auc_prc = pred$AUC$AUPRC,
           roc_band = pred$Curve_bands$ROC,
           prc_band = pred$Curve_bands$PRC,
           file.name = "Example_classification_curves")
```

## **SHAP values**

[`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md)
explains which features drive the predictions of the selected model. It
only needs the trained model:

``` r

shap <- compute_shap_values(model_trained = res$Model,
                            task_type = "classification",
                            seed = 123)
head(shap)
```

See the **Interpreting models with SHAP values** tutorial for how to
read and plot them.

## **Training and prediction in one step**

When the test set is already prepared,
[`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
trains the models on the training set and applies the selected one to
the test set. The outcome is read from `coldata`, whose row names must
be the sample names of the feature tables; `trait` is the outcome
column:

``` r

res_onestep <- compute_features.ML(features_train = X_train,
                                   features_test = X_test,
                                   coldata = data,
                                   task_type = "classification",
                                   trait = "target",
                                   trait.positive = "1",
                                   metric = "AUROC",
                                   k_folds = 5,
                                   n_rep = 2,
                                   ncores = 2,
                                   file_name = "Example_onestep",
                                   return = FALSE)
```

`res_onestep$Model` is the output of
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md),
and the test set results are in `AUC`, `Metrics`, `Prediction` and
`Curve_bands`:

``` r

res_onestep$Model$Model
res_onestep$AUC
```
