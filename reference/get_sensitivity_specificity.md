# Calculate Sensitivity and Specificity Values

This function calculates sensitivity (recall), specificity, and other
related metrics (accuracy, precision, recall, F1 score, MCC) at each
probability threshold, from the predicted probabilities and the true
class labels.

## Usage

``` r
get_sensitivity_specificity(predictions, observed, ml.model)
```

## Arguments

- predictions:

  A data frame with a column `yes`: the predicted probability of the
  positive class of each sample.

- observed:

  A vector of true class labels (`"yes"` / `"no"`), in the same order as
  the rows of `predictions`.

- ml.model:

  Character. Name of the model, stored in the column `model` of the
  output.

## Value

A data frame with one row per sample, sorted by decreasing predicted
probability. Each row gives the metrics obtained when that sample and
all the samples above it are predicted as positive: `yes` (the
probability used as threshold), `model`, `Sensitivity`, `Specificity`,
`fpr`, `Accuracy`, `Precision`, `Recall`, `F1` and `MCC`.
