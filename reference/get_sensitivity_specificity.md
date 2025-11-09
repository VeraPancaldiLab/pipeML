# Calculate Sensitivity and Specificity Values

This function calculates sensitivity (recall), specificity, and other
related metrics (accuracy, precision, recall, F1 score, MCC) from
predicted and true class labels.

## Usage

``` r
get_sensitivity_specificity(predictions, observed, ml.model)
```

## Arguments

- predictions:

  A vector of predicted class labels or probabilities.

- observed:

  A vector of true class labels.

- ml.model:

  The trained machine learning model used to generate predictions.

## Value

A data frame containing sensitivity, specificity, precision, recall, F1
score, MCC, and other metrics.
