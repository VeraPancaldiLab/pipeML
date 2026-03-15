# DEPRECATED: Perform Platt Scaling for Probability Calibration

Fits a logistic regression model to calibrate predicted probabilities.

## Usage

``` r
compute_platt.scaling(obs, yes)
```

## Arguments

- obs:

  Vector of observed binary outcomes.

- yes:

  Vector of predicted probabilities for the positive class.

## Value

Numeric vector of calibrated probabilities.
