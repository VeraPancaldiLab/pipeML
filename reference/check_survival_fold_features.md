# Check the survival features built by a custom fold construction function

Internal helper. The fold construction function receives the survival
outcome in the `time` and `event` columns of `data`; the features it
returns must keep these two columns and must not contain a copy of them
(e.g. because the outcome was not removed before building the features).

## Usage

``` r
check_survival_fold_features(df, where)
```

## Arguments

- df:

  Data frame of features plus `time` and `event` (training data of a
  fold, or final training set).

- where:

  Character. Where `df` comes from, used in the error message.

## Value

Invisibly `TRUE`; stops with an informative error otherwise.
