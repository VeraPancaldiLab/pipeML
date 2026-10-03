# Evaluate a Precision-Recall Curve on a Fixed Recall Grid (Internal)

Returns, for each grid value, the precision of the precision-recall
curve at that recall, with straight lines between consecutive points:
the curve as
[`get_curves()`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md)
draws it (`geom_line()`) and as
[`calculate_auprc()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auprc.md)
integrates it, starting at recall 0 with the precision of the first
point. At a grid value equal to the recall of a point, it is the
precision at the first threshold reaching that recall (so the final drop
to the prevalence at recall 1 is not in the band).

## Usage

``` r
prc_at_grid(recall, precision, grid)
```

## Arguments

- recall:

  Numeric vector of recall values, non-decreasing (as produced by
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md)).

- precision:

  Numeric vector of precision values, same length as `recall`.

- grid:

  Numeric vector of recall values in \[0, 1\].

## Value

Numeric vector of precisions, one per grid value (`NA` if the curve is
undefined, e.g. a resample without positives).
