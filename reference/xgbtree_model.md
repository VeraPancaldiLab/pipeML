# caret model definition for xgbTree compatible with xgboost \>= 3

Internal helper returning a copy of caret's `"xgbTree"` model definition
that works with every xgboost version. Since xgboost 3, a booster is an
ALTREP object to which caret cannot add fields (`modelFit$xNames <- ...`
fails with "ALTLIST classes must provide a Set_elt method"), passing
`objective` outside `params` is deprecated and `ntreelimit` was removed.
Here the booster is stored in `modelFit$booster`, the objective goes in
`params` and sub-models are predicted with `iterationrange`.

## Usage

``` r
xgbtree_model()
```

## Value

A caret model definition (list) to pass as `method` to
[`caret::train()`](https://rdrr.io/pkg/caret/man/train.html).
