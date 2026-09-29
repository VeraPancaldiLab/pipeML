# Train a caret model, with pipeML's xgbTree definition

Internal wrapper around
[`caret::train()`](https://rdrr.io/pkg/caret/man/train.html) that
replaces `method = "xgbTree"` with
[`xgbtree_model()`](https://verapancaldilab.github.io/pipeML/reference/xgbtree_model.md)
and keeps `"xgbTree"` as the `method` of the returned object.

## Usage

``` r
caret_train(..., method)
```

## Arguments

- ...:

  Arguments passed to
  [`caret::train()`](https://rdrr.io/pkg/caret/man/train.html).

- method:

  Character. caret method name.

## Value

A caret `train` object.
