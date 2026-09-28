# Ensure the caret package is attached

Internal helper that attaches caret to the search path if it is not
already attached. caret's model code is evaluated from the global
environment, and the `"knn"` model calls `knn3()` without a namespace
prefix. When caret is only loaded (e.g. `pipeML::` calls without
[`library(caret)`](https://github.com/topepo/caret/)), fitting a knn
model with `trainControl(method = "none")` fails with
`could not find function "knn3"`.

## Usage

``` r
ensure_caret()
```
