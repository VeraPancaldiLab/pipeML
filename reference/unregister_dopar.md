# Reset foreach backend to sequential

Internal helper that unregisters any active parallel backend used by
foreach and restores the sequential backend via
[`foreach::registerDoSEQ()`](https://rdrr.io/pkg/foreach/man/registerDoSEQ.html).
This is used to ensure that parallel clusters created during model
training are properly released and that subsequent operations run
sequentially.

## Usage

``` r
unregister_dopar()
```
