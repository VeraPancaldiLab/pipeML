# Ensure the censored package is loaded

Internal helper that loads the censored namespace, which registers the
"censored regression" engines used by parsnip. Loading the namespace is
sufficient; the package does not need to be attached.

## Usage

``` r
ensure_censored()
```
