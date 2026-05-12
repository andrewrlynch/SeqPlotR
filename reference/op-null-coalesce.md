# Null-coalescing operator

Returns `e1` if it is not `NULL`, otherwise returns `e2`.

## Usage

``` r
e1 %||% e2
```

## Arguments

- e1:

  Left-hand side value.

- e2:

  Right-hand side fallback value.

## Value

`e1` if not `NULL`, else `e2`.

## Examples

``` r
NULL %||% "fallback"
#> [1] "fallback"
"a" %||% "b"
#> [1] "a"
```
