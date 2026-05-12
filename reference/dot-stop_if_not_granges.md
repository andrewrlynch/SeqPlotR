# Assert that an object is a GRanges

Internal helper. Errors if `x` does not inherit from `"GRanges"`.

## Usage

``` r
.stop_if_not_granges(x, arg_name)
```

## Arguments

- x:

  Object to check.

- arg_name:

  Name of the argument, used in the error message.

## Value

Invisibly `TRUE` if the check passes; otherwise stops.
