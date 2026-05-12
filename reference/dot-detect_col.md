# Detect a column name from a list of candidates

Detect a column name from a list of candidates

## Usage

``` r
.detect_col(gr, hint, candidates, predicate, what)
```

## Arguments

- gr:

  A `GRanges` object.

- hint:

  Optional explicit column name.

- candidates:

  Character vector of candidate names, in preference order.

- predicate:

  A function `function(x)` accepting an mcols column and returning
  `TRUE` when the column is usable.

- what:

  Label for the error/warning message (e.g. `"CN state"`).

## Value

The detected column name.
