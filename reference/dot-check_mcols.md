# Check that required mcols columns exist on a GRanges

Internal helper. Given a GRanges `gr` and a character vector of column
names `cols`, errors if any are missing from `mcols(gr)`.

## Usage

``` r
.check_mcols(gr, cols)
```

## Arguments

- gr:

  A `GRanges` object.

- cols:

  Character vector of required column names in `mcols(gr)`.

## Value

Invisibly `TRUE` if the check passes; otherwise stops.
