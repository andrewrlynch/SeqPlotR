# Coerce Hi-C input into a sparse contact `GRanges`

Accepts either:

- a `GRanges` whose mcols already carry `i_start`, `i_end`, `j_start`,
  `j_end`, `score`;

- a numeric matrix / data.frame whose `rownames` and `colnames` name bin
  positions as `"chr:start-end"` or integer bin indices.

## Usage

``` r
.hic_to_granges(data, windows)
```

## Arguments

- data:

  Input as described above.

- windows:

  A `GRanges` giving the genomic view (used to resolve integer-indexed
  matrix rownames / colnames).

## Value

A `GRanges` with `i_start`, `i_end`, `j_start`, `j_end`, `score` mcols.
