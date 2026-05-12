# Compute a y-range from a scale object (helper shared by y1 and y2)

Compute a y-range from a scale object (helper shared by y1 and y2)

## Usage

``` r
.scale_to_yrange(sy, y_windows = NULL)
```

## Arguments

- sy:

  A `SeqPositionScale` or `NULL`.

- y_windows:

  Optional `GRanges` — used when the scale is genomic.

## Value

Numeric length-2 vector, or `NULL` when nothing can be derived.
