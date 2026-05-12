# Build layout bounds for the row/column positional mode

Iterates `rows` (a list of lists of `SeqTrackR6`), assigning each row a
y-band and each track within a row its proportional x-band. Within each
track, the genomic windows are split into per-window panels.

## Usage

``` r
.build_positional_layout(rows, aesthetics)
```

## Arguments

- rows:

  A list of lists of `SeqTrackR6` objects.

- aesthetics:

  Merged plot-level aesthetics (already containing `margins`,
  `trackGaps`, `windowGaps`).

## Value

A list with `panelBounds` (integer-indexed) and `trackBounds`.
