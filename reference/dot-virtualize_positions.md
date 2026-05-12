# Convert genomic (seqname, position) pairs to virtual positions

Uses the map returned by
[`.build_virtual_map()`](http://andrewlynch.io/SeqPlotR/reference/dot-build_virtual_map.md).
Each position is resolved against the window whose seqname matches and
whose genomic range contains the position. Positions that match no
window get `NA_real_`.

## Usage

``` r
.virtualize_positions(seqnames, positions, vmap)
```

## Arguments

- seqnames:

  Character vector of seqnames.

- positions:

  Numeric vector of genomic positions, same length.

- vmap:

  A virtualization map.

## Value

Numeric vector of virtual positions, `NA` for unmatched.
