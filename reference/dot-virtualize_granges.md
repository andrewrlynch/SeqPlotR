# Virtualize a `GRanges` according to a virtualization map

Returns a new `GRanges` with seqnames replaced by the combined sentinel,
and ranges replaced by virtual coordinates. Rows whose original
(seqname, position) pair does not fall in any window are silently
dropped. Optional companion mcols giving a *secondary* position pair
(e.g. `j_start`/`j_end` with an optional `j_chrom` column for Hi-C) are
also virtualized in place.

## Usage

``` r
.virtualize_granges(
  gr,
  vmap,
  j_start_col = NULL,
  j_end_col = NULL,
  j_chrom_col = NULL
)
```

## Arguments

- gr:

  A `GRanges`.

- vmap:

  A virtualization map.

- j_start_col, j_end_col, j_chrom_col:

  Optional mcols column names for a secondary position pair (with
  optional companion seqnames). When `j_chrom_col` is `NULL`, the j
  seqname is taken to be the same as the primary GRanges seqname.

## Value

A virtualized `GRanges`, possibly shorter than `gr`.
