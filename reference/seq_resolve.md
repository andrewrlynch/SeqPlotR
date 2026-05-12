# Compose wrapper-built plots into a parent plot

Each wrapper (e.g.
[`seq_copynumber()`](http://andrewlynch.io/SeqPlotR/reference/seq_copynumber.md),
[`seq_hic()`](http://andrewlynch.io/SeqPlotR/reference/seq_hic.md),
[`seq_chip()`](http://andrewlynch.io/SeqPlotR/reference/seq_chip.md))
returns a `seq_plot`. `seq_resolve()` unpacks those child plots into a
parent `seq_plot`, letting you combine heterogeneous track stacks in a
single figure.

## Usage

``` r
seq_resolve(parent, ..., direction = "under")
```

## Arguments

- parent:

  A `seq_plot` object to add tracks into.

- ...:

  One or more `seq_plot` objects produced by wrapper functions.

- direction:

  `"under"` (default) or `"right"`. See Details.

## Value

The modified `parent`, invisibly (R6 reference semantics).

## Details

Tracks from each child are appended to `parent` respecting their
original row grouping. The `direction` argument controls where the
**first** row of each child lands relative to the current bottom of
`parent`:

- `"under"` — start a new row below whatever is already in the parent
  (the common case).

- `"right"` — append to the current row (put the child beside what is
  already there).

Subsequent rows within a child are always placed `"under"` (i.e. the
child's internal row structure is preserved). Multiple children supplied
in one call are stacked in order, each with `direction` applied between
children.

Plot-level links and annotations carried by the children are transferred
to `parent` in document order.

Duplicate `track_id`s across children cause an error — pass unique
`track_id` values to each wrapper call to resolve the conflict.

## Examples

``` r
if (FALSE) { # \dontrun{
win <- GRanges("chr1", IRanges(1, 1e6))
cn  <- seq_copynumber(cn_gr, windows = win, track_id = "CN")
hic <- seq_hic(hic_gr, windows = win, style = "triangle",
               track_id = "HiC")
seq_resolve(seq_plot(), cn, hic)$plot()
} # }
```
