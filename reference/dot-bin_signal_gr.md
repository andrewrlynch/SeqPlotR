# Bin a GRanges signal track into fixed-width genomic bins

Assigns each range in `gr` (weighted by overlap fraction) to bins of
width `bin_size` covering `region`, then aggregates using `agg_fun`.

## Usage

``` r
.bin_signal_gr(gr, region, bin_size, agg_fun)
```

## Arguments

- gr:

  GRanges with a numeric `score` mcols column.

- region:

  Single-range GRanges defining the window.

- bin_size:

  Integer bin width in bp.

- agg_fun:

  Resolved aggregation function (from
  [`.resolve_bin_fun()`](http://andrewlynch.io/SeqPlotR/reference/dot-resolve_bin_fun.md)).

## Value

data.frame with columns: seqnames, start, end, score.
