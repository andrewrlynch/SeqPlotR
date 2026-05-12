# Infer a shared x-axis scale factor for a track

Uses the narrowest window to determine the appropriate unit (Mb/kb/bp)
so that all windows in the track share the same scale.

## Usage

``` r
.infer_track_scale_factor(widths)
```

## Arguments

- widths:

  Integer vector of genomic window widths in bp.

## Value

A single numeric scale factor: `1e-6` (Mb), `1e-3` (kb), or `1` (bp).
