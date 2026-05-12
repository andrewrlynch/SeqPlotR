# Resolve a SeqAes + SeqScale into concrete per-observation vectors

Resolve a SeqAes + SeqScale into concrete per-observation vectors

## Usage

``` r
.resolve_aes(data_mcols, aes_obj, scale_obj, n, default_color = "#1C1B1A")
```

## Arguments

- data_mcols:

  data.frame of GRanges metadata columns.

- aes_obj:

  A `SeqAes` object (from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md)), or
  `NULL`.

- scale_obj:

  A `SeqScale` object (from `seq_scale_*()`), or `NULL`.

- n:

  Number of observations.

- default_color:

  Fallback color when no mapping is specified.

## Value

Named list with resolved vectors for `color`, `fill`, `alpha`, `size`,
`shape` (only those present in `aes_obj`).
