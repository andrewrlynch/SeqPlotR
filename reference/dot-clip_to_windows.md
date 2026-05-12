# Clip genomic ranges to an x-axis scale

Returns the subset of paired ranges that intersect the closed interval
`[xscale[1], xscale[2]]`, with each range trimmed to fit inside that
interval. Out-of-bounds ranges are silently dropped — no warning.

## Usage

``` r
.clip_to_windows(x0, x1, xscale)
```

## Arguments

- x0:

  Numeric vector of range starts.

- x1:

  Numeric vector of range ends (same length as `x0`).

- xscale:

  Numeric vector of length 2: `c(min, max)`.

## Value

A list with elements:

- `x0`:

  clipped starts

- `x1`:

  clipped ends

- `mask`:

  logical vector aligned with the input that marks which rows survived
  clipping (used to subset companion vectors).
