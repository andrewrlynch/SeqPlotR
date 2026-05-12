# Build a virtualization map for a multi-region GRanges

Concatenates the windows into a single virtual coordinate system in
insertion order. Each window contributes `width(window)` to the virtual
axis, with optional gaps (`virtual_gap`) inserted between windows so a
separator can be drawn at the boundaries.

## Usage

``` r
.build_virtual_map(windows, virtual_gap = 0)
```

## Arguments

- windows:

  A `GRanges` of one or more windows.

- virtual_gap:

  Numeric (bp) gap inserted between consecutive windows in virtual
  coordinates. Default `0`.

## Value

A list with components:

- `seqnames`:

  character vector of original seqnames, length `length(windows)`.

- `genomic_start`, `genomic_end`:

  numeric vectors of original genomic ranges.

- `virtual_start`, `virtual_end`:

  numeric vectors of the corresponding virtual ranges.

- `virtual_total`:

  the total virtual extent (right edge of the last window).

- `virtual_gap`:

  the gap size used (echoed back).

- `combined_window`:

  a single-range `GRanges` covering the full virtual extent, used in
  place of `windows` for layout.
