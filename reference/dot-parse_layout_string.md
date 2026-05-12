# Parse a patchwork layout string

Converts a multiline character string into a structured layout
description. Each non-`#` letter denotes a region; `#` cells are
rendered as blank.

## Usage

``` r
.parse_layout_string(s)
```

## Arguments

- s:

  A single multiline character string. Rows are split on `"\n"`,
  trimmed, and dropped if empty. All rows must have the same number of
  characters.

## Value

A list with elements:

- `nrow`:

  number of rows

- `ncol`:

  number of columns

- `regions`:

  named list, one entry per unique non-`#` letter, each a list
  `list(r0, r1, c0, c1)` (1-indexed, inclusive)

- `blank_cells`:

  list of `list(row, col)` for every `#` cell
