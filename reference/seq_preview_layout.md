# Preview the grid layout of a SeqPlotR plot

Renders a schematic toy representation of the intended plot layout — no
data required. Each track region is drawn as a filled, labeled
rectangle. Useful for verifying a layout string or operator chain before
adding data.

## Usage

``` r
seq_preview_layout(
  plot_obj = NULL,
  layout = NULL,
  labels = TRUE,
  colors = NULL,
  margins = TRUE,
  margin_size = 0.05
)
```

## Arguments

- plot_obj:

  A `seq_plot` object built with the operator chain, or `NULL`. If
  provided, the layout is extracted from the object.

- layout:

  A multiline layout string (e.g. `"##AA\n##AA\nBBBC\nBBBD"`), or
  `NULL`. Used when `plot_obj` is not provided.

- labels:

  Logical. If `TRUE` (default), draw the track ID label centered in each
  region.

- colors:

  Optional named character vector mapping track IDs to fill colors (e.g.
  `c(A = "#FF0000", B = "#00FF00")`). Unspecified IDs receive automatic
  colors.

- margins:

  Logical. If `TRUE` (default), draw a thin border rectangle indicating
  the outer canvas margin.

- margin_size:

  Numeric in the range 0 to 0.5. Fractional canvas margin on each side.
  Default `0.05`.

## Value

A named list of npc bounding boxes (invisibly): each element is
`list(x0, x1, y0, y1)` keyed by track ID.

## Examples

``` r
if (FALSE) { # \dontrun{
# Preview a patchwork layout string directly
layout <- "
##AA
##AA
BBBC
BBBD
"
seq_preview_layout(layout = layout)

# Preview a layout built with operators
p <- seq_plot() |>
  (\(x) x %|% seq_track(track_id = "Signal",  track_width = 3))() |>
  (\(x) x %|% seq_track(track_id = "Ideogram", track_width = 1))() |>
  (\(x) x %__% seq_track(track_id = "Genes"))()
seq_preview_layout(plot_obj = p)
} # }
```
