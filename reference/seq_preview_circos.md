# Preview a SeqPlotR plot as a circos layout

Renders a schematic circular layout — no data required. Each track is
shown as a colored arc sector. Rows (defined via `direction = "under"` /
`%__%`) map to concentric rings from outer to inner; tracks within a row
(defined via `direction = "right"` / `%|%`) map to arc sectors with
angular width proportional to `track_width`. Ring radial thickness is
proportional to the maximum `track_height` in that row.

Only positional layouts are supported — a plot built with a patchwork
layout string has no natural circos interpretation and is rejected.

## Usage

``` r
seq_preview_circos(
  plot_obj = NULL,
  labels = TRUE,
  colors = NULL,
  start_angle = 90,
  end_angle = -270,
  gap_degrees = 2,
  ring_gap = 0.02,
  outer_radius = 0.45,
  inner_radius = 0.08,
  cx = 0.5,
  cy = 0.5
)
```

## Arguments

- plot_obj:

  A `seq_plot` object built with the operator chain. Required.

- labels:

  Logical. Draw `track_id` labels. Default `TRUE`.

- colors:

  Optional named character vector mapping track IDs to fill colors.
  Unspecified IDs receive automatic colors.

- start_angle:

  Numeric clock-face degrees where the first sector starts. Default `90`
  (top). Sectors sweep clockwise (decreasing angle).

- end_angle:

  Numeric clock-face degrees where the last sector ends. Default `-270`
  (full circle = `start_angle - 360`).

- gap_degrees:

  Numeric blank gap in degrees between sectors of the same ring. Default
  `2`.

- ring_gap:

  Numeric npc blank gap between concentric rings. Default `0.02`.

- outer_radius:

  Numeric npc radius of the outermost ring's outer edge. Default `0.45`.

- inner_radius:

  Numeric npc radius of the innermost ring's inner edge. Default `0.08`.

- cx, cy:

  Numeric npc coordinates of the circle centre. Default `(0.5, 0.5)`.

## Value

Named list of polar bounding boxes (invisibly). Each entry is
`list(theta0, theta1, r0, r1)` in clock-face degrees and npc radii,
keyed by track ID.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- seq_plot() %|%
  seq_track(track_id = "Chr1", track_width = 3) %|%
  seq_track(track_id = "Chr2", track_width = 2) %|%
  seq_track(track_id = "Chr3", track_width = 1) %__%
  seq_track(track_id = "Signal") %__%
  seq_track(track_id = "CopyNum")
seq_preview_circos(plot_obj = p)
} # }
```
