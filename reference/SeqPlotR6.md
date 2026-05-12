# SeqPlot R6 class

SeqPlot R6 class

SeqPlot R6 class

## Details

Internal R6 class backing
[`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md).
Public users should go through the snake_case constructor.

## Public fields

- `rows`:

  List of lists of SeqTrackR6 — positional layout. `NULL` in patchwork
  mode.

- `layout_str`:

  Raw layout string — patchwork mode (authoritative).

- `aesthetics`:

  SeqAes for plot-wide aesthetics.

- `tracks`:

  Flat list of SeqTrackR6 — populated in patchwork mode.

- `plot_links`:

  Plot-level SeqLink objects — deferred, drawn last.

- `plot_annotations`:

  Plot-level SeqAnnotation objects — deferred.

- `layout`:

  Layout metadata produced by `$layoutGrid()` — list with `panelBounds`
  and `trackBounds`.

- `flat_theme`:

  The merged flat theme map (plot-level default theme

  - user overrides). Populated at `layoutGrid()` time.

- `show_legend`:

  Logical. Global legend switch. When `FALSE`, no legend output is
  produced for any track. Default `TRUE`.

- `windows`:

  Optional plot-wide `GRanges`. Inherited by any track whose `windows`
  is `NULL`. Genomic scales without explicit `windows` also inherit from
  their enclosing track.

- `plot_margin`:

  Scalar (or length-4 `c(b, l, t, r)`) used as the outermost canvas
  margin. Defaults to 0.02.

- `highlight_margins`:

  Logical. When `TRUE`, draw translucent overlays on every margin band —
  a development aid.

## Methods

### Public methods

- [`SeqPlotR6$new()`](#method-SeqPlot-new)

- [`SeqPlotR6$addTrack()`](#method-SeqPlot-addTrack)

- [`SeqPlotR6$allTracks()`](#method-SeqPlot-allTracks)

- [`SeqPlotR6$trackIds()`](#method-SeqPlot-trackIds)

- [`SeqPlotR6$layoutGrid()`](#method-SeqPlot-layoutGrid)

- [`SeqPlotR6$drawGrid()`](#method-SeqPlot-drawGrid)

- [`SeqPlotR6$drawGridlines()`](#method-SeqPlot-drawGridlines)

- [`SeqPlotR6$drawAxes()`](#method-SeqPlot-drawAxes)

- [`SeqPlotR6$drawElements()`](#method-SeqPlot-drawElements)

- [`SeqPlotR6$drawLegends()`](#method-SeqPlot-drawLegends)

- [`SeqPlotR6$plot()`](#method-SeqPlot-plot)

- [`SeqPlotR6$drawMarginHighlights()`](#method-SeqPlot-drawMarginHighlights)

- [`SeqPlotR6$clone()`](#method-SeqPlot-clone)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqPlotR6.

#### Usage

    SeqPlotR6$new(
      layout = NULL,
      aesthetics = aes(),
      tracks = list(),
      show_legend = TRUE,
      legend = NULL,
      windows = NULL,
      plot_margin = 0.02,
      highlight_margins = FALSE,
      ...
    )

#### Arguments

- `layout`:

  Either NULL (positional layout, default) or a multiline character
  string defining a patchwork layout.

- `aesthetics`:

  A SeqAes object.

- `tracks`:

  Optional list of `SeqTrackR6` objects to pre-populate the plot.
  Elements can also be added with `%+%`. When supplied with no `layout`
  string, tracks are placed in one positional row.

- `show_legend`:

  Logical. Global legend switch. Default `TRUE`.

- `legend`:

  Convenience argument. `legend = FALSE` is sugar for
  `show_legend = FALSE`. All other values are ignored.

- `windows`:

  Optional plot-wide `GRanges`. Tracks that omit `windows` inherit from
  this value.

- `plot_margin`:

  Numeric. Outermost canvas margin (default 0.02).

- `highlight_margins`:

  Logical. Draw translucent margin overlays for development (default
  FALSE).

- `...`:

  Additional arguments (currently ignored).

------------------------------------------------------------------------

### Method `addTrack()`

Add a track. Uses `direction` in positional mode; appends to the flat
tracks list in patchwork mode.

#### Usage

    SeqPlotR6$addTrack(track)

#### Arguments

- `track`:

  A SeqTrackR6 instance.

------------------------------------------------------------------------

### Method `allTracks()`

Return all tracks as a flat list regardless of mode.

#### Usage

    SeqPlotR6$allTracks()

------------------------------------------------------------------------

### Method `trackIds()`

Return a character vector of every track's `track_id` in the order
tracks were added. Tracks without a `track_id` contribute
`NA_character_`.

#### Usage

    SeqPlotR6$trackIds()

------------------------------------------------------------------------

### Method `layoutGrid()`

Compute the layout grid for all tracks and windows. Errors if any track
has `windows = NULL` (every track must have windows defined in SeqPlotR
— there are no global plot-level windows). Builds `scale_x` from
`seq_scale_genomic(windows)` when missing, and auto-infers `scale_y`
from element data when missing. Then dispatches to either
[`.build_positional_layout()`](http://andrewlynch.io/SeqPlotR/reference/dot-build_positional_layout.md)
or
[`.build_patchwork_layout()`](http://andrewlynch.io/SeqPlotR/reference/dot-build_patchwork_layout.md)
and opens a `grid` viewport.

#### Usage

    SeqPlotR6$layoutGrid()

#### Returns

The plot, invisibly.

------------------------------------------------------------------------

### Method `drawGrid()`

Draw track backgrounds / borders and per-window panel chrome. Delegates
to
[`.draw_track_chrome()`](http://andrewlynch.io/SeqPlotR/reference/dot-draw_track_chrome.md)
for each track using the track's resolved theme.

#### Usage

    SeqPlotR6$drawGrid()

------------------------------------------------------------------------

### Method `drawGridlines()`

Draw x and y gridlines at axis break positions for all track windows.
Gridlines sit after window backgrounds and before elements. Enable per
axis via `axis.x.gridline = TRUE` (or `= aes(color, lwd, lty, alpha)`)
and `axis.y.gridline = TRUE` in
[`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md) or
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
aesthetics. Styling inherits from `axis.gridline.*` defaults in the
theme hierarchy.

#### Usage

    SeqPlotR6$drawGridlines()

#### Returns

Renders gridlines to the graphics device; returns invisibly.

------------------------------------------------------------------------

### Method `drawAxes()`

Draw x and y axes for every track. Delegates to
[`.draw_track_axes()`](http://andrewlynch.io/SeqPlotR/reference/dot-draw_track_axes.md),
which reads the track's `resolved_theme` and renders up to four axes
(x1, x2, y1, y2) with hierarchical aesthetic control, break computation,
and cap modes.

#### Usage

    SeqPlotR6$drawAxes()

------------------------------------------------------------------------

### Method `drawElements()`

Draw all elements in drawing order: within-track non-links, then
within-track links, then plot-level deferred links, then plot-level
annotations.

#### Usage

    SeqPlotR6$drawElements()

------------------------------------------------------------------------

### Method `drawLegends()`

Draw legends for all tracks. Dispatches on `position` in each
`SeqLegendSpec` found on element `legend` fields. Bare `LegendKey` or
list-of-`LegendKey` on an element are automatically wrapped in a default
`"inside"` spec. Phase 1 handles `"inside"` and `"track_margin"` per
element; Phase 2 aggregates all `"canvas_margin"` specs and draws once
per side.

Call after `drawElements()`.

#### Usage

    SeqPlotR6$drawLegends()

#### Returns

Renders legends to the graphics device; returns invisibly.

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Run the full plot pipeline: `layoutGrid()` -\> `drawGrid()` -\>
`drawAxes()` -\> `drawElements()` -\> `drawLegends()`. Resets
base-graphics `par(mar, oma, mai, omi)` to zero so that any residual
margins inherited from the active device (knitr's PNG device, RStudio's
plot pane, etc.) do not leave whitespace around the grid viewport.

#### Usage

    SeqPlotR6$plot()

------------------------------------------------------------------------

### Method `drawMarginHighlights()`

Overlay translucent coloured rects on every margin band for development.
Triggered by `seq_plot(highlight_margins = TRUE)`. Colours: plot margin
red, outer track dark blue, inner track light blue, outer window dark
green, inner window light green, plotting area pink. All drawn at
`alpha = 0.5`.

#### Usage

    SeqPlotR6$drawMarginHighlights()

#### Returns

Invisibly `NULL`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPlotR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
