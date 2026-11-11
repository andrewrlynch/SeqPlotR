# SeqConnector R6 class

Internal R6 generator backing
[`seq_connector()`](http://andrewlynch.io/SeqPlotR/reference/seq_connector.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Draws a cubic Bezier from one edge of an interval to one edge of
another, on the same or different tracks.

## Super classes

`SeqElement` -\> `SeqLink` -\> `SeqConnector`

## Public fields

- `styles`:

  Named list of `list(col, lty, lwd)` keyed by style.

- `height`:

  Block height as a fraction of the panel inner height.

- `handle`:

  Handle length as a fraction of the track inner width.

- `bulge_height`:

  Same-track bulge as a fraction of the panel inner height.

## Methods

### Public methods

- [`SeqConnector$new()`](#method-SeqConnector-initialize)

- [`SeqConnector$prep()`](#method-SeqConnector-prep)

- [`SeqConnector$draw()`](#method-SeqConnector-draw)

- [`SeqConnector$clone()`](#method-SeqConnector-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqLink$.resolve_track_ref()`
- `SeqLink$resolve()`

------------------------------------------------------------------------

### `SeqConnector$new()`

Construct a SeqConnectorR6.

#### Usage

    SeqConnector$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      styles = NULL,
      height = 0.5,
      handle = 0.04,
      bulge_height = 0.45,
      layer = c("above", "below"),
      legend = NULL,
      show_legend = TRUE,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame`.

- `mapping`:

  Optional `SeqMap`.

- `t0, t1`:

  Track identifiers.

- `aesthetics`:

  Optional `SeqAes`.

- `styles`:

  Optional named list of `list(col, lty, lwd)`.

- `height, handle, bulge_height`:

  Geometry fractions.

- `layer`:

  `"above"` or `"below"`.

- `legend`:

  Optional legend.

- `show_legend`:

  Logical.

- `...`:

  Reserved.

------------------------------------------------------------------------

### `SeqConnector$prep()`

Resolve both anchors, map the chosen edges to canvas npc and store four
Bezier control points per connector in `self$coordCanvas`.

#### Usage

    SeqConnector$prep(
      layout_all_tracks,
      track_windows_list,
      plot_track_index = NULL
    )

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback track reference.

------------------------------------------------------------------------

### `SeqConnector$draw()`

Draw each prepared connector with
[`grid::grid.bezier()`](https://rdrr.io/r/grid/grid.bezier.html).

#### Usage

    SeqConnector$draw()

------------------------------------------------------------------------

### `SeqConnector$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqConnector$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
