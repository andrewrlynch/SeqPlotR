# SeqZoom R6 class

SeqZoom R6 class

SeqZoom R6 class

## Details

Internal R6 generator backing
[`seq_zoom()`](http://andrewlynch.io/SeqPlotR/reference/seq_zoom.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Draws a filled quadrilateral projecting a genomic region in track `t0`
onto the same (or a different) region in track `t1`. Typically used to
connect an overview / ideogram track to a zoomed detail track.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqZoom`

## Methods

### Public methods

- [`SeqZoomR6$new()`](#method-SeqZoom-new)

- [`SeqZoomR6$resolve()`](#method-SeqZoom-resolve)

- [`SeqZoomR6$prep()`](#method-SeqZoom-prep)

- [`SeqZoomR6$draw()`](#method-SeqZoom-draw)

- [`SeqZoomR6$clone()`](#method-SeqZoom-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqZoomR6.

#### Usage

    SeqZoomR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame` carrying the region(s).

- `mapping`:

  Optional `SeqMap`. Required: `x0`, `x0_end`. Optional: `x1`, `x1_end`
  (default to `x0`, `x0_end` when absent), `chrom0`, `chrom1`.

- `t0, t1`:

  Track identifiers. Must be specified explicitly for plot-level
  placement.

- `aesthetics`:

  Optional `SeqAes`. Recognised: `fill`, `color`, `alpha`, `linewidth`,
  `stemOffset` (npc, default `0.01`).

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `resolve()`

Resolve the mapping and populate region fields. Unlike other `SeqLink`
subclasses, `seq_zoom` only requires `x0` / `x0_end` — the `x1` /
`x1_end` defaults mirror the `x0` region when absent.

#### Usage

    SeqZoomR6$resolve(track_data = NULL, track_mapping = NULL)

#### Arguments

- `track_data`:

  Optional `GRanges`/`data.frame` from the parent track.

- `track_mapping`:

  Optional `SeqMap` from the parent track.

#### Returns

The link, invisibly.

------------------------------------------------------------------------

### Method `prep()`

Project the region(s) onto both tracks and build a `data.frame` of
4-corner polygon coordinates in `self$coordCanvas`.

#### Usage

    SeqZoomR6$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Unused — `seq_zoom` requires explicit `t0`/`t1`.

------------------------------------------------------------------------

### Method `draw()`

Draw each prepared polygon with
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html).

#### Usage

    SeqZoomR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqZoomR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
