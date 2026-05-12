# SeqArch R6 class

SeqArch R6 class

SeqArch R6 class

## Details

Internal R6 generator backing
[`seq_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_arch.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Same anchor pattern as
[`SeqArcR6`](http://andrewlynch.io/SeqPlotR/reference/SeqArcR6.md) but
draws full vertical stems from the baseline to the arch endpoints and
renders partial-window stubs when only one anchor is visible.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqArch`

## Public fields

- `stubs`:

  List of partial-window stub descriptions populated by `prep()`. Each
  entry: `list(x, y0, y1, dir, partner, color, width)`.

## Methods

### Public methods

- [`SeqArchR6$new()`](#method-SeqArch-new)

- [`SeqArchR6$prep()`](#method-SeqArch-prep)

- [`SeqArchR6$draw()`](#method-SeqArch-draw)

- [`SeqArchR6$clone()`](#method-SeqArch-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)
- [`SeqPlotR::SeqLink$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqArchR6.

#### Usage

    SeqArchR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame` with both anchors.

- `mapping`:

  Optional `SeqMap`. Required: `x0`, `x1`. Optional: `chrom0`, `chrom1`,
  `strand0`, `strand1`, `y0`, `y1`, `height`.

- `t0, t1`:

  Track identifiers. Locked to the parent track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
  via `%+%`.

- `aesthetics`:

  Optional `SeqAes`. Recognised: `color`, `linewidth`, `arcColor`,
  `stemColor`, `arcWidth`, `stemWidth`, `orientation`, `curve`,
  `height`, `plotStubs`, `stubAngle`, `stubLength`.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Resolve both anchors against their referenced tracks, assign each anchor
to a window via `findOverlaps()`, and build `coordCanvas` for
fully-visible links plus `stubs` for half-visible links.

#### Usage

    SeqArchR6$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback identifier for `t0` / `t1` when both are `NULL` (within-track
  placement).

------------------------------------------------------------------------

### Method `draw()`

Draw arches with stems, then partial-window stubs (when
`aes(plotStubs = TRUE)`).

#### Usage

    SeqArchR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqArchR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
