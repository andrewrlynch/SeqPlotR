# SeqArch R6 class

R6 generator backing
[`seq_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_arch.md),
exported for subclassing (for example
[`seq_cross_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_cross_arch.md)).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Same anchor pattern as
[`seq_arc()`](http://andrewlynch.io/SeqPlotR/reference/seq_arc.md) but
draws full vertical stems from the baseline to the arch endpoints and
renders partial-window stubs when only one anchor is visible.

## Super classes

`SeqElement` -\> `SeqLink` -\> `SeqArch`

## Public fields

- `stubs`:

  List of partial-window stub descriptions populated by `prep()`. Each
  entry: `list(x, y0, y1, dir, partner, color, width)`.

## Methods

### Public methods

- [`SeqArch$new()`](#method-SeqArch-initialize)

- [`SeqArch$prep()`](#method-SeqArch-prep)

- [`SeqArch$draw()`](#method-SeqArch-draw)

- [`SeqArch$clone()`](#method-SeqArch-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqLink$.resolve_track_ref()`
- `SeqLink$resolve()`

------------------------------------------------------------------------

### `SeqArch$new()`

Construct a SeqArchR6.

#### Usage

    SeqArch$new(
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

### `SeqArch$prep()`

Resolve both anchors against their referenced tracks, assign each anchor
to a window via `findOverlaps()`, and build `coordCanvas` for
fully-visible links plus `stubs` for half-visible links.

#### Usage

    SeqArch$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback identifier for `t0` / `t1` when both are `NULL` (within-track
  placement).

------------------------------------------------------------------------

### `SeqArch$draw()`

Draw arches with stems, then partial-window stubs (when
`aes(plotStubs = TRUE)`).

#### Usage

    SeqArch$draw()

------------------------------------------------------------------------

### `SeqArch$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqArch$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
