# SeqLollipop R6 class

SeqLollipop R6 class

SeqLollipop R6 class

## Details

Internal R6 generator backing
[`seq_lollipop()`](http://andrewlynch.io/SeqPlotR/reference/seq_lollipop.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqLollipop`

## Methods

### Public methods

- [`SeqLollipopR6$new()`](#method-SeqLollipop-new)

- [`SeqLollipopR6$prep()`](#method-SeqLollipop-prep)

- [`SeqLollipopR6$draw()`](#method-SeqLollipop-draw)

- [`SeqLollipopR6$clone()`](#method-SeqLollipop-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqLollipopR6.

#### Usage

    SeqLollipopR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Recognised: `x`, `y`.

- `aesthetics`:

  Optional `SeqAes`. `baseline` sets the stem's lower y-value (default
  0). Also supports `color`, `linewidth`, `size`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc. Builds per-window
`list(x, y0, y1)` used by both the stem and the head.

#### Usage

    SeqLollipopR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw stems via
[`grid::grid.segments()`](https://rdrr.io/r/grid/grid.segments.html) and
heads via
[`grid::grid.points()`](https://rdrr.io/r/grid/grid.points.html).

#### Usage

    SeqLollipopR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqLollipopR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
