# SeqPath R6 class

SeqPath R6 class

SeqPath R6 class

## Details

Internal R6 generator backing
[`seq_path()`](http://andrewlynch.io/SeqPlotR/reference/seq_path.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqPath`

## Methods

### Public methods

- [`SeqPathR6$new()`](#method-SeqPath-new)

- [`SeqPathR6$prep()`](#method-SeqPath-prep)

- [`SeqPathR6$draw()`](#method-SeqPath-draw)

- [`SeqPathR6$clone()`](#method-SeqPath-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqPathR6.

#### Usage

    SeqPathR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `x`, `y`. Optional: `group`.

- `aesthetics`:

  Optional `SeqAes`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip, transform, partition by `group`.

#### Usage

    SeqPathR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw polylines via
[`grid::grid.polyline()`](https://rdrr.io/r/grid/grid.lines.html) with
`id`.

#### Usage

    SeqPathR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPathR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
