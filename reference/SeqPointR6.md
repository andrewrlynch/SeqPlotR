# SeqPoint R6 class

SeqPoint R6 class

SeqPoint R6 class

## Details

Internal R6 generator backing
[`seq_point()`](http://andrewlynch.io/SeqPlotR/reference/seq_point.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqPoint`

## Methods

### Public methods

- [`SeqPointR6$new()`](#method-SeqPoint-new)

- [`SeqPointR6$prep()`](#method-SeqPoint-prep)

- [`SeqPointR6$draw()`](#method-SeqPoint-draw)

- [`SeqPointR6$clone()`](#method-SeqPoint-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqPointR6.

#### Usage

    SeqPointR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`.

- `aesthetics`:

  Optional `SeqAes`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, apply auto-scales for non-concrete color/fill/shape,
clip to windows, transform genomic coordinates to canvas npc.

#### Usage

    SeqPointR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw points via
[`grid::grid.points()`](https://rdrr.io/r/grid/grid.points.html) for
each window.

#### Usage

    SeqPointR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPointR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
