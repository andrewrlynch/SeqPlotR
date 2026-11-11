# SeqPoint R6 class

R6 generator backing
[`seq_point()`](http://andrewlynch.io/SeqPlotR/reference/seq_point.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).
Exported so that other packages can subclass it.

## Super class

`SeqElement` -\> `SeqPoint`

## Methods

### Public methods

- [`SeqPoint$new()`](#method-SeqPoint-initialize)

- [`SeqPoint$prep()`](#method-SeqPoint-prep)

- [`SeqPoint$draw()`](#method-SeqPoint-draw)

- [`SeqPoint$clone()`](#method-SeqPoint-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$resolve()`

------------------------------------------------------------------------

### `SeqPoint$new()`

Construct a SeqPointR6.

#### Usage

    SeqPoint$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

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

### `SeqPoint$prep()`

Resolve mappings, apply auto-scales for non-concrete color/fill/shape,
clip to windows, transform genomic coordinates to canvas npc.

#### Usage

    SeqPoint$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### `SeqPoint$draw()`

Draw points via
[`grid::grid.points()`](https://rdrr.io/r/grid/grid.points.html) for
each window.

#### Usage

    SeqPoint$draw()

------------------------------------------------------------------------

### `SeqPoint$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPoint$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
