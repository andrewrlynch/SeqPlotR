# SeqPoly R6 class

SeqPoly R6 class

SeqPoly R6 class

## Details

Internal R6 generator backing
[`seq_poly()`](http://andrewlynch.io/SeqPlotR/reference/seq_poly.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqPoly`

## Methods

### Public methods

- [`SeqPolyR6$new()`](#method-SeqPoly-new)

- [`SeqPolyR6$prep()`](#method-SeqPoly-prep)

- [`SeqPolyR6$draw()`](#method-SeqPoly-draw)

- [`SeqPolyR6$clone()`](#method-SeqPoly-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqPolyR6.

#### Usage

    SeqPolyR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

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

    SeqPolyR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw polygons via
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html) with
`id`.

#### Usage

    SeqPolyR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPolyR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
