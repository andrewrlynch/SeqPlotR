# SeqSegment R6 class

SeqSegment R6 class

SeqSegment R6 class

## Details

Internal R6 generator backing
[`seq_segment()`](http://andrewlynch.io/SeqPlotR/reference/seq_segment.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqSegment`

## Methods

### Public methods

- [`SeqSegmentR6$new()`](#method-SeqSegment-new)

- [`SeqSegmentR6$prep()`](#method-SeqSegment-prep)

- [`SeqSegmentR6$draw()`](#method-SeqSegment-draw)

- [`SeqSegmentR6$clone()`](#method-SeqSegment-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqSegmentR6.

#### Usage

    SeqSegmentR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required fields: `x`, `x_end`, `y`, `y_end`.

- `aesthetics`:

  Optional `SeqAes`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, apply auto-scale for non-concrete color, clip to
windows, transform to npc.

#### Usage

    SeqSegmentR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw segments via
[`grid::grid.segments()`](https://rdrr.io/r/grid/grid.segments.html).

#### Usage

    SeqSegmentR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqSegmentR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
