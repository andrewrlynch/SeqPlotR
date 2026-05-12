# SeqLine R6 class

SeqLine R6 class

SeqLine R6 class

## Details

Internal R6 generator backing
[`seq_line()`](http://andrewlynch.io/SeqPlotR/reference/seq_line.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqLine`

## Methods

### Public methods

- [`SeqLineR6$new()`](#method-SeqLine-new)

- [`SeqLineR6$prep()`](#method-SeqLine-prep)

- [`SeqLineR6$draw()`](#method-SeqLine-draw)

- [`SeqLineR6$clone()`](#method-SeqLine-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqLineR6.

#### Usage

    SeqLineR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`.

- `aesthetics`:

  Optional `SeqAes`. `type = "step"` (or `"s"`) enables step-line
  rendering.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc.

#### Usage

    SeqLineR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw the polyline with
[`grid::grid.lines()`](https://rdrr.io/r/grid/grid.lines.html).

#### Usage

    SeqLineR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqLineR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
