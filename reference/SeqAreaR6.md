# SeqArea R6 class

SeqArea R6 class

SeqArea R6 class

## Details

Internal R6 generator backing
[`seq_area()`](http://andrewlynch.io/SeqPlotR/reference/seq_area.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqArea`

## Methods

### Public methods

- [`SeqAreaR6$new()`](#method-SeqArea-new)

- [`SeqAreaR6$prep()`](#method-SeqArea-prep)

- [`SeqAreaR6$draw()`](#method-SeqArea-draw)

- [`SeqAreaR6$clone()`](#method-SeqArea-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqAreaR6.

#### Usage

    SeqAreaR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `x`, `y`.

- `aesthetics`:

  Optional `SeqAes`. `baseline` sets the closing y-value (default 0).

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc, then build a closed
polygon by appending a reversed baseline path.

#### Usage

    SeqAreaR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw filled polygons via
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html).

#### Usage

    SeqAreaR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqAreaR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
