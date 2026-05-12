# SeqText R6 class

SeqText R6 class

SeqText R6 class

## Details

Internal R6 generator backing
[`seq_text()`](http://andrewlynch.io/SeqPlotR/reference/seq_text.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqText`

## Methods

### Public methods

- [`SeqTextR6$new()`](#method-SeqText-new)

- [`SeqTextR6$prep()`](#method-SeqText-prep)

- [`SeqTextR6$draw()`](#method-SeqText-draw)

- [`SeqTextR6$clone()`](#method-SeqText-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqTextR6.

#### Usage

    SeqTextR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `x`, `y`, `label`.

- `aesthetics`:

  Optional `SeqAes`. Supports `size`, `color`, `angle`, `hjust`,
  `vjust`, and a constant `label`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip, transform to npc.

#### Usage

    SeqTextR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw text via
[`grid::grid.text()`](https://rdrr.io/r/grid/grid.text.html).

#### Usage

    SeqTextR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqTextR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
