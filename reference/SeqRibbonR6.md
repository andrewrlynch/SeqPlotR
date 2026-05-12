# SeqRibbon R6 class

SeqRibbon R6 class

SeqRibbon R6 class

## Details

Internal R6 generator backing
[`seq_ribbon()`](http://andrewlynch.io/SeqPlotR/reference/seq_ribbon.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqRibbon`

## Methods

### Public methods

- [`SeqRibbonR6$new()`](#method-SeqRibbon-new)

- [`SeqRibbonR6$prep()`](#method-SeqRibbon-prep)

- [`SeqRibbonR6$draw()`](#method-SeqRibbon-draw)

- [`SeqRibbonR6$clone()`](#method-SeqRibbon-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqRibbonR6.

#### Usage

    SeqRibbonR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `x`, `y_min`, `y_max`.

- `aesthetics`:

  Optional `SeqAes`. Defaults to grey fill, alpha 0.8.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc, and build a closed
polygon bounded above by `y_max` and below by `y_min`.

#### Usage

    SeqRibbonR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw the ribbon polygon via
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html).

#### Usage

    SeqRibbonR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqRibbonR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
