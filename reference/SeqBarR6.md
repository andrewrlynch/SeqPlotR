# SeqBar R6 class

SeqBar R6 class

SeqBar R6 class

## Details

Internal R6 generator backing
[`seq_bar()`](http://andrewlynch.io/SeqPlotR/reference/seq_bar.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqBar`

## Methods

### Public methods

- [`SeqBarR6$new()`](#method-SeqBar-new)

- [`SeqBarR6$prep()`](#method-SeqBar-prep)

- [`SeqBarR6$draw()`](#method-SeqBar-draw)

- [`SeqBarR6$clone()`](#method-SeqBar-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqBarR6.

#### Usage

    SeqBarR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Recognised: `x`, `y`, `group`, `fill`.

- `aesthetics`:

  Optional `SeqAes`. Supports `fill`, `color`, `linewidth`, `alpha`.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc. When a `group`
mapping is present, bars stack at identical x positions.

#### Usage

    SeqBarR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw bars via
[`grid::grid.rect()`](https://rdrr.io/r/grid/grid.rect.html).

#### Usage

    SeqBarR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqBarR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
