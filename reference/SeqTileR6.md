# SeqTile R6 class

SeqTile R6 class

SeqTile R6 class

## Details

Internal R6 generator backing
[`seq_tile()`](http://andrewlynch.io/SeqPlotR/reference/seq_tile.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqTile`

## Public fields

- `data2`:

  Optional `GRanges` giving the y-axis genomic coordinates for rotated
  mode (one range per row of the primary `data`).

## Methods

### Public methods

- [`SeqTileR6$new()`](#method-SeqTile-new)

- [`SeqTileR6$prep()`](#method-SeqTile-prep)

- [`SeqTileR6$draw()`](#method-SeqTile-draw)

- [`SeqTileR6$clone()`](#method-SeqTile-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqTileR6.

#### Usage

    SeqTileR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      data2 = NULL,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Recognised: `x`, `y`, `fill`.

- `aesthetics`:

  Optional `SeqAes`. `rotate` toggles diamond mode; `fill` sets a
  constant fill color when no `fill` mapping is given.

- `data2`:

  Optional `GRanges` for rotated mode (y genomic ranges).

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, clip to windows, transform to npc.

#### Usage

    SeqTileR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw tiles. Rectangles via
[`grid::grid.rect()`](https://rdrr.io/r/grid/grid.rect.html) for the
unrotated mode; diamond polygons via
`grid::grid.polygon(id.lengths = ...)` for the rotated mode.

#### Usage

    SeqTileR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqTileR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
