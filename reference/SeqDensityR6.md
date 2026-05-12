# SeqDensity R6 class

SeqDensity R6 class

SeqDensity R6 class

## Details

Internal R6 generator backing
[`seq_density()`](http://andrewlynch.io/SeqPlotR/reference/seq_density.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqDensity`

## Methods

### Public methods

- [`SeqDensityR6$new()`](#method-SeqDensity-new)

- [`SeqDensityR6$prep()`](#method-SeqDensity-prep)

- [`SeqDensityR6$draw()`](#method-SeqDensity-draw)

- [`SeqDensityR6$clone()`](#method-SeqDensity-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqDensityR6.

#### Usage

    SeqDensityR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `y`.

- `aesthetics`:

  Optional `SeqAes`. `bw` controls the kernel bandwidth (passed to
  [`stats::density()`](https://rdrr.io/r/stats/density.html)). Default
  fill `"grey60"`, alpha 0.8.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Compute [`stats::density()`](https://rdrr.io/r/stats/density.html) of
the resolved y, then build a closed polygon. Density evaluation points
are mapped to canvas x via the panel's `yscale`; heights are normalised
to `[0, 1]` of panel height.

#### Usage

    SeqDensityR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw the density polygon via
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html).

#### Usage

    SeqDensityR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqDensityR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
