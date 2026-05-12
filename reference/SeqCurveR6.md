# SeqCurve R6 class

SeqCurve R6 class

SeqCurve R6 class

## Details

Internal R6 generator backing
[`seq_curve()`](http://andrewlynch.io/SeqPlotR/reference/seq_curve.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqCurve`

## Methods

### Public methods

- [`SeqCurveR6$new()`](#method-SeqCurve-new)

- [`SeqCurveR6$prep()`](#method-SeqCurve-prep)

- [`SeqCurveR6$draw()`](#method-SeqCurve-draw)

- [`SeqCurveR6$clone()`](#method-SeqCurve-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqCurveR6.

#### Usage

    SeqCurveR6$new(data = NULL, mapping = NULL, aesthetics = aes(), ...)

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Required: `x`, `y`, `x_end`, `y_end`.

- `aesthetics`:

  Optional `SeqAes`. `curvature` (default 0.3) sets the fractional
  y-offset of the bezier control points.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, compute bezier grobs per window.

#### Usage

    SeqCurveR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw bezier curves via
[`grid::bezierGrob()`](https://rdrr.io/r/grid/grid.bezier.html).

#### Usage

    SeqCurveR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqCurveR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
