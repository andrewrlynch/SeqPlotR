# SeqIdeogram R6 class

SeqIdeogram R6 class

SeqIdeogram R6 class

## Details

Internal R6 generator backing
[`seq_ideogram()`](http://andrewlynch.io/SeqPlotR/reference/seq_ideogram.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).
Draws chromosome ideograms from a `GRanges` of cytogenetic bands. Each
band becomes a filled rectangle; consecutive `acen` bands within a
window render as two red triangles meeting at the centromere.

## Super class

`SeqPlotR::SeqElement` -\> `SeqIdeogram`

## Public fields

- `centroPolys`:

  Per-window list of centromere triangle polygons populated by `prep()`.
  Each entry is `NULL` or `list(list(x, y), list(x, y))`.

- `highlightBoxes`:

  Per-window list of highlight rectangles `list(x0, x1, y0, y1)`.
  Populated only when `scope = "full"`.

- `scope`:

  Character. `"window"` (default) or `"full"`.

- `style`:

  Character. `"block"` (default) or `"rounded"`.

- `highlight_range`:

  Optional `GRanges`. When set (and `scope = "full"`), the highlight
  rectangle covers this range instead of the parent track's windows.

## Methods

### Public methods

- [`SeqIdeogramR6$new()`](#method-SeqIdeogram-new)

- [`SeqIdeogramR6$prep()`](#method-SeqIdeogram-prep)

- [`SeqIdeogramR6$draw()`](#method-SeqIdeogram-draw)

- [`SeqIdeogramR6$clone()`](#method-SeqIdeogram-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqIdeogramR6.

#### Usage

    SeqIdeogramR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      scope = "window",
      style = "block",
      highlight_range = NULL,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` of cytogenetic bands with a `gieStain` mcol (or a
  `stain` mapping in `mapping`). Falls back to the parent track's data.

- `mapping`:

  Optional `SeqMap`. Recognised: `stain` (defaults to the `gieStain`
  mcol when absent).

- `aesthetics`:

  Optional `SeqAes`. Recognised: `color` (band border, default
  `"black"`), `linewidth` (band border width, default `0.1`), `outline`
  (nested aes for the chromosome's perimeter outline; sub-keys: `col`,
  `lwd`, `visible`), `highlight` (nested aes for the `scope = "full"`
  highlight box; sub-keys: `fill`, `col`, `lwd`, `alpha`),
  `telomere.radius` (corner radius as a fraction of band height for
  `style = "rounded"`; default `0.4`).

- `scope`:

  Character. `"window"` (default — only bands overlapping the track
  windows are drawn) or `"full"` (the whole chromosome is drawn rescaled
  to the panel; the current window is overlaid as a highlight
  rectangle).

- `style`:

  Character. `"block"` (default — rectangular bands) or `"rounded"`
  (rounded caps on the leftmost and rightmost bands).

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Resolve mapping, find band overlaps with each window, and populate
`coordCanvas` (non-centromeric bands) and `centroPolys` (per-window
centromere triangles). When `scope = "full"`, the whole chromosome is
mapped into the panel and `highlightBoxes` is populated with one entry
per window.

#### Usage

    SeqIdeogramR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw non-centromeric bands with
[`grid::grid.rect()`](https://rdrr.io/r/grid/grid.rect.html) and
centromere triangles with
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html). With
`style = "rounded"`, the leftmost / rightmost bands draw with rounded
telomere caps. With `scope = "full"`, the current window region is
overlaid as a highlight rectangle.

#### Usage

    SeqIdeogramR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqIdeogramR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
