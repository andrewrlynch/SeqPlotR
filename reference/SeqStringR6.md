# SeqString R6 class

SeqString R6 class

SeqString R6 class

## Details

Internal R6 generator backing
[`seq_string()`](http://andrewlynch.io/SeqPlotR/reference/seq_string.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Draws a smooth cubic-Bezier "string" between two anchors, typically used
to connect breakpoints between two tracks. C vs. S curve shape is
inferred from the resolved `strand0`/`strand1` fields when
`aes(type = "auto")` (the default); `aes(type = "c")` or
`aes(type = "s")` forces the shape.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqString`

## Methods

### Public methods

- [`SeqStringR6$new()`](#method-SeqString-new)

- [`SeqStringR6$prep()`](#method-SeqString-prep)

- [`SeqStringR6$draw()`](#method-SeqString-draw)

- [`SeqStringR6$clone()`](#method-SeqString-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)
- [`SeqPlotR::SeqLink$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqStringR6.

#### Usage

    SeqStringR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame` carrying both anchors.

- `mapping`:

  Optional `SeqMap`. Required: `x0`, `x1`. Optional: `chrom0`, `chrom1`,
  `strand0`, `strand1`, `y0`, `y1`.

- `t0, t1`:

  Track identifiers. Locked to the parent track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
  via `%+%`.

- `aesthetics`:

  Optional `SeqAes`. Recognised: `color`, `linewidth`, `alpha`, `type`
  (`"auto"`, `"c"`, `"s"`), `bulge`, `orientation`.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Resolve both anchors, find overlaps with each referenced track's
windows, and store per-link canvas coordinates + resolved curve type in
`self$coordCanvas`.

#### Usage

    SeqStringR6$prep(
      layout_all_tracks,
      track_windows_list,
      plot_track_index = NULL
    )

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback track reference when `t0`/`t1` are `NULL`.

------------------------------------------------------------------------

### Method `draw()`

Draw each prepared string with
[`drawSeqString()`](http://andrewlynch.io/SeqPlotR/reference/drawSeqString.md).

#### Usage

    SeqStringR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqStringR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
