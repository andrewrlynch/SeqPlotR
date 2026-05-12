# SeqArc R6 class

SeqArc R6 class

SeqArc R6 class

## Details

Internal R6 generator backing
[`seq_arc()`](http://andrewlynch.io/SeqPlotR/reference/seq_arc.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Draws a single Bezier arch between two genomic loci within one track.
Both anchors live on the same track (`t0 == t1`), and `%+%` locks them
to the parent track when added inside a
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqArc`

## Methods

### Public methods

- [`SeqArcR6$new()`](#method-SeqArc-new)

- [`SeqArcR6$prep()`](#method-SeqArc-prep)

- [`SeqArcR6$draw()`](#method-SeqArc-draw)

- [`SeqArcR6$clone()`](#method-SeqArc-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)
- [`SeqPlotR::SeqLink$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqArcR6.

#### Usage

    SeqArcR6$new(
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
  `y0`, `y1`, `height`.

- `t0, t1`:

  Track identifiers. Both default to the parent track when the arc is
  added via `%+%` inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

- `aesthetics`:

  Optional `SeqAes`. Recognised: `color`, `linewidth`, `orientation`,
  `curve`, `height`.

- `...`:

  Unused — accepted for forward compatibility.

------------------------------------------------------------------------

### Method `prep()`

Resolve the arc's mapping against the referenced track's panels, find
overlaps with the track's windows, and store per-link canvas coordinates
for `draw()`.

#### Usage

    SeqArcR6$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback track reference when `t0` is `NULL`.

------------------------------------------------------------------------

### Method `draw()`

Draw each prepared arch with
[`drawSeqArch()`](http://andrewlynch.io/SeqPlotR/reference/drawSeqArch.md).
No stems are rendered.

#### Usage

    SeqArcR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqArcR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
