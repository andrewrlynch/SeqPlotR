# SeqSynteny R6 class

SeqSynteny R6 class

SeqSynteny R6 class

## Details

Internal R6 generator backing
[`seq_synteny()`](http://andrewlynch.io/SeqPlotR/reference/seq_synteny.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Connects a genomic region in track `t0` to a homologous region in track
`t1` with a filled quadrilateral. Each region is defined by a pair of
edges (`x0`/`x0_end` and `x1`/`x1_end`) plus an optional data-scale y
for each side.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqSynteny`

## Methods

### Public methods

- [`SeqSyntenyR6$new()`](#method-SeqSynteny-new)

- [`SeqSyntenyR6$prep()`](#method-SeqSynteny-prep)

- [`SeqSyntenyR6$draw()`](#method-SeqSynteny-draw)

- [`SeqSyntenyR6$clone()`](#method-SeqSynteny-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)
- [`SeqPlotR::SeqLink$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqSyntenyR6.

#### Usage

    SeqSyntenyR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame` carrying both region pairs.

- `mapping`:

  Optional `SeqMap`. Required: `x0`, `x1` (and, for `data.frame`,
  `chrom0`, `chrom1`). Optional: `x0_end`, `x1_end`, `y0`, `y1`,
  `color`, `fill`.

- `t0, t1`:

  Track identifiers. Locked to the parent track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
  via `%+%`.

- `aesthetics`:

  Optional `SeqAes`. Recognised: `fill`, `color`, `alpha`, `linewidth`.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Resolve both region endpoints, find which windows each region's left
anchor falls into, and store a 4-corner polygon per region in
`self$coordCanvas`.

#### Usage

    SeqSyntenyR6$prep(
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

  Fallback track reference.

------------------------------------------------------------------------

### Method `draw()`

Draw each prepared polygon with
[`grid::grid.polygon()`](https://rdrr.io/r/grid/grid.polygon.html).

#### Usage

    SeqSyntenyR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqSyntenyR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
