# SeqRecon R6 class

SeqRecon R6 class

SeqRecon R6 class

## Details

Internal R6 generator backing
[`seq_recon()`](http://andrewlynch.io/SeqPlotR/reference/seq_recon.md).
Inherits from
[`SeqArchR6`](http://andrewlynch.io/SeqPlotR/reference/SeqArchR6.md).
Classifies each link by the strand pair on its two anchors:

- `+/+` head-to-head inversion

- `-/-` tail-to-tail inversion

- `-/+` tandem duplication

- `+/-` deletion

- different chromosomes — translocation

Each class draws on its own vertical tier with a class-specific color
and a guide line + label per tier.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqPlotR::SeqArch`
-\> `SeqRecon`

## Public fields

- `col_h2h`:

  Color for `+/+` head-to-head inversions.

- `col_t2t`:

  Color for `-/-` tail-to-tail inversions.

- `col_dup`:

  Color for `-/+` tandem duplications.

- `col_del`:

  Color for `+/-` deletions.

- `col_trans`:

  Color for cross-chromosome translocations.

- `drawClasses`:

  Class tiers to draw.

- `tierMultipliers`:

  Per-link tier value populated by `prep()`.

- `last_arc_track`:

  Cached panels list of the arc track, used by `draw()` to render tier
  guides and class labels.

## Methods

### Public methods

- [`SeqReconR6$new()`](#method-SeqRecon-new)

- [`SeqReconR6$prep()`](#method-SeqRecon-prep)

- [`SeqReconR6$draw()`](#method-SeqRecon-draw)

- [`SeqReconR6$clone()`](#method-SeqRecon-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)
- [`SeqPlotR::SeqLink$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqReconR6.

#### Usage

    SeqReconR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      drawClasses = c("Inversion", "Dup/Del", "Translocation"),
      ...
    )

#### Arguments

- `data`:

  BEDPE-like `GRanges` or `data.frame`.

- `mapping`:

  `SeqMap`. Required: `x0`, `x1`, `strand0`, `strand1`. For data.frame
  `data`, also required: `chrom0`, `chrom1`.

- `t0, t1`:

  Track identifiers.

- `aesthetics`:

  Optional `SeqAes`. Override defaults via `h2hColor`, `t2tColor`,
  `dupColor`, `delColor`, `transColor`.

- `drawClasses`:

  Character vector of tiers to render.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Resolve both anchors, classify each link by strand and chromosome, set
per-link tier height + color, then delegate the coordinate work to
`SeqArchR6$prep()`.

#### Usage

    SeqReconR6$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Fallback identifier for `t0`/`t1`.

------------------------------------------------------------------------

### Method `draw()`

Draw class tier guides + labels, then the arches via `SeqArchR6$draw()`
(only for the `Inversion` tier so the labels are layered correctly above
the arches, matching THEfunc `SeqRecon`).

#### Usage

    SeqReconR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqReconR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
