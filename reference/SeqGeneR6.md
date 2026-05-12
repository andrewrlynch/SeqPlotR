# SeqGene R6 class

SeqGene R6 class

SeqGene R6 class

## Details

Internal R6 generator backing
[`seq_gene()`](http://andrewlynch.io/SeqPlotR/reference/seq_gene.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqGene`

## Public fields

- `exon_height`:

  Proportion of per-tier height allocated to the exon box.

- `label_pad`:

  Minimum padding (in bp) around each gene's extent to keep labels from
  overlapping neighbouring genes.

- `label_cex`:

  Expansion factor for gene label text.

- `label_offset`:

  Horizontal offset (in npc units) applied when placing labels flush
  against the gene backbone.

- `backbone_type`:

  Character. One of `"arrow"`, `"solid"`, `"dashed"`. Controls how the
  gene backbone line is rendered. Default `"arrow"`.

- `show_start`:

  Logical. When `TRUE`, a TSS flag arrow is drawn above the first exon
  (by row order) of each gene. Default `FALSE`.

- `tss_position`:

  Named list of length-2 numeric vectors. Keys are gene IDs (matching
  `group` values); values are `c(start, end)` genomic coordinates
  overriding the auto-detected first-exon TSS position. Default `NULL`
  (auto from first exon by row order).

- `separate_strands`:

  Logical. When `TRUE`, genes are placed in two horizontal sub-bands by
  strand (`"+"` top, `"-"` bottom), each labelled. Silently ignored when
  all strands are `"*"` or only one strand is present. Default `FALSE`.

- `style_type`:

  Character. One of `"exon"`, `"gene"`, `"point"`. Selects the per-gene
  rendering style: `"exon"` draws backbone + exon/UTR boxes (default);
  `"gene"` draws a single chevron-shaped polygon spanning the full gene
  extent; `"point"` draws a single filled circle at the TSS. Default
  `"exon"`.

## Methods

### Public methods

- [`SeqGeneR6$new()`](#method-SeqGene-new)

- [`SeqGeneR6$prep()`](#method-SeqGene-prep)

- [`SeqGeneR6$draw()`](#method-SeqGene-draw)

- [`SeqGeneR6$clone()`](#method-SeqGene-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqGeneR6.

#### Usage

    SeqGeneR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      backbone_type = "arrow",
      show_start = FALSE,
      tss_position = NULL,
      separate_strands = FALSE,
      style_type = c("exon", "gene", "point"),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`. Recognised: `group`, `strand`, `label`, `type`,
  `color`.

- `aesthetics`:

  Optional `SeqAes`. Supports `color` (default `"gray30"` when no
  `color` mapping is given), `linewidth`, `alpha`.

- `backbone_type`:

  Character. One of `"arrow"`, `"solid"`, `"dashed"`. Default `"arrow"`.
  Ignored when `style_type` is `"gene"` or `"point"`.

- `show_start`:

  Logical. Draw TSS flag arrow above first exon. Default `FALSE`.
  Ignored when `style_type = "point"` (the point itself marks the TSS).

- `tss_position`:

  Named list overriding per-gene TSS genomic positions. Default `NULL`.

- `separate_strands`:

  Logical. Split track into `"+"` and `"-"` sub-bands. Default `FALSE`.

- `style_type`:

  Character. One of `"exon"` (default), `"gene"`, or `"point"`. Controls
  per-gene rendering: full backbone + exon boxes, single chevron
  polygon, or a single TSS point.

- `...`:

  Unused.

------------------------------------------------------------------------

### Method `prep()`

Resolve mappings, stack genes into non-overlapping tiers, and compute
canvas coordinates for backbones, exons/UTRs, arrows, and labels.
Populates `self$coordCanvas` with a single `data.frame` whose rows are
per-feature draw primitives.

#### Usage

    SeqGeneR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw backbones, exons/UTRs, directional arrows, and labels.

#### Usage

    SeqGeneR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqGeneR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
