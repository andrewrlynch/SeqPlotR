# Draw gene models

Composite element rendering gene models with backbone lines, exon/UTR
boxes, directional arrows, and labels. Entirely format-agnostic — supply
the relevant column names via
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).

## Usage

``` r
seq_gene(
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
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Recognised fields: `group` (feature-to-gene
  grouping; features sharing a value are one gene), `strand` (default
  `"+"`), `label` (defaults to the group value), `type` (`"exon"`
  full-height box, `"UTR"` 80-percent-height box, anything else no box),
  and `color` (per-feature color).

- aesthetics:

  Optional `SeqAes`. Supports `color` (default `"gray30"`), `linewidth`,
  `alpha`.

- backbone_type:

  Character. One of `"arrow"` (default), `"solid"`, or `"dashed"`.
  Controls backbone line style and whether chevron arrows are drawn.
  Ignored when `style_type` is `"gene"` or `"point"` (no backbone exists
  in those modes).

- show_start:

  Logical. When `TRUE`, a TSS flag arrow is drawn above the first exon
  of each gene. Default `FALSE`. Ignored when `style_type = "point"`
  (the point itself marks the TSS).

- tss_position:

  Named list overriding per-gene TSS genomic positions. Keys are gene
  IDs; values are `c(start, end)` genomic coordinates. Default `NULL`
  (auto-detected from first exon by row order).

- separate_strands:

  Logical. When `TRUE`, genes are split into `"+"` (top) and `"-"`
  (bottom) sub-bands, each labelled. Silently ignored when only one
  strand is present. Default `FALSE`.

- style_type:

  Character. One of `"exon"` (default), `"gene"`, or `"point"`. Selects
  the per-gene rendering style: `"exon"` draws the full backbone with
  exon/UTR boxes; `"gene"` draws a single chevron-shaped polygon
  spanning the gene extent (no exon detail); `"point"` draws a single
  filled circle at the TSS. Labels are drawn in all three modes.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqGeneR6` instance.

## Examples

``` r
seq_gene(map(group = gene_id, strand = strand, label = gene_name))
#> <SeqGene>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     backbone_type: arrow
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: NULL
#>     draw: function () 
#>     exon_height: 0.8
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), backbone_type = "arrow", 
#>     label_cex: 0.6
#>     label_offset: 0.01
#>     label_pad: 50000
#>     legend: NULL
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     separate_strands: FALSE
#>     show_legend: TRUE
#>     show_start: FALSE
#>     style_type: exon
#>     tss_position: NULL
```
