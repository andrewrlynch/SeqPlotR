# Ideogram-slice blocks

Draws one block per range: a rectangle (or an arrow strip) spanning the
exact track positions `[x, x_end]`, filled either with a solid colour or
with the cytobands of a reference interval `chrom:ref_start-ref_end`
mapped into the block, mirrored when the block is flipped. A centromere
interval is drawn as a red bowtie with the outline on the arms only.
Labels sit inside the block when requested and they fit, otherwise
outside (above or below) with a three-level stagger. Useful for synteny,
rearrangement and assembly figures where pieces of reference chromosomes
are laid out along new coordinates.

## Usage

``` r
seq_block(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  fill = c("solid", "cytoband"),
  cytobands = NULL,
  palette = NULL,
  tints = c(0.15, 0.5, 0.3, 0.65, 0.4),
  block_colors = NULL,
  shape = c("rect", "arrow"),
  height = 0.5,
  label_side = c("below", "above"),
  flip_style = c("rotate", "overline", "minus"),
  annotations = NULL,
  legend = NULL,
  show_legend = TRUE,
  ...
)
```

## Arguments

- data:

  `GRanges` whose seqnames are the track windows the blocks sit in and
  whose ranges are the block positions in track coordinates.

- mapping:

  Optional `SeqMap`. Fields (all optional): `x`, `x_end` (exact numeric
  block edges; default `start`, `end`); `chrom`, `ref_start`, `ref_end`
  (reference interval drawn inside the block, half-open; default
  `seqnames`, `x`, `x_end`); `flipped` (mirror the reference interval);
  `label`; `key` (fill key or a colour; default `chrom`); `tint`
  (integer tint index; default 1); `label_inside`; `dashed` (dashed
  outline); `centromere_start`, `centromere_end` (reference coordinates;
  `NA` means none).

- aesthetics:

  Optional `SeqAes`. Keys: `outline.col` (default `"grey25"`),
  `outline.lwd` (0.7), `label.size` (cex, 0.7), `label.col` (`"black"`),
  `centromere.fill` (`"#D33A2C"`), `fade` (alpha of a white overlay, 0),
  `arrow.head` (maximum head length as a fraction of block width, 0.35),
  `arrow.head_max` (maximum head length as a fraction of the track
  width, 0.012).

- fill:

  `"solid"` (default) or `"cytoband"`. Blocks with no matching cytobands
  fall back to solid with a message.

- cytobands:

  Cytoband table: a `GRanges` as returned by
  [`load_cytobands()`](http://andrewlynch.io/SeqPlotR/reference/load_cytobands.md)
  or a data frame with `chrom`, `start` (0-based), `end`, `stain` (UCSC
  `chromStart`, `chromEnd`, `gieStain` also work). `NULL` loads the
  bundled hg38 bands. Chromosome names are matched with or without a
  `chr` prefix.

- palette:

  `NULL` (default hues), `"soft"` (a softer set), or a character vector
  of colours. A named vector is looked up by key; an unnamed one is
  recycled over the distinct keys in data order.

- tints:

  Fractions mixed toward white, indexed by `tint`.

- block_colors:

  Optional named colours overriding the fill by label, then by key (not
  tinted).

- shape:

  `"rect"` (default) or `"arrow"` (solid strips pointing toward
  increasing reference coordinate, reversed when flipped).

- height:

  Block height as a fraction of the panel height. Default 0.5.

- label_side:

  `"below"` (default) or `"above"` for outside labels.

- flip_style:

  How flipped blocks are marked: `"rotate"` (label rotated 180 degrees),
  `"overline"` (line over the label) or `"minus"` (a leading minus
  sign).

- annotations:

  Optional data frame of point or interval annotations with columns
  `chrom`, `start`, optional `end`, `label`, and optional `side`
  (`"above"`/`"below"`), `color`, `italic`, `dot`, `strand` (`"+"`/`"-"`
  draws a direction arrow). Each row is drawn on every block whose
  reference interval contains its midpoint, mapped through flips.

- legend:

  Optional `LegendKey` list or `SeqLegendSpec`. No legend is generated
  automatically.

- show_legend:

  Logical. Default `TRUE`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqBlockR6` instance.

## Examples

``` r
gr <- GenomicRanges::GRanges("der1", IRanges::IRanges(c(1, 3e7), c(3e7, 6e7)),
                             label = c("A", "B"), flipped = c(FALSE, TRUE))
seq_block(gr, map(label = label, flipped = flipped, chrom = "chr1"))
#> <SeqBlock>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     annotations: NULL
#>     auto_legend: NULL
#>     block_colors: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     cytobands: NULL
#>     data: GRanges
#>     draw: function () 
#>     fill: solid
#>     flip_style: rotate
#>     height: 0.5
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), fill = c("solid", 
#>     label_side: below
#>     legend: NULL
#>     mapping: SeqMap
#>     palette: NULL
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     shape: rect
#>     show_legend: TRUE
#>     tints: 0.15 0.5 0.3 0.65 0.4
#>     y_fields: y y_end y_min y_max y0 y1 height
#>     y_transform: NULL
#>   Private:
#>     bands: NULL
```
