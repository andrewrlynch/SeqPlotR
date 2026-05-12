# SeqReads R6 class

SeqReads R6 class

SeqReads R6 class

## Details

Internal R6 generator backing
[`seq_reads()`](http://andrewlynch.io/SeqPlotR/reference/seq_reads.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqReads`

## Public fields

- `bam`:

  Path to the BAM file.

- `region`:

  Sorted GRanges of windows to load.

- `gr_reads`:

  GRanges of loaded reads with `qname` and `y` mcols.

- `link_spans`:

  GRanges of mate-pair spans with `y` mcols.

- `nrows`:

  Maximum row index used for read packing.

- `min_mapq`:

  Minimum mapping quality threshold.

- `max_reads`:

  Maximum number of reads to load.

- `max_width`:

  Guardrail on per-window width in bp.

- `sort_by`:

  Row-packing sort order (`"insert_length"` or `"start"`).

- `link_mates`:

  Whether to draw mate-pair link lines.

- `show_strand`:

  Whether to render strand chevrons.

## Methods

### Public methods

- [`SeqReadsR6$new()`](#method-SeqReads-new)

- [`SeqReadsR6$getYLimits()`](#method-SeqReads-getYLimits)

- [`SeqReadsR6$.infer_scale_y()`](#method-SeqReads-.infer_scale_y)

- [`SeqReadsR6$prep()`](#method-SeqReads-prep)

- [`SeqReadsR6$draw()`](#method-SeqReads-draw)

- [`SeqReadsR6$clone()`](#method-SeqReads-clone)

Inherited methods

- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqReadsR6 from a BAM file and region.

#### Usage

    SeqReadsR6$new(
      bam,
      region,
      min_mapq = 0L,
      max_reads = 20000L,
      max_width = 100000L,
      sort_by = c("insert_length", "start"),
      link_mates = TRUE,
      show_strand = TRUE,
      aesthetics = list(),
      ...
    )

#### Arguments

- `bam`:

  Character. Path to an indexed BAM file.

- `region`:

  `GRanges` of windows to load.

- `min_mapq`:

  Integer. Minimum mapping quality. Default `0`.

- `max_reads`:

  Integer. Maximum total reads to load. Default `20000`.

- `max_width`:

  Integer. Maximum allowed window width in bp. Default `100000`.

- `sort_by`:

  Character. `"insert_length"` (default) or `"start"`.

- `link_mates`:

  Logical. Default `TRUE`.

- `show_strand`:

  Logical. Default `TRUE`.

- `aesthetics`:

  Named list of optional aesthetics. See
  [`seq_reads()`](http://andrewlynch.io/SeqPlotR/reference/seq_reads.md).

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `getYLimits()`

Return the y-axis limits for the row-packed reads.

#### Usage

    SeqReadsR6$getYLimits()

------------------------------------------------------------------------

### Method `.infer_scale_y()`

Infer a continuous y scale covering the row-packed reads. Used by
`seq_plot$layoutGrid()` to pick a primary y scale when the parent track
does not define one.

#### Usage

    SeqReadsR6$.infer_scale_y()

------------------------------------------------------------------------

### Method `prep()`

Convert reads and link spans to canvas NPC coordinates, one entry per
window.

#### Usage

    SeqReadsR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Render reads (as chevron polygons) and mate-pair links.

#### Usage

    SeqReadsR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqReadsR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
