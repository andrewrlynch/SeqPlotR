# IGV-style read alignment track element

Renders read alignments from an indexed BAM file as chevron polygons
(pointing right for + strand reads, left for −). Mate pairs are
optionally connected by a thin horizontal line. Reads are loaded once at
construction time and cached on the element; only the requested windows
are pulled.

## Usage

``` r
seq_reads(
  bam,
  region,
  min_mapq = 0L,
  max_reads = 20000L,
  max_width = 100000L,
  sort_by = c("insert_length", "start"),
  link_mates = TRUE,
  show_strand = TRUE,
  aesthetics = list()
)
```

## Arguments

- bam:

  Character. Path to an indexed BAM file.

- region:

  `GRanges` of windows to load. Every window must be ≤ `max_width` bp.

- min_mapq:

  Integer. Minimum mapping quality. Default `0`.

- max_reads:

  Integer. Maximum total reads to load (qname-aware downsampling
  preserves mate pairs). Default `20000`.

- max_width:

  Integer. Maximum allowed window width in bp. Default `100000`.
  Prevents accidental loading of entire chromosomes.

- sort_by:

  Character. Row-packing sort order: `"insert_length"` (default —
  longest inserts at the top) or `"start"` (leftmost reads at the top).

- link_mates:

  Logical. Draw a horizontal line connecting each mated pair's outer
  extent. Default `TRUE`.

- show_strand:

  Logical. Draw chevron arrowheads indicating strand. Default `TRUE`.

- aesthetics:

  Named list. Supported keys: `fill_plus`, `fill_minus`,
  `fill_unstranded`, `col`, `lwd`, `row_gap`, `tip_mm`,
  `tip_min_body_mm`, `link_col`, `link_lwd`, `link_lty`.

## Value

A `SeqReadsR6` instance.
