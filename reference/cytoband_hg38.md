# UCSC hg38 cytoband table

Cytogenetic band annotations for the hg38 human reference assembly,
pulled from the UCSC cytoBand track. Used by
[`seq_ideogram()`](http://andrewlynch.io/SeqPlotR/reference/seq_ideogram.md)
and
[`load_cytobands()`](http://andrewlynch.io/SeqPlotR/reference/load_cytobands.md).

## Usage

``` r
cytoband_hg38
```

## Format

A data frame with 1,549 rows and 5 columns:

- chrom:

  Chromosome name (UCSC-style, with `"chr"` prefix).

- chromStart:

  0-based start coordinate.

- chromEnd:

  End coordinate (half-open).

- name:

  Cytogenetic band name (e.g. `"p36.33"`).

- gieStain:

  Giemsa stain intensity code (`"gneg"`, `"gpos25"`, `"gpos50"`,
  `"gpos75"`, `"gpos100"`, `"acen"`, `"stalk"`, `"gvar"`).

## Source

UCSC Genome Browser —
<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz>
