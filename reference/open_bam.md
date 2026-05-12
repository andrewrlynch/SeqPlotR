# Open a BAM file connection

Creates a lightweight connection object that validates the BAM index,
probes the header for sequence names and lengths, and exposes a
`$fetch(region, ...)` method returning a `GRanges` of alignments.

## Usage

``` r
open_bam(path, max_fetch_bp = 100000L)
```

## Arguments

- path:

  Character. Path to an indexed BAM file (`.bam`).

- max_fetch_bp:

  Integer. Maximum genomic span per `$fetch()` call. Default `1e5` (100
  kb). BAM fetches over wide regions can be very slow — keep this tight.

## Value

A `SeqBam` S3 object with `seqnames`, `seq_lengths`, `max_fetch_bp`, and
a `fetch(region, min_mapq, max_reads)` closure.
