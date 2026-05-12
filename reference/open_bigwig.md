# Open a bigWig file connection

Creates a lightweight connection object that probes the bigWig header
(sequence info only — no data) and exposes a `$fetch(region)` method
returning a `GRanges` of signal values restricted to the requested
genomic span.

## Usage

``` r
open_bigwig(path, max_fetch_bp = 50000000L)
```

## Arguments

- path:

  Character. Path to a `.bw` or `.bigwig` file.

- max_fetch_bp:

  Integer. Maximum genomic span allowed per `$fetch()` call. Default
  `5e7` (50 Mb). Prevents accidental whole-genome loads.

## Value

A `SeqBigWig` S3 object with fields `path`, `seqnames`, `seqinfo`,
`max_fetch_bp`, and a `fetch(region)` closure.
