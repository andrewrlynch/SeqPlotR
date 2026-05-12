# Open a `.hic` file connection (Juicer format)

Creates a lightweight connection object that probes the Juicer header
(chromosomes and resolutions) and exposes a `$fetch(region, ...)` method
returning a `data.frame` of contact pairs. The `region` argument may
contain multiple genomic ranges; each range yields its own intra-region
contact submatrix and the results are concatenated.

## Usage

``` r
open_hic(path, resolution = NULL, max_fetch_bp = 280000000L)
```

## Arguments

- path:

  Character. Path to a Juicer `.hic` file.

- resolution:

  Integer. Default bin resolution in bp. May be omitted at construction
  and supplied to `$fetch()`.

- max_fetch_bp:

  Integer. Maximum genomic span per range in a single `$fetch()` call.
  Default `2.8e8` (280 Mb). Each range in the input GRanges is checked
  independently.

## Value

A `SeqHic` S3 object with `chromosomes`, `resolutions`, `max_fetch_bp`,
and a `fetch()` closure.
