# Open a .cool or .mcool HDF5 contact matrix file

Supports both `.cool` (single-resolution) and `.mcool`
(multi-resolution) files produced by cooler / HiCExplorer.

## Usage

``` r
open_h5(path, resolution = NULL, max_fetch_bp = 5000000L)
```

## Arguments

- path:

  Character. Path to the `.cool` or `.mcool` file.

- resolution:

  Integer or `NULL`. For `.mcool` files, which resolution to open. When
  `NULL`, the coarsest available resolution is used for the probe; you
  must set `resolution` explicitly in `$fetch()`.

- max_fetch_bp:

  Integer. Maximum genomic span per `$fetch()` call. Default `5e6` (5
  Mb).

## Value

A `SeqH5` S3 object with `$fetch()` and `$fetch_binned()` methods.
