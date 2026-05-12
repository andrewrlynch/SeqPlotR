# Test whether `x` is a SeqPlotR file connection

Connection objects from
[`open_bigwig()`](http://andrewlynch.io/SeqPlotR/reference/open_bigwig.md),
[`open_bam()`](http://andrewlynch.io/SeqPlotR/reference/open_bam.md),
[`open_hic()`](http://andrewlynch.io/SeqPlotR/reference/open_hic.md),
and [`open_h5()`](http://andrewlynch.io/SeqPlotR/reference/open_h5.md)
all carry the `SeqFileConn` class.

## Usage

``` r
is_seq_file_conn(x)
```

## Arguments

- x:

  An object.

## Value

`TRUE` if `x` inherits from `SeqFileConn`, otherwise `FALSE`.
