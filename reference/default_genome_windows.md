# Default whole-genome windows

Returns a `GRanges` spanning all hg38 autosomes plus X and Y, one range
per chromosome, covering the full chromosome length. Useful as a default
`windows` argument for whole-genome ideogram or summary tracks.

## Usage

``` r
default_genome_windows(add_chr = TRUE)
```

## Arguments

- add_chr:

  Logical. If `TRUE` (default), seqnames carry the `"chr"` prefix
  (`"chr1"`, `"chr2"`, ...). If `FALSE`, the prefix is stripped (`"1"`,
  `"2"`, ...).

## Value

A `GRanges` of length 24 (chr1–22, X, Y).

## Examples

``` r
default_genome_windows()
#> GRanges object with 24 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]     chr1 1-248956422      *
#>    [2]     chr2 1-242193529      *
#>    [3]     chr3 1-198295559      *
#>    [4]     chr4 1-190214555      *
#>    [5]     chr5 1-181538259      *
#>    ...      ...         ...    ...
#>   [20]    chr20  1-64444167      *
#>   [21]    chr21  1-46709983      *
#>   [22]    chr22  1-50818468      *
#>   [23]     chrX 1-156040895      *
#>   [24]     chrY  1-57227415      *
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
default_genome_windows(add_chr = FALSE)
#> GRanges object with 24 ranges and 0 metadata columns:
#>        seqnames      ranges strand
#>           <Rle>   <IRanges>  <Rle>
#>    [1]        1 1-248956422      *
#>    [2]        2 1-242193529      *
#>    [3]        3 1-198295559      *
#>    [4]        4 1-190214555      *
#>    [5]        5 1-181538259      *
#>    ...      ...         ...    ...
#>   [20]       20  1-64444167      *
#>   [21]       21  1-46709983      *
#>   [22]       22  1-50818468      *
#>   [23]        X 1-156040895      *
#>   [24]        Y  1-57227415      *
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```
