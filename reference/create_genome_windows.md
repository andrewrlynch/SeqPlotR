# Build genome windows from region strings

Parses a vector of region strings — either `"chr:start-end"` ranges or
bare chromosome names — into a `GRanges` of whole-chromosome or
sub-chromosomal windows. Bare chromosome names expand to the full length
of that chromosome in the supplied genome (via
[`default_genome_windows()`](http://andrewlynch.io/SeqPlotR/reference/default_genome_windows.md)
when `genome = "hg38"`). Overlapping / adjacent regions are merged via
[`GenomicRanges::reduce()`](https://rdrr.io/pkg/GenomicRanges/man/inter-range-methods.html)
before return.

## Usage

``` r
create_genome_windows(regions, padding = 0, genome = "hg38", add_chr = TRUE)
```

## Arguments

- regions:

  Character vector of region strings. Each element is either
  `"chr1:1000-2000"` / `"1:1,000-2,000"` (commas tolerated) or a bare
  chromosome name (`"chr1"`, `"X"`, etc.). The `"chr"` prefix is
  optional on input and controlled on output by `add_chr`.

- padding:

  Integer. Base pairs to extend on each side of every region. Clipped so
  `start >= 1`. Default `0`.

- genome:

  Reference genome for bare-chromosome expansion. Only `"hg38"` is
  supported.

- add_chr:

  Logical. If `TRUE` (default), output seqnames carry the `"chr"`
  prefix.

## Value

A sorted, non-overlapping `GRanges` of windows.

## Examples

``` r
create_genome_windows(c("chr1", "chr2:1000000-2000000"))
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames          ranges strand
#>          <Rle>       <IRanges>  <Rle>
#>   [1]     chr2 1000000-2000000      *
#>   [2]     chr1     1-248956422      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
create_genome_windows("chr3", padding = 5e5)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames      ranges strand
#>          <Rle>   <IRanges>  <Rle>
#>   [1]     chr3 1-198795559      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
