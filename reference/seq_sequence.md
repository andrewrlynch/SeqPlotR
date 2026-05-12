# IGV-style nucleotide sequence track element

Renders per-base coloured rectangles (and optionally letters) for
genomic windows up to 200 bp wide. Wider windows emit a message and show
nothing.

## Usage

``` r
seq_sequence(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  genome = NULL,
  sequence = NULL,
  show_letters = FALSE,
  rect_height = NULL,
  colors = NULL,
  ...
)
```

## Arguments

- data:

  Ignored; sequence is fetched from `genome` or `sequence`.

- mapping:

  Ignored.

- aesthetics:

  Optional `SeqAes`. Supports `color` (letter color, defaults to
  matching the rectangle fill), `background` (letter background
  rectangle color, default `NA`).

- genome:

  Character. BSgenome package name (e.g.
  `"BSgenome.Hsapiens.UCSC.hg38"`). Required when `sequence` is `NULL`.

- sequence:

  Character string of nucleotides spanning the first window. When
  provided, `genome` is ignored.

- show_letters:

  Logical. Show nucleotide letters when window width is \<= 80 bp.
  Default `FALSE`.

- rect_height:

  Numeric in (0, 1\] or `NULL`. Fraction of track height used for the
  rectangles. Default `NULL` (full track height).

- colors:

  Named character vector mapping nucleotide codes to hex colors.
  Defaults to UCSC standard: `A="#00AA00"`, `T="#FF0000"`,
  `C="#0000FF"`, `G="#FFB300"`.

- ...:

  Reserved.

## Value

A `SeqSequenceR6` instance.

## Examples

``` r
# From BSgenome (requires BSgenome package)
if (FALSE) { # \dontrun{
  seq_sequence(genome = "BSgenome.Hsapiens.UCSC.hg38")
} # }

# From a string
seq_sequence(sequence = "ATCGATCGATCG", show_letters = TRUE)
#> <SeqSequence>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     colors: #00AA00 #FF0000 #0000FF #FFB300 #AAAAAA
#>     coordCanvas: NULL
#>     data: NULL
#>     draw: function () 
#>     genome: NULL
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), genome = NULL, 
#>     legend: NULL
#>     mapping: NULL
#>     prep: function (layout_track, track_windows) 
#>     rect_height: NULL
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     sequence: ATCGATCGATCG
#>     show_legend: TRUE
#>     show_letters: TRUE
```
