# Giemsa stain code → fill color

Maps UCSC `gieStain` codes (`gneg`, `gpos25`..`gpos100`, `acen`,
`stalk`, `gvar`) to fill colors. Unknown codes fall through to
`"#CCCCCC"`.

## Usage

``` r
.ideogram_fill_colors(stain)
```

## Arguments

- stain:

  Character vector of `gieStain` codes.

## Value

Character vector of fill colors, same length as `stain`.
