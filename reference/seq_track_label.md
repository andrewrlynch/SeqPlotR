# Track label with an optional second line

Draws a text label beside a track's window area, vertically centred on
the track, with an optional smaller second line beneath it. Useful for
naming rows in multi-track figures (samples, haplotypes, assembled
sequences) without relying on axis titles.

## Usage

``` r
seq_track_label(
  label,
  sublabel = NULL,
  side = c("left", "right"),
  aesthetics = aes()
)
```

## Arguments

- label:

  Label text (a single string). `NA` or `""` draws nothing.

- sublabel:

  Optional second line (a single string). `NULL` or `NA` draws no second
  line.

- side:

  `"left"` (default) places the label right-aligned just left of the
  window area; `"right"` places it left-aligned just right of it.

- aesthetics:

  Optional `SeqAes`. Keys: `size` (cex, default 0.75), `col` (default
  `"#1C1B1A"`), `offset` (gap from the window area in lines, default
  0.6), `sublabel.size` (default 0.65), `sublabel.col` (default
  `"grey35"`), `sublabel.face` (default `"italic"`).

## Value

A `SeqTrackLabelR6` instance.

## Examples

``` r
seq_track_label("sample 1", sublabel = "haplotype A")
#> <SeqTrackLabel>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: NULL
#>     draw: function () 
#>     initialize: function (label, sublabel = NULL, side = c("left", "right"), 
#>     label: sample 1
#>     legend: NULL
#>     mapping: NULL
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: FALSE
#>     side: left
#>     sublabel: haplotype A
#>     y_fields: y y_end y_min y_max y0 y1 height
#>     y_transform: NULL
```
