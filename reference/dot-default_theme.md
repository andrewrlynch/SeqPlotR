# The built-in default theme for SeqPlotR

Flat dotted-key list combining layout parameters (`track_gaps`,
`window_gaps`, `margins`) with the hierarchical axis and track chrome
keys. Users override via `seq_plot(aesthetics = aes(...))` or
`seq_track(aesthetics = aes(...))`.

## Usage

``` r
.default_theme()
```

## Details

Layout-only keys (`track_gaps`, `window_gaps`, `margins`) are not
subject to axis inheritance — they're looked up by exact name.
