# Merge a track's flat theme over the plot's flat theme

Last-write-wins semantics. Track keys override plot keys for identical
dotted paths.

## Usage

``` r
.merge_themes(plot_theme, track_theme)
```

## Arguments

- plot_theme:

  Flattened plot-level theme.

- track_theme:

  Flattened track-level theme.

## Value

Merged flat theme.
