# Caption and banner under the figure

Adds a caption line centred under the whole figure and an optional
banner line beneath it (for example a warning that the figure is
approximate). Both strings are drawn exactly as given. The bottom canvas
margin grows by a fixed `allowance` (in lines) whether or not a banner
is given, so track heights never depend on the caption content.

## Usage

``` r
seq_caption(text, banner = NULL, allowance = 2.8, aesthetics = aes())
```

## Arguments

- text:

  Caption text (a single string). `NULL`, `NA` or `""` draws nothing.

- banner:

  Optional banner text, drawn one line under the caption.

- allowance:

  Bottom space reserved, in lines (default `2.8`).

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `caption.size` (cex, `0.6`), `caption.col` (`"grey30"`), `banner.size`
  (`0.7`), `banner.col` (`"#D33A2C"`), `banner.face` (`"bold"`).

## Value

A `SeqCaptionR6` instance.

## Details

Add to a plot with `%+%`; a caption cannot be added to a track.

## Examples

``` r
w <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1e6))
p <- seq_plot() %+% seq_track(windows = w) %+%
  seq_caption("Model output, not truth", banner = "APPROXIMATE")
```
