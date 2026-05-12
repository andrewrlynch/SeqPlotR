# SeqScalebar R6 class

SeqScalebar R6 class

SeqScalebar R6 class

## Details

Internal R6 generator backing
[`seq_scalebar()`](http://andrewlynch.io/SeqPlotR/reference/seq_scalebar.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqPlotR::SeqElement` -\> `SeqScalebar`

## Public fields

- `length_bp`:

  Length of the bar in base pairs.

- `label`:

  Character label rendered above the bar.

- `hjust`:

  Horizontal anchor of the bar's right edge in panel NPC.

- `vjust`:

  Vertical anchor of the bar in panel NPC.

- `ticks`:

  Whether to draw tick caps at the bar's ends.

- `tick_height`:

  Tick height as a fraction of the panel height.

- `bar_lwd`:

  Line width of the bar and ticks.

- `bar_col`:

  Bar / tick color.

- `label_cex`:

  Label `cex`.

- `label_col`:

  Label color.

- `label_pad`:

  Vertical NPC gap between bar top and label bottom.

## Methods

### Public methods

- [`SeqScalebarR6$new()`](#method-SeqScalebar-new)

- [`SeqScalebarR6$prep()`](#method-SeqScalebar-prep)

- [`SeqScalebarR6$draw()`](#method-SeqScalebar-draw)

- [`SeqScalebarR6$clone()`](#method-SeqScalebar-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqScalebarR6.

#### Usage

    SeqScalebarR6$new(
      length_bp,
      label = NULL,
      hjust = 0.95,
      vjust = 0.05,
      ticks = TRUE,
      tick_height = 0.15,
      bar_lwd = 1.2,
      bar_col = "#1C1B1A",
      label_cex = 0.65,
      label_col = "#1C1B1A",
      label_pad = 0.005,
      ...
    )

#### Arguments

- `length_bp`:

  Numeric, positive.

- `label`:

  Optional character label. When `NULL`, auto-formatted.

- `hjust`:

  Numeric in 0–1. Default `0.95`.

- `vjust`:

  Numeric in 0–1. Default `0.05`.

- `ticks`:

  Logical. Default `TRUE`.

- `tick_height`:

  Numeric. Default `0.15`.

- `bar_lwd`:

  Numeric. Default `1.2`.

- `bar_col`:

  Character. Default `"#1C1B1A"`.

- `label_cex`:

  Numeric. Default `0.65`.

- `label_col`:

  Character. Default `"#1C1B1A"`.

- `label_pad`:

  Numeric. Default `0.005`.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Compute per-window NPC coordinates for the bar, ticks, and label.

#### Usage

    SeqScalebarR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw the bar, tick caps, and label.

#### Usage

    SeqScalebarR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqScalebarR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
