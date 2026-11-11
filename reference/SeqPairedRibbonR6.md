# SeqPairedRibbon R6 class

R6 generator backing
[`seq_paired_ribbon()`](http://andrewlynch.io/SeqPlotR/reference/seq_paired_ribbon.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqElement` -\> `SeqPairedRibbon`

## Public fields

- `y_fields`:

  Resolved fields in y axis units.

- `legend_arg`:

  The `legend` argument as given to the constructor.

- `legend_labels`:

  Character relabelling the two auto keys.

- `legend_extra`:

  List of `LegendKey` appended to the auto keys.

- `legend_position`:

  Position of the element legend.

- `legend_side`:

  Side of the element legend.

## Methods

### Public methods

- [`SeqPairedRibbon$new()`](#method-SeqPairedRibbon-initialize)

- [`SeqPairedRibbon$prep()`](#method-SeqPairedRibbon-prep)

- [`SeqPairedRibbon$draw()`](#method-SeqPairedRibbon-draw)

- [`SeqPairedRibbon$clone()`](#method-SeqPairedRibbon-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$resolve()`

------------------------------------------------------------------------

### `SeqPairedRibbon$new()`

Construct a SeqPairedRibbonR6.

#### Usage

    SeqPairedRibbon$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      legend = NULL,
      legend_labels = NULL,
      legend_extra = NULL,
      legend_position = "canvas_margin",
      legend_side = "right",
      show_legend = TRUE,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges`.

- `mapping`:

  Optional `SeqMap`.

- `aesthetics`:

  Optional `SeqAes`.

- `legend`:

  `NULL` (auto keys), `FALSE`, a `SeqLegendSpec`, or a list of
  `LegendKey`.

- `legend_labels`:

  Character relabelling the two auto keys.

- `legend_extra`:

  List of `LegendKey` appended to the legend keys.

- `legend_position, legend_side`:

  Placement of the element legend.

- `show_legend`:

  Logical.

- `...`:

  Unused.

------------------------------------------------------------------------

### `SeqPairedRibbon$prep()`

Resolve mappings, clip spans to windows and fill `coordCanvas` with one
data.frame per window.

#### Usage

    SeqPairedRibbon$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

#### Returns

`coordCanvas`, invisibly.

------------------------------------------------------------------------

### `SeqPairedRibbon$draw()`

Draw all second-value rects, then all first-value rects.

#### Usage

    SeqPairedRibbon$draw()

------------------------------------------------------------------------

### `SeqPairedRibbon$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqPairedRibbon$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
