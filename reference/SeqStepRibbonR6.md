# SeqStepRibbon R6 class

R6 generator backing
[`seq_step_ribbon()`](http://andrewlynch.io/SeqPlotR/reference/seq_step_ribbon.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqElement` -\> `SeqStepRibbon`

## Public fields

- `y_fields`:

  Resolved fields in y axis units.

- `legend_arg`:

  The `legend` argument as given to the constructor.

- `legend_labels`:

  Named character relabelling auto keys by thickness.

- `legend_extra`:

  List of `LegendKey` appended to the auto keys.

- `legend_position`:

  Position of the element legend.

- `legend_side`:

  Side of the element legend.

## Methods

### Public methods

- [`SeqStepRibbon$new()`](#method-SeqStepRibbon-initialize)

- [`SeqStepRibbon$prep()`](#method-SeqStepRibbon-prep)

- [`SeqStepRibbon$draw()`](#method-SeqStepRibbon-draw)

- [`SeqStepRibbon$clone()`](#method-SeqStepRibbon-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$resolve()`

------------------------------------------------------------------------

### `SeqStepRibbon$new()`

Construct a SeqStepRibbonR6.

#### Usage

    SeqStepRibbon$new(
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

  Named character relabelling auto keys.

- `legend_extra`:

  List of `LegendKey` appended to the legend keys.

- `legend_position, legend_side`:

  Placement of the element legend.

- `show_legend`:

  Logical.

- `...`:

  Unused.

------------------------------------------------------------------------

### `SeqStepRibbon$prep()`

Resolve mappings, clip spans to windows and fill `coordCanvas` with one
data.frame per window.

#### Usage

    SeqStepRibbon$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

#### Returns

`coordCanvas`, invisibly.

------------------------------------------------------------------------

### `SeqStepRibbon$draw()`

Draw rects, interval lines, triangles and labels.

#### Usage

    SeqStepRibbon$draw()

------------------------------------------------------------------------

### `SeqStepRibbon$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqStepRibbon$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
