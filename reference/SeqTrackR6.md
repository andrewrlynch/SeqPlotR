# SeqTrack R6 class

SeqTrack R6 class

SeqTrack R6 class

## Details

Internal R6 class backing
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).
Public users should go through the snake_case constructor; this class is
exported only so that other internal files can reference it.

## Public fields

- `data`:

  GRanges default data for elements in this track.

- `mapping`:

  SeqMap default mapping for elements in this track.

- `aesthetics`:

  SeqAes default constant aesthetics.

- `windows`:

  GRanges windows for this track (only place windows are set).

- `track_id`:

  Character unique identifier for this track.

- `track_width`:

  Relative width unit within its row.

- `track_height`:

  Relative height unit within its column.

- `direction`:

  One of "right" or "under".

- `elements`:

  List of SeqElement / SeqLink added via `%+%`.

- `scale_x`:

  Optional position scale for the primary x axis.

- `scale_y`:

  Optional position scale for the primary y axis.

- `scale_x2`:

  Optional position scale for the secondary x axis. Populated
  automatically at `layoutGrid()` time from elements whose
  `map(axis.x = 2)` targets the secondary x axis.

- `scale_y2`:

  Optional position scale for the secondary y axis.

- `uses_genomic_y`:

  TRUE when mapping\$y resolves to a genomic special.

- `y_windows`:

  Optional GRanges windows for the primary y axis.

- `y_windows2`:

  Optional GRanges windows for the secondary y axis.

- `has_axis_x2`:

  TRUE when any element targets axis.x = 2, when scale_x2 is explicitly
  set, or when the theme forces visibility. Set at layoutGrid() time.

- `has_axis_y2`:

  Symmetric to `has_axis_x2` for the y direction.

- `resolved_theme`:

  Populated at layoutGrid(): nested list with per-axis specs
  (x1/x2/y1/y2), track chrome, and per-window flag.

- `track_outer_margin`:

  Named list of track-level outer margin (npc) — reserves the band where
  axis titles sit.

- `track_inner_margin`:

  Named list of track-level inner margin (npc) — separates axis titles
  from the window row.

- `window_outer_margin`:

  Named list of per-window outer margin (npc).

- `window_inner_margin`:

  Named list of per-window inner margin (npc) — holds per-window axis
  ticks and labels.

- `window_margin`:

  Deprecated. Set `aes("window.gap.width" = <value>)` on the plot or
  track instead. Stored for backward compatibility but has no effect on
  layout.

- `combine_windows`:

  Logical. When `TRUE`, multi-region windows are concatenated into a
  single virtual panel, with per-original- window axis labels and a thin
  separator at each window boundary. Used to draw cross-window data
  (e.g. inter-chromosomal Hi-C contacts) within one continuous track.

- `combine_y_windows`:

  Logical. Symmetric to `combine_windows` for the y-axis (genomic y
  tracks only).

- `flip_x`:

  Logical. When `TRUE`, mirror the x-axis: low data values render at the
  right edge of the panel and high values at the left. Tick labels
  follow the same orientation.

- `flip_y`:

  Logical. Symmetric to `flip_x` for the y-axis. For Hi-C `triangle`
  style this produces a downward-pointing triangle; for `diagonal` it
  shows the lower diagonal.

- `show_legend`:

  Logical. When `FALSE`, this track contributes no legend keys
  regardless of the `legend` fields on its elements. Default `TRUE`.

- `window_scale`:

  Numeric vector or `NULL`. Per-window x-axis scale factors (e.g. `1e-6`
  for Mb, `1e-3` for kb, `1` for bp). When `NULL` (default), the scale
  is inferred from the narrowest window. When length 1, the value is
  applied to all windows. When length equals the number of windows,
  values are applied positionally. Any other length triggers a warning
  and recycles with [`rep_len()`](https://rdrr.io/r/base/rep.html).

## Methods

### Public methods

- [`SeqTrackR6$new()`](#method-SeqTrack-new)

- [`SeqTrackR6$addElement()`](#method-SeqTrack-addElement)

- [`SeqTrackR6$collect_legend_keys()`](#method-SeqTrack-collect_legend_keys)

- [`SeqTrackR6$clone()`](#method-SeqTrack-clone)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqTrackR6.

#### Usage

    SeqTrackR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      windows = NULL,
      track_id = NULL,
      direction = "right",
      track_width = 1,
      track_height = 1,
      scale_x = NULL,
      scale_y = NULL,
      scale_x2 = NULL,
      scale_y2 = NULL,
      y_windows = NULL,
      y_windows2 = NULL,
      track_outer_margin = 0.02,
      track_inner_margin = 0.02,
      window_outer_margin = 0,
      window_inner_margin = 0.02,
      window_margin = NULL,
      combine_windows = FALSE,
      combine_y_windows = FALSE,
      flip_x = FALSE,
      flip_y = FALSE,
      elements = list(),
      show_legend = TRUE,
      window_scale = NULL,
      ...
    )

#### Arguments

- `data`:

  GRanges default data.

- `mapping`:

  SeqMap default mapping.

- `aesthetics`:

  SeqAes default constant aesthetics.

- `windows`:

  GRanges windows for this track.

- `track_id`:

  Character unique identifier.

- `direction`:

  One of "right" or "under".

- `track_width`:

  Relative width unit.

- `track_height`:

  Relative height unit.

- `scale_x`:

  Optional primary x position scale.

- `scale_y`:

  Optional primary y position scale.

- `scale_x2`:

  Optional secondary x position scale.

- `scale_y2`:

  Optional secondary y position scale.

- `y_windows`:

  Optional GRanges y-axis windows for the primary axis.

- `y_windows2`:

  Optional GRanges y-axis windows for the secondary axis.

- `track_outer_margin`:

  Scalar, length-4 `c(bottom, left, top, right)`, or named list. Default
  `0.02`. Band where axis titles (derived from the track `mapping`) are
  drawn.

- `track_inner_margin`:

  Same form. Default `0.02`. Separates the title zone from the window
  row.

- `window_outer_margin`:

  Same form. Default `0`. Spacer around each window inside the track
  plot region.

- `window_inner_margin`:

  Same form. Default `0.02`. Holds each window's axis ticks and labels.

- `window_margin`:

  Deprecated. Use `aes("window.gap.width" = <value>)` on the plot or
  track aesthetics to control inter-window gap. Passing a non-`NULL`
  value emits a deprecation warning and is otherwise ignored.

- `elements`:

  Optional list of `SeqElement` / `SeqLink` objects to pre-populate this
  track. Elements can also be added later via `addElement()` or the
  `%+%` operator.

- `show_legend`:

  Logical. When `FALSE`, this track contributes no legend keys
  regardless of element `legend` fields. Default `TRUE`.

- `window_scale`:

  Numeric vector or `NULL`. Per-window x-axis scale factors (e.g. `1e-6`
  for Mb, `1e-3` for kb, `1` for bp). When `NULL` (default), the scale
  is inferred from the narrowest window. Length 1 applies the value to
  all windows; length equal to the number of windows uses values
  positionally; any other length triggers a warning and recycles with
  [`rep_len()`](https://rdrr.io/r/base/rep.html).

- `...`:

  Additional arguments (currently ignored).

------------------------------------------------------------------------

### Method `addElement()`

Append an element (SeqElement or SeqLink) to this track.

#### Usage

    SeqTrackR6$addElement(elem)

#### Arguments

- `elem`:

  The element to add.

------------------------------------------------------------------------

### Method `collect_legend_keys()`

Collect legend keys from all elements in this track.

Iterates every element (including `SeqLink` subclasses) and calls
`$collect_legend_keys()` on each. Returns a flat list of entries, where
each entry is a named list with fields `title`, `key`, and
`element_class` (as produced by `SeqElement$collect_legend_keys()`).
Returns `NULL` when `show_legend` is `FALSE` or no element contributes
any keys.

#### Usage

    SeqTrackR6$collect_legend_keys()

#### Returns

A list of legend-key entries, or `NULL`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqTrackR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
