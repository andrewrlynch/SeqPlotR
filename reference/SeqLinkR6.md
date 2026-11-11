# SeqLink R6 base class

R6 generator for cross-track and within-track link elements, exported
for subclassing. Extends
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md)
with track identifiers and a pair of anchor GRanges synthesised from the
resolved [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md)
fields. Subclasses implement
`prep(layout_all_tracks, track_windows_list, plot_track_index)` and
`draw()`;
[`seq_to_canvas()`](http://andrewlynch.io/SeqPlotR/reference/seq_to_canvas.md)
maps data values into any track's panel.

## Super class

`SeqElement` -\> `SeqLink`

## Public fields

- `t0`:

  `track_id` (or integer index) of the first anchor's track.

- `t1`:

  `track_id` (or integer index) of the second anchor's track.

- `anchor0_gr`:

  Synthetic point GRanges for anchor 0, populated by `resolve()` from
  the resolved `chrom0`, `x0`, and `strand0` fields.

- `anchor1_gr`:

  Synthetic point GRanges for anchor 1.

- `anchor1_missing`:

  Logical per row, set by `resolve()`: `TRUE` where `x1` or `chrom1` is
  `NA` and anchor 1 falls back to anchor 0.

## Methods

### Public methods

- [`SeqLink$new()`](#method-SeqLink-initialize)

- [`SeqLink$resolve()`](#method-SeqLink-resolve)

- [`SeqLink$prep()`](#method-SeqLink-prep)

- [`SeqLink$.resolve_track_ref()`](#method-SeqLink-.resolve_track_ref)

- [`SeqLink$clone()`](#method-SeqLink-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$draw()`

------------------------------------------------------------------------

### `SeqLink$new()`

Construct a new SeqLinkR6.

#### Usage

    SeqLink$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      legend = NULL,
      show_legend = TRUE,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame`. A single argument carries both
  anchors; there is no `data2`.

- `mapping`:

  Optional `SeqMap`. Must define `x0` and `x1` (and, for data.frame
  `data`, `chrom0` and `chrom1`).

- `t0, t1`:

  Track identifiers (string or integer index). Locked to the parent
  track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
  via `%+%`.

- `aesthetics`:

  Optional `SeqAes`.

- `legend`:

  Optional `LegendKey` or list of `LegendKey` objects.

- `show_legend`:

  Logical. Set to `FALSE` to suppress legend output.

- `...`:

  Unused — accepted so subclasses can pass extra arguments.

------------------------------------------------------------------------

### `SeqLink$resolve()`

Resolve the mapping against the effective `data`, then synthesise
`anchor0_gr` and `anchor1_gr` from the resolved fields. Errors if `x0`
or `x1` is missing, or if `chrom0`/`chrom1` cannot be derived (only
auto-fillable for GRanges `data`). Rows with an `NA` `x1` or `chrom1`
take anchor 0 as anchor 1 and are flagged in `anchor1_missing`.

#### Usage

    SeqLink$resolve(track_data = NULL, track_mapping = NULL)

#### Arguments

- `track_data`:

  Optional `GRanges`/`data.frame` from the parent track (used as a
  fallback when the link has no own `data`).

- `track_mapping`:

  Optional `SeqMap` from the parent track.

#### Returns

The link, invisibly.

------------------------------------------------------------------------

### `SeqLink$prep()`

Override in subclasses. Default implementation errors.

#### Usage

    SeqLink$prep(layout_all_tracks, track_windows_list, plot_track_index = NULL)

#### Arguments

- `layout_all_tracks`:

  Named list of panel bounds keyed by `track_id`.

- `track_windows_list`:

  Named list of `GRanges` windows keyed by `track_id`.

- `plot_track_index`:

  Optional integer index of the parent track, used only to set defaults
  for within-track links.

------------------------------------------------------------------------

### `SeqLink$.resolve_track_ref()`

Look up a track reference (`track_id` string or integer index) in
`layout_all_tracks` and return its panel-bounds list. Errors with a
clear message when the reference cannot be resolved.

#### Usage

    SeqLink$.resolve_track_ref(ref, layout_all_tracks)

#### Arguments

- `ref`:

  A `track_id` (character) or integer index.

- `layout_all_tracks`:

  Named list of per-track panel bounds (the `panelBounds` produced by
  [`SeqPlotR6`](http://andrewlynch.io/SeqPlotR/reference/SeqPlotR6.md)'s
  `layoutGrid()`).

#### Returns

The selected track's panel-bounds list.

------------------------------------------------------------------------

### `SeqLink$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqLink$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
