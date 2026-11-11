# SeqTrackLabel R6 class

Internal R6 generator backing
[`seq_track_label()`](http://andrewlynch.io/SeqPlotR/reference/seq_track_label.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqElement` -\> `SeqTrackLabel`

## Public fields

- `label`:

  Row label text.

- `sublabel`:

  Optional second-line text.

- `side`:

  `"left"` or `"right"`.

## Methods

### Public methods

- [`SeqTrackLabel$new()`](#method-SeqTrackLabel-initialize)

- [`SeqTrackLabel$prep()`](#method-SeqTrackLabel-prep)

- [`SeqTrackLabel$draw()`](#method-SeqTrackLabel-draw)

- [`SeqTrackLabel$clone()`](#method-SeqTrackLabel-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$resolve()`

------------------------------------------------------------------------

### `SeqTrackLabel$new()`

Construct a SeqTrackLabelR6.

#### Usage

    SeqTrackLabel$new(
      label,
      sublabel = NULL,
      side = c("left", "right"),
      aesthetics = aes(),
      ...
    )

#### Arguments

- `label`:

  Label text.

- `sublabel`:

  Optional second line.

- `side`:

  `"left"` or `"right"`.

- `aesthetics`:

  Optional `SeqAes`.

- `...`:

  Unused.

------------------------------------------------------------------------

### `SeqTrackLabel$prep()`

Place the label beside the track's window area.

#### Usage

    SeqTrackLabel$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges (unused).

#### Returns

The `coordCanvas` list, invisibly.

------------------------------------------------------------------------

### `SeqTrackLabel$draw()`

Draw the label and optional sublabel.

#### Usage

    SeqTrackLabel$draw()

------------------------------------------------------------------------

### `SeqTrackLabel$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqTrackLabel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
