# SeqElement R6 base class

SeqElement R6 base class

SeqElement R6 base class

## Details

Internal R6 generator for the base class of every drawable element.
Subclasses override `prep()` and `draw()`. Not user-facing — concrete
elements expose snake_case constructors (e.g.
[`seq_point()`](http://andrewlynch.io/SeqPlotR/reference/seq_point.md))
defined in their own files.

## Public fields

- `data`:

  GRanges — element's own data (overrides track data when set).

- `mapping`:

  SeqMap — element's own mapping (merged field-by-field with the parent
  track's mapping).

- `aesthetics`:

  SeqAes — constant aesthetics applied to all glyphs.

- `resolved`:

  Named list populated by `resolve()`. Includes resolved mapping fields
  and `.data` / `.mapping` for downstream use in `prep()`.

- `coordCanvas`:

  List populated by `prep()` with canvas npc coordinates ready for
  `grid` drawing primitives.

- `legend`:

  A single `LegendKey`, a named list of `LegendKey` objects, or `NULL`.
  Defines the legend entry/entries contributed by this element.

- `show_legend`:

  Logical. When `FALSE`, this element contributes no keys to any legend
  regardless of the `legend` field. Default `TRUE`.

- `auto_legend`:

  Auto-generated legend spec (set by `prep()`). A `SeqLegendSpec`,
  `GradientLegendSpec`, list of those, or `NULL`. Only consulted when
  `legend` is `NULL`. Users should not set this directly — use the
  `legend` field instead.

## Methods

### Public methods

- [`SeqElementR6$new()`](#method-SeqElement-new)

- [`SeqElementR6$resolve()`](#method-SeqElement-resolve)

- [`SeqElementR6$prep()`](#method-SeqElement-prep)

- [`SeqElementR6$draw()`](#method-SeqElement-draw)

- [`SeqElementR6$.infer_scale_y()`](#method-SeqElement-.infer_scale_y)

- [`SeqElementR6$collect_legend_keys()`](#method-SeqElement-collect_legend_keys)

- [`SeqElementR6$clone()`](#method-SeqElement-clone)

------------------------------------------------------------------------

### Method `new()`

Construct a new SeqElementR6.

#### Usage

    SeqElementR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      legend = NULL,
      show_legend = TRUE,
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` for the element.

- `mapping`:

  Optional `SeqMap`.

- `aesthetics`:

  Optional `SeqAes`. Defaults to an empty
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md).

- `legend`:

  Optional `LegendKey` or list of `LegendKey` objects.

- `show_legend`:

  Logical. Set to `FALSE` to suppress legend output.

- `...`:

  Unused — accepted so subclasses can pass extra arguments.

------------------------------------------------------------------------

### Method `resolve()`

Resolve the effective data + mapping for this element by merging in the
parent track's defaults. Element fields take priority; missing fields
are inherited from the track.

#### Usage

    SeqElementR6$resolve(track_data = NULL, track_mapping = NULL)

#### Arguments

- `track_data`:

  Optional `GRanges` from the parent track.

- `track_mapping`:

  Optional `SeqMap` from the parent track.

#### Returns

The element, invisibly.

------------------------------------------------------------------------

### Method `prep()`

Override in subclasses. Default implementation errors.

#### Usage

    SeqElementR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  The current track's panel metadata list.

- `track_windows`:

  The current track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Override in subclasses. Default implementation errors.

#### Usage

    SeqElementR6$draw()

------------------------------------------------------------------------

### Method `.infer_scale_y()`

Infer a `SeqPositionScale` for the y-axis from the resolved `y` mapping.
Used by `SeqPlotR6$layoutGrid()` before any `prep()` call so that
[`.compute_track_yscale()`](http://andrewlynch.io/SeqPlotR/reference/dot-compute_track_yscale.md)
has something to work with. Returns `NULL` when no usable y vector is
available.

#### Usage

    SeqElementR6$.infer_scale_y()

------------------------------------------------------------------------

### Method `collect_legend_keys()`

Collect legend keys contributed by this element. Returns a list of
entries, each a named list with fields:

- title:

  Character or `NULL`. The legend group title, taken from the
  `LegendKey`'s `title` field when present.

- key:

  A `LegendKey` object.

- element_class:

  Character. The R6 class name of the contributing element.

Returns `NULL` when `show_legend` is `FALSE` or `legend` is `NULL`.

#### Usage

    SeqElementR6$collect_legend_keys()

#### Returns

A list of entries, or `NULL`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqElementR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
