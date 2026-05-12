# SeqSequence R6 class

SeqSequence R6 class

SeqSequence R6 class

## Details

Internal R6 generator backing
[`seq_sequence()`](http://andrewlynch.io/SeqPlotR/reference/seq_sequence.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).
Renders per-base coloured rectangles (and optionally letters) for
genomic windows up to 200 bp wide. Wider windows emit a message and
render nothing.

## Super class

`SeqPlotR::SeqElement` -\> `SeqSequence`

## Public fields

- `genome`:

  Character. BSgenome package name (e.g.
  `"BSgenome.Hsapiens.UCSC.hg38"`). Required when `sequence` is `NULL`.

- `sequence`:

  Character string of nucleotides. When provided, `genome` is ignored.
  The string is assumed to span the first window.

- `show_letters`:

  Logical. When `TRUE` and window width \<= 80 bp, draw the nucleotide
  letter centred on each rectangle. Default `FALSE`.

- `rect_height`:

  Numeric in (0, 1\] or `NULL`. Fraction of track height allocated to
  each nucleotide rectangle. `NULL` (default) uses the full track
  height.

- `colors`:

  Named character vector mapping nucleotide characters to hex color
  strings. Defaults to UCSC standard.

## Methods

### Public methods

- [`SeqSequenceR6$new()`](#method-SeqSequence-new)

- [`SeqSequenceR6$prep()`](#method-SeqSequence-prep)

- [`SeqSequenceR6$draw()`](#method-SeqSequence-draw)

- [`SeqSequenceR6$clone()`](#method-SeqSequence-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqElement$resolve()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-resolve)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqSequenceR6.

#### Usage

    SeqSequenceR6$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      genome = NULL,
      sequence = NULL,
      show_letters = FALSE,
      rect_height = NULL,
      colors = NULL,
      ...
    )

#### Arguments

- `data`:

  Ignored; sequence is fetched from `genome` or `sequence`.

- `mapping`:

  Ignored.

- `aesthetics`:

  Optional `SeqAes`. Supports `color` (letter color, defaults to
  matching the rectangle fill), `background` (letter background
  rectangle color, default `NA`).

- `genome`:

  Character. BSgenome package name. Required when `sequence` is `NULL`.

- `sequence`:

  Character string of nucleotides spanning the first window. When
  provided, `genome` is ignored.

- `show_letters`:

  Logical. Show nucleotide letters when window width is \<= 80 bp.
  Default `FALSE`.

- `rect_height`:

  Numeric in (0, 1\] or `NULL`. Fraction of track height used for the
  rectangles. Default `NULL` (full track height).

- `colors`:

  Named character vector mapping nucleotide codes to hex colors.
  Defaults to UCSC standard.

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `prep()`

Fetch or validate the nucleotide sequence, check window length, and
compute per-base canvas coordinates.

#### Usage

    SeqSequenceR6$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  The current track's `windows` GRanges.

------------------------------------------------------------------------

### Method `draw()`

Draw per-base coloured rectangles and optional letters.

#### Usage

    SeqSequenceR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqSequenceR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
