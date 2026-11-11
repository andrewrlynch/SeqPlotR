# SeqBlock R6 class

Internal R6 generator backing
[`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md).
Inherits from
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md).

## Super class

`SeqElement` -\> `SeqBlock`

## Public fields

- `fill`:

  `"solid"` or `"cytoband"`.

- `cytobands`:

  Cytoband table (`GRanges` or data frame) or `NULL`.

- `palette`:

  Palette spec passed to the fill colour helper.

- `tints`:

  Tint fractions mixed toward white.

- `block_colors`:

  Named colour overrides by label or key.

- `shape`:

  `"rect"` or `"arrow"`.

- `height`:

  Block height as a fraction of the panel height.

- `label_side`:

  `"below"` or `"above"`.

- `flip_style`:

  `"rotate"`, `"overline"` or `"minus"`.

- `annotations`:

  Optional annotation data frame.

## Methods

### Public methods

- [`SeqBlock$new()`](#method-SeqBlock-initialize)

- [`SeqBlock$prep()`](#method-SeqBlock-prep)

- [`SeqBlock$draw()`](#method-SeqBlock-draw)

- [`SeqBlock$clone()`](#method-SeqBlock-clone)

Inherited methods

- `SeqElement$.infer_scale_y()`
- `SeqElement$collect_legend_keys()`
- `SeqElement$resolve()`

------------------------------------------------------------------------

### `SeqBlock$new()`

Construct a SeqBlockR6.

#### Usage

    SeqBlock$new(
      data = NULL,
      mapping = NULL,
      aesthetics = aes(),
      fill = c("solid", "cytoband"),
      cytobands = NULL,
      palette = NULL,
      tints = c(0.15, 0.5, 0.3, 0.65, 0.4),
      block_colors = NULL,
      shape = c("rect", "arrow"),
      height = 0.5,
      label_side = c("below", "above"),
      flip_style = c("rotate", "overline", "minus"),
      annotations = NULL,
      legend = NULL,
      show_legend = TRUE,
      ...
    )

#### Arguments

- `data`:

  `GRanges` of block positions in track coordinates.

- `mapping`:

  Optional `SeqMap`.

- `aesthetics`:

  Optional `SeqAes`.

- `fill, cytobands, palette, tints, block_colors, shape, height`:

  See
  [`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md).

- `label_side, flip_style, annotations`:

  See
  [`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md).

- `legend`:

  Optional legend.

- `show_legend`:

  Logical.

- `...`:

  Unused.

------------------------------------------------------------------------

### `SeqBlock$prep()`

Resolve mappings, assign blocks to windows and compute npc geometry for
fills, bands, centromeres, labels and annotations.

#### Usage

    SeqBlock$prep(layout_track, track_windows)

#### Arguments

- `layout_track`:

  Per-window panel metadata list.

- `track_windows`:

  Track's `windows` GRanges.

#### Returns

The `coordCanvas` list (one entry per window), invisibly.

------------------------------------------------------------------------

### `SeqBlock$draw()`

Draw fills, bands, centromeres, fade, outlines, labels and annotations
with grid.

#### Usage

    SeqBlock$draw()

------------------------------------------------------------------------

### `SeqBlock$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqBlock$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
