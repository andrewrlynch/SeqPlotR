# SeqHighlight R6 class

SeqHighlight R6 class

SeqHighlight R6 class

## Details

Internal R6 generator backing
[`seq_highlight()`](http://andrewlynch.io/SeqPlotR/reference/seq_highlight.md).
Inherits from
[`SeqLinkR6`](http://andrewlynch.io/SeqPlotR/reference/SeqLinkR6.md).
Draws a continuous filled highlight band that passes through every track
from `t0` down to `t1` (inclusive). Per-track widths in NPC follow each
track's own genomic scale, so the band fans / compresses across scale
changes. Adjacent tracks are bridged by trapezoids whose top / bottom
edges match the corresponding track's projected `xL` / `xR`.

## Super classes

`SeqPlotR::SeqElement` -\> `SeqPlotR::SeqLink` -\> `SeqHighlight`

## Methods

### Public methods

- [`SeqHighlightR6$new()`](#method-SeqHighlight-new)

- [`SeqHighlightR6$resolve()`](#method-SeqHighlight-resolve)

- [`SeqHighlightR6$prep()`](#method-SeqHighlight-prep)

- [`SeqHighlightR6$draw()`](#method-SeqHighlight-draw)

- [`SeqHighlightR6$clone()`](#method-SeqHighlight-clone)

Inherited methods

- [`SeqPlotR::SeqElement$.infer_scale_y()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-.infer_scale_y)
- [`SeqPlotR::SeqElement$collect_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/SeqElement.html#method-collect_legend_keys)
- [`SeqPlotR::SeqLink$.resolve_track_ref()`](http://andrewlynch.io/SeqPlotR/reference/SeqLink.html#method-.resolve_track_ref)

------------------------------------------------------------------------

### Method `new()`

Construct a SeqHighlightR6.

#### Usage

    SeqHighlightR6$new(
      data = NULL,
      mapping = NULL,
      t0 = NULL,
      t1 = NULL,
      aesthetics = aes(),
      ...
    )

#### Arguments

- `data`:

  Optional `GRanges` or `data.frame` carrying highlight region(s).

- `mapping`:

  Optional `SeqMap`. Required: `x0`, `x0_end`. Optional: `chrom0`
  (auto-derived from GRanges seqnames). The same genomic region is
  projected onto every track in `[t0..t1]`.

- `t0, t1`:

  Track identifiers (string or integer index) bracketing the inclusive
  run of tracks the band passes through. When `t1` is `NULL`, the
  highlight is restricted to `t0` only.

- `aesthetics`:

  Optional `SeqAes`. Recognised: `fill` (default `"grey50"`), `alpha`
  (default `0.25`), `color` (border, default `NA`), `linewidth` (default
  `0.5`).

- `...`:

  Reserved.

------------------------------------------------------------------------

### Method `resolve()`

Resolve the mapping. Like
[`SeqZoomR6`](http://andrewlynch.io/SeqPlotR/reference/SeqZoomR6.md),
this only requires `x0` / `x0_end` — there is no `x1` / `x1_end`, since
the same genomic region is reused on every track in the band.

#### Usage

    SeqHighlightR6$resolve(track_data = NULL, track_mapping = NULL)

#### Arguments

- `track_data`:

  Optional `GRanges`/`data.frame` from the parent track (used as a
  fallback when the link has no own `data`).

- `track_mapping`:

  Optional `SeqMap` from the parent track.

#### Returns

The link, invisibly.

------------------------------------------------------------------------

### Method `prep()`

Project each highlight region onto every track in `[t0..t1]` and stash a
per-rectangle `data.frame` in `self$coordCanvas` with columns
`region_id`, `track_pos`, `window_idx`, `xL`, `xR`, `y0_npc`, `y1_npc`.

#### Usage

    SeqHighlightR6$prep(
      layout_all_tracks,
      track_windows_list,
      plot_track_index = NULL
    )

#### Arguments

- `layout_all_tracks`:

  Named list of per-track panel-bounds lists.

- `track_windows_list`:

  Named list of per-track `GRanges` windows.

- `plot_track_index`:

  Unused — `seq_highlight` requires explicit `t0` (and optionally `t1`).

------------------------------------------------------------------------

### Method `draw()`

Render the highlight band: one filled rectangle per per-track /
per-window slice, and a connecting trapezoid in the gap between adjacent
tracks (matched by `window_idx`).

#### Usage

    SeqHighlightR6$draw()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqHighlightR6$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
