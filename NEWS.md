# SeqPlotR 0.1.0.9000

## Gene models (`seq_gene()`)

* **Breaking:** the new `oob` argument defaults to `"truncate"`, so genes that
  straddle a window edge are now clipped to the window. Pass `oob = "keep"` to
  restore the previous behaviour of drawing out-of-window portions;
  `oob = "remove"` drops any gene not entirely within the window.
* `label_position = "flush"` packing now accounts for label width: each gene's
  packing interval is widened on the label side by its measured rendered label
  width plus `label_pad` (previously unused), so genes whose labels would
  collide are pushed onto separate tiers. `label_pad` defaults to `50000`.
* `expanded = TRUE` now ranks genes **per window** (a gene absent from a window
  no longer reserves an empty row there) while keeping glyph height constant
  across windows.
* New styling arguments: `backbone_lwd`, `exon_lwd`, `exon_border_col`,
  `arrow_size`, and TSS-marker controls `tss_arrow_size`, `tss_stem_lwd`,
  `tss_rise`, `tss_run`.
* The `show_start` TSS marker is now drawn as a solid filled triangle rather
  than a grid arrow (which rendered as an unfilled outline on many devices),
  and it leads slightly ahead of the arm.

## Lollipops (`seq_lollipop()`)

* Fixed: negative values now place the head at the value rather than at the
  baseline.
* `map(color = ...)` (and `fill`/`shape`) are now honoured per point and
  generate a legend, matching `seq_point()`.

## Highlights (`seq_highlight()`)

* Each highlight region is now drawn as a single closed polygon per disjoint
  piece, so a border colour strokes only the outer boundary instead of every
  internal rectangle and trapezoid seam.

## Layout

* New canonical theme key `aes("track.gap.width" = <npc>)` controls the gap
  between tracks, mirroring `window.gap.width`. The `track_gaps` key remains a
  supported alias.
