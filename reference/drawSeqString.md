# Draw a single cubic-Bezier string between two anchors

Internal helper used by `seq_string`. Draws a smooth horizontal Bezier
curve connecting `(x0, y0)` and `(x1, y1)` in canvas npc coordinates.
The shape (C vs. S) is determined by the strand pair — matching strands
bend both control points the same direction, opposing strands bend them
in opposite directions.

## Usage

``` r
drawSeqString(
  x0,
  y0,
  x1,
  y1,
  strand1 = "*",
  strand2 = "*",
  orientation = "*",
  type = c("c", "s"),
  bulge = 0.04,
  lwd = 1.5,
  col = "red",
  alpha = 1
)
```

## Arguments

- x0, y0:

  Canvas npc coordinates of the first anchor.

- x1, y1:

  Canvas npc coordinates of the second anchor.

- strand1, strand2:

  Strand character for each anchor (`"+"`, `"-"`, `"*"`).

- orientation:

  Optional explicit direction (`"+"` / `"-"`); when unset, direction is
  inferred from `strand1`.

- type:

  `"c"` or `"s"`. Not consulted directly — retained for API
  compatibility; the C vs. S shape is driven by the strand pair.

- bulge:

  Horizontal control-point offset in npc units. Clamped to `[0, 0.35]`.

- lwd, col, alpha:

  Line aesthetics.

## Value

Invisible `NULL`.
