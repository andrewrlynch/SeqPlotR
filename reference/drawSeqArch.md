# Draw a single cubic-Bezier arch and its vertical stems

Internal helper used by `seq_arc`, `seq_arch`, and `seq_recon`. Draws
the Bezier arch connecting two stem tops at canvas npc coordinates, then
the two vertical stems. Pass `stemWidth = 0` to suppress stem rendering.

## Usage

``` r
drawSeqArch(
  x0,
  y0,
  x1,
  y1,
  top0,
  top1,
  orientation = "*",
  curve = "length",
  stemWidth = 1,
  arcWidth = 1,
  arcColor = "black",
  stemColor = "black"
)
```

## Arguments

- x0, y0:

  Canvas npc coordinates of the first stem base.

- x1, y1:

  Canvas npc coordinates of the second stem base.

- top0, top1:

  Canvas npc y-coordinates where the arch meets the stem tops at `x0`
  and `x1`.

- orientation:

  Character. Controls which side the arch bulges to. `"+"` / `"*"`
  curves up; `"-"` curves down. Mixed strings such as `"+/-"` control
  each end independently.

- curve:

  `"length"` (offset scales with span), `"equal"` (fixed `0.2` offset),
  or a numeric value used directly as the offset.

- stemWidth, stemColor:

  Stem aesthetic.

- arcWidth, arcColor:

  Arch aesthetic.

## Value

Invisible `NULL`.
