# Convert polar coordinates to npc Cartesian

Mathematical convention: 0 degrees = 3 o'clock, counter-clockwise
positive.

## Usage

``` r
.polar_to_npc(r, theta_deg, cx = 0.5, cy = 0.5)
```

## Arguments

- r:

  Numeric npc radius.

- theta_deg:

  Numeric angle in degrees.

- cx, cy:

  Numeric npc coordinates of the circle centre.

## Value

`list(x, y)` in npc.
