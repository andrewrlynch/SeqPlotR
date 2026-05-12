# Build a filled annular-sector polygon as xy vertex vectors

Samples `n_pts` points along the outer arc (radius `r1`) from
`theta0_deg` to `theta1_deg`, then `n_pts` points along the inner arc
(radius `r0`) reversed to close the polygon.

## Usage

``` r
.arc_polygon(theta0_deg, theta1_deg, r0, r1, cx = 0.5, cy = 0.5, n_pts = 60)
```

## Arguments

- theta0_deg, theta1_deg:

  Sector angular bounds (degrees).

- r0, r1:

  Inner and outer radii (npc).

- cx, cy:

  Circle centre (npc).

- n_pts:

  Samples per arc edge. Default `60`.

## Value

`list(x, y)` of length `2 * n_pts` each.
