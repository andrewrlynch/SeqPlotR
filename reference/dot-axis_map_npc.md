# Map a data value to a canvas npc coordinate inside a panel range

Map a data value to a canvas npc coordinate inside a panel range

## Usage

``` r
.axis_map_npc(val, scale_lim, npc_lo, npc_hi, flip = FALSE)
```

## Arguments

- val:

  Numeric values to map.

- scale_lim:

  Length-2 data-range (from which `val` comes).

- npc_lo, npc_hi:

  Canvas npc endpoints of the panel span.

- flip:

  Logical. When `TRUE`, mirror the npc output around the midpoint of
  `[npc_lo, npc_hi]` — i.e. low data values render at `npc_hi` and high
  values at `npc_lo`.

## Value

Numeric vector of canvas npc coordinates.
