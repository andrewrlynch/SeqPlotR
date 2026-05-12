# Evaluate a SeqMap against GRanges or data.frame data

Internal helper. Builds an evaluation environment from `data` and
evaluates each expression in `mapping` against it.

## Usage

``` r
.resolve_mapping(data, mapping, env = parent.frame())
```

## Arguments

- data:

  A `GRanges` or `data.frame`.

- mapping:

  A `SeqMap` object (or `NULL`).

- env:

  The enclosing environment for expression evaluation.

## Value

A named list of resolved vectors, one per mapping field.

## Details

For a `GRanges`, the env is the union of:

- positional specials: `start`, `end`, `width`, `mid`

- GRanges accessors: `seqnames`, `strand` (coerced to character)

- the columns of `mcols(data)`

For a `data.frame`, the env is just `as.list(data)` — no specials are
injected, so column names must match the user's
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) references
exactly. (Link elements rely on this: BEDPE-like inputs reference their
own column names via `map(x0 = start1, x1 = start2, ...)`.)

If `mapping` is `NULL` or `data` is `NULL`, returns
[`list()`](https://rdrr.io/r/base/list.html).
