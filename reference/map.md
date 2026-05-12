# Capture unevaluated mapping expressions

Captures R expressions to be evaluated later against GRanges data. Each
expression may be a bare column name (resolved against `mcols(data)`),
one of the genomic specials (`start`, `end`, `width`, `mid`), or an
arbitrary R expression (e.g. `log2(score + 1)`, `(start + end) / 2`).

## Usage

``` r
map(...)
```

## Arguments

- ...:

  Named expressions (e.g. `x = start`, `y = log2(score + 1)`).

## Value

A `SeqMap` object: a named list of unevaluated language objects.

## Details

Evaluation is deferred until `prep()` time inside each element, where
[`.resolve_mapping()`](http://andrewlynch.io/SeqPlotR/reference/dot-resolve_mapping.md)
is called.

## Examples

``` r
m <- map(x = start, y = score)
```
