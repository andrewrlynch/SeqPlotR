# Infer C vs. S curve type from a pair of strands

Opposing strands (`+/-` or `-/+`) produce an S curve; matching or
unknown strands fall back to the supplied default (typically `"c"`).

## Usage

``` r
.string_type_from_strand(strand0, strand1, default = "c")
```

## Arguments

- strand0, strand1:

  Character vectors of strands (`"+"`, `"-"`, `"*"`).

- default:

  Fallback curve type when strand info is missing or matching.

## Value

Character vector of `"c"` or `"s"` values.
