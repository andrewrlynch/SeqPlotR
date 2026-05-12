# Build the inheritance chain for an axis theme key

Given a dotted key like `axis.x1.line.col`, returns the sequence of
candidate keys to try in order of most to least specific:

## Usage

``` r
.axis_inheritance_chain(key_path)
```

## Arguments

- key_path:

  Dotted key string.

## Value

Character vector of candidate keys.

## Details

axis.x1.line.col → axis.x.line.col → axis.line.col

For non-`axis.*` keys, returns just the key itself (no inheritance).
