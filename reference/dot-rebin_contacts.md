# Rebin a contact data.frame to a coarser bin size

Aggregates an existing BEDPE-style contact data.frame (as returned by
`$fetch()`) into larger bins using `agg_fun`.

## Usage

``` r
.rebin_contacts(contacts, bin_size, agg_fun)
```

## Arguments

- contacts:

  data.frame with columns seqnames1, start1, end1, seqnames2, start2,
  end2, score.

- bin_size:

  Integer. New (coarser) bin size in bp.

- agg_fun:

  Resolved aggregation function.

## Value

data.frame with the same columns, at the new resolution.
