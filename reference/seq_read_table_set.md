# Read a table set

Reads a versioned set of related tables, either from a directory or from
an in-memory list of data frames. A directory holds
`<prefix>.<table>.tsv` files, a `<prefix>.columns.tsv` manifest (columns
`table`, `column`, `type`) and a `<prefix>.schema_version` file whose
first line is `<set>`, a tab and an integer version. Column types are
`character`, `integer`, `numeric` or `logical`. Character cells are
percent-decoded (`%09`, `%0A`, `%0D`, `%4E`, then `%25`) and `NA` is
read as a missing value.

## Usage

``` r
seq_read_table_set(x, set, supported = 1L, required = NULL, prefix = set)
```

## Arguments

- x:

  A directory path, or a list of data frames carrying a `schema_version`
  element or attribute. Data tables are accepted.

- set:

  Name of the table set, checked against the version file.

- supported:

  Integer schema versions that are accepted. Any other version is
  rejected.

- required:

  Named list mapping table names to required column names.

- prefix:

  File prefix of the set in a directory.

## Value

A named list of data frames with attribute `schema_version`.
