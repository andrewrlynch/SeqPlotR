#' Read a table set
#'
#' Reads a versioned set of related tables, either from a directory or from an
#' in-memory list of data frames. A directory holds `<prefix>.<table>.tsv`
#' files, a `<prefix>.columns.tsv` manifest (columns `table`, `column`, `type`)
#' and a `<prefix>.schema_version` file whose first line is `<set>`, a tab and
#' an integer version. Column types are `character`, `integer`, `numeric` or
#' `logical`. Character cells are percent-decoded (`%09`, `%0A`, `%0D`, `%4E`,
#' then `%25`) and `NA` is read as a missing value.
#'
#' @param x A directory path, or a list of data frames carrying a
#'   `schema_version` element or attribute. Data tables are accepted.
#' @param set Name of the table set, checked against the version file.
#' @param supported Integer schema versions that are accepted. Any other
#'   version is rejected.
#' @param required Named list mapping table names to required column names.
#' @param prefix File prefix of the set in a directory.
#' @return A named list of data frames with attribute `schema_version`.
#' @export
seq_read_table_set <- function(x, set, supported = 1L, required = NULL, prefix = set) {
  # argument checks ----
  if (!is.character(set) || length(set) != 1L || is.na(set)) {
    stop("seq_read_table_set(): set must be a single string", call. = FALSE)
  }
  fail <- function(...) stop(set, ": ", ..., call. = FALSE)
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) fail("prefix must be a single string")
  supported <- suppressWarnings(as.integer(supported))
  if (length(supported) == 0L || anyNA(supported)) fail("supported must be whole numbers")
  if (!is.null(required) && !(is.list(required) && !is.null(names(required)) &&
        all(vapply(required, is.character, logical(1))))) {
    fail("required must be a named list of column names")
  }
  check_version <- function(ver) {
    if (!(ver %in% supported)) {
      fail("schema version ", ver, " is not supported (supported: ", paste(supported, collapse = ", "), ")")
    }
  }

  if (is.character(x) && length(x) == 1L) {
    # directory input ----
    if (!dir.exists(x)) fail("'", x, "' is not a directory")
    vp <- file.path(x, paste0(prefix, ".schema_version"))
    mp <- file.path(x, paste0(prefix, ".columns.tsv"))
    if (!file.exists(vp) || !file.exists(mp)) {
      fail("no table set with prefix '", prefix, "' in '", x, "' (need ", prefix,
           ".schema_version and ", prefix, ".columns.tsv)")
    }
    fields <- strsplit(readLines(vp, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
    ver <- if (length(fields) == 2L) suppressWarnings(as.numeric(fields[2])) else NA_real_
    if (length(fields) != 2L || is.na(ver) || ver != round(ver)) {
      fail("cannot read the schema version from '", vp, "'")
    }
    ver <- as.integer(ver)
    if (fields[1] != set) fail("'", x, "' holds table set '", fields[1], "', not '", set, "'")
    check_version(ver)

    manifest <- utils::read.delim(mp, colClasses = "character", quote = "", comment.char = "",
                                  na.strings = character(0), stringsAsFactors = FALSE)
    if (!identical(names(manifest), c("table", "column", "type"))) {
      fail("'", mp, "' must have columns table, column, type")
    }
    for (i in seq_len(nrow(manifest))) {
      if (!(manifest$type[i] %in% c("character", "integer", "numeric", "logical"))) {
        fail("unknown column type '", manifest$type[i], "' for ", manifest$table[i], ".",
             manifest$column[i], " in ", prefix, ".columns.tsv")
      }
    }

    out <- list()
    for (tb in unique(manifest$table)) {
      p <- file.path(x, paste0(prefix, ".", tb, ".tsv"))
      if (!file.exists(p)) fail("table file '", prefix, ".", tb, ".tsv' listed in ", prefix, ".columns.tsv is missing")
      d <- utils::read.delim(p, header = TRUE, colClasses = "character", na.strings = "NA",
                             quote = "", comment.char = "", check.names = FALSE,
                             strip.white = FALSE, blank.lines.skip = FALSE, fill = FALSE,
                             stringsAsFactors = FALSE, encoding = "UTF-8")
      m <- manifest[manifest$table == tb, ]
      if (!identical(names(d), m$column)) fail("columns of '", prefix, ".", tb, ".tsv' do not match ", prefix, ".columns.tsv")
      for (j in seq_len(nrow(m))) {
        v <- d[[m$column[j]]]
        d[[m$column[j]]] <- switch(
          m$type[j],
          character = {
            v <- gsub("%09", "\t", v, fixed = TRUE)
            v <- gsub("%0A", "\n", v, fixed = TRUE)
            v <- gsub("%0D", "\r", v, fixed = TRUE)
            v <- gsub("%4E", "N", v, fixed = TRUE)
            gsub("%25", "%", v, fixed = TRUE)
          },
          {
            cv <- suppressWarnings(switch(m$type[j], integer = as.integer(v),
                                          numeric = as.numeric(v), logical = as.logical(v)))
            if (any(is.na(cv) & !is.nan(cv) & !is.na(v))) fail("column ", tb, ".", m$column[j], " holds values that are not ", m$type[j])
            cv
          }
        )
      }
      out[[tb]] <- d
    }
  } else if (is.list(x) && !is.data.frame(x)) {
    # in-memory input ----
    ver <- if (!is.null(x[["schema_version"]])) x[["schema_version"]] else attr(x, "schema_version")
    if (is.null(ver)) fail("x has no schema_version (element or attribute)")
    if (!is.numeric(ver) || length(ver) != 1L || is.na(ver) || ver != round(ver)) {
      fail("schema_version must be a single whole number")
    }
    ver <- as.integer(ver)
    check_version(ver)
    out <- list()
    for (nm in setdiff(names(x), "schema_version")) {
      t <- x[[nm]]
      if (is.null(t)) next
      if (!is.data.frame(t)) fail("element '", nm, "' of x is not a data frame")
      n <- nrow(t)
      d <- lapply(as.list(t), function(v) if (is.factor(v)) as.character(v) else v)
      out[[nm]] <- structure(d, names = names(t), row.names = c(NA_integer_, -n), class = "data.frame")
    }
  } else {
    fail("x must be a directory path or a list of data frames")
  }

  # required tables and columns ----
  for (tb in names(required)) {
    if (is.null(out[[tb]])) fail("required table '", tb, "' is missing")
    miss <- setdiff(required[[tb]], names(out[[tb]]))
    if (length(miss) > 0L) fail("table '", tb, "' lacks required column(s) ", paste(miss, collapse = ", "))
  }
  attr(out, "schema_version") <- ver
  out
}
