# ── map() ─────────────────────────────────────────────────────────────────────

#' Capture unevaluated mapping expressions
#'
#' Captures R expressions to be evaluated later against GRanges data. Each
#' expression may be a bare column name (resolved against `mcols(data)`), one of
#' the genomic specials (`start`, `end`, `width`, `mid`), or an arbitrary R
#' expression (e.g. `log2(score + 1)`, `(start + end) / 2`).
#'
#' Evaluation is deferred until `prep()` time inside each element, where
#' `.resolve_mapping()` is called.
#'
#' @param ... Named expressions (e.g. `x = start`, `y = log2(score + 1)`).
#' @return A `SeqMap` object: a named list of unevaluated language objects.
#' @examples
#' m <- map(x = start, y = score)
#' @export
map <- function(...) {
  structure(as.list(substitute(list(...)))[-1], class = "SeqMap")
}

# ── aes() ─────────────────────────────────────────────────────────────────────

#' Capture evaluated constant aesthetics
#'
#' Captures static aesthetic values (colors, line widths, sizes) that do not
#' vary per observation. For per-observation, data-driven aesthetics, use
#' [map()] instead.
#'
#' @param ... Named values (e.g. `color = "blue"`, `linewidth = 1.5`).
#' @return A `SeqAes` object: a named list of evaluated values.
#' @examples
#' a <- aes(color = "blue", linewidth = 1.5)
#' @export
aes <- function(...) {
  structure(list(...), class = "SeqAes")
}

# ── .resolve_mapping() ────────────────────────────────────────────────────────

#' Evaluate a SeqMap against GRanges or data.frame data
#'
#' Internal helper. Builds an evaluation environment from `data` and
#' evaluates each expression in `mapping` against it.
#'
#' For a `GRanges`, the env is the union of:
#'   * positional specials: `start`, `end`, `width`, `mid`
#'   * GRanges accessors:   `seqnames`, `strand` (coerced to character)
#'   * the columns of `mcols(data)`
#'
#' For a `data.frame`, the env is just `as.list(data)` — no specials are
#' injected, so column names must match the user's `map()` references
#' exactly. (Link elements rely on this: BEDPE-like inputs reference their
#' own column names via `map(x0 = start1, x1 = start2, ...)`.)
#'
#' If `mapping` is `NULL` or `data` is `NULL`, returns `list()`.
#'
#' @param data A `GRanges` or `data.frame`.
#' @param mapping A `SeqMap` object (or `NULL`).
#' @param env The enclosing environment for expression evaluation.
#' @return A named list of resolved vectors, one per mapping field.
#' @keywords internal
.resolve_mapping <- function(data, mapping, env = parent.frame()) {
  if (is.null(mapping) || is.null(data)) return(list())

  if (inherits(data, "GRanges")) {
    n <- length(data)
    specials <- list(
      start    = BiocGenerics::start(data),
      end      = BiocGenerics::end(data),
      width    = BiocGenerics::width(data),
      seqnames = as.character(GenomicRanges::seqnames(data)),
      strand   = as.character(BiocGenerics::strand(data))
    )
    specials$mid <- (specials$start + specials$end) / 2
    eval_env <- c(specials, as.list(S4Vectors::mcols(data)))
  } else if (is.data.frame(data)) {
    n <- nrow(data)
    eval_env <- as.list(data)
  } else {
    stop("data must be a GRanges or a data.frame; got '",
         class(data)[1], "'.", call. = FALSE)
  }

  lapply(mapping, function(expr) {
    v <- eval(expr, envir = eval_env, enclos = env)
    # Broadcast scalar literals (map(x = 0), map(color = "red"), etc.) to
    # the data length so downstream [idx] indexing works. Only broadcast
    # atomic scalars — a closure or environment that bubbles up from the
    # caller is passed through unchanged (downstream code ignores it).
    if (length(v) == 1L && n > 1L && (is.atomic(v) || is.factor(v)))
      rep(v, n)
    else
      v
  })
}

# ── .aes_to_gpar() ────────────────────────────────────────────────────────────

#' Translate a SeqAes object to grid::gpar()
#'
#' Internal helper. Maps SeqAes field names onto grid graphical parameter
#' names. `NULL` values are silently dropped by `grid::gpar()` — this is
#' intentional, so absent fields fall through to grid defaults.
#'
#' @param a A `SeqAes` object.
#' @return A `gpar` object.
#' @keywords internal
.aes_to_gpar <- function(a) {
  grid::gpar(
    col      = a$color %||% a$col,
    fill     = a$fill,
    lwd      = a$linewidth %||% a$lwd,
    lty      = a$linetype,
    alpha    = a$alpha,
    cex      = a$size,
    fontsize = a$fontsize
  )
}
