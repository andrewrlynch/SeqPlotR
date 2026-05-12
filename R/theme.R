# ── Hierarchical theme system ────────────────────────────────────────────────
#
# SeqPlotR accepts ggplot2-style hierarchical aesthetic keys inside `aes()`:
#
#   aes(axis.line.col = "grey30",          # all axes
#       axis.x.line.lwd = 1.2,              # both x axes
#       axis.y2.line.col = "steelblue",     # only secondary y axis
#       axis.x1.scale.cap = "capped",       # cap axis line to breaks
#       track.background.fill = "#faf8f5")
#
# Users can also write nested aes() values:
#
#   aes(axis.x = aes(position = "top",
#                    line = aes(col = "red")))
#
# Both forms flatten to the same dotted-key map. At `layoutGrid()` time the
# plot- and track-level `aesthetics` are flattened once and merged (track
# overrides plot). Draw helpers resolve specific keys via `.resolve_theme()`,
# which walks the inheritance chain:
#
#   axis.x1.line.col  →  axis.x.line.col  →  axis.line.col  →  default
#
# The `track.*` namespace is flat — `track.window.*` does NOT inherit from
# `track.*` (semantically different rectangles).

#' Flatten a (possibly nested) SeqAes to a dotted-key map
#'
#' Recursively walks a `SeqAes` or named list, producing a single-level
#' named list where each name is a dot-separated path to a leaf value.
#' Leaves are anything that is not itself a `SeqAes` or a named `list`:
#' scalars, numeric vectors, functions, palettes, `SeqScale*` objects,
#' `NA`, etc. Unnamed elements are dropped.
#'
#' @param x A `SeqAes`, named list, or `NULL`.
#' @param prefix Internal — used during recursion to accumulate the path.
#' @return A named list of leaf values, keyed by dotted paths.
#' @keywords internal
.flatten_theme <- function(x, prefix = "") {
  if (is.null(x) || length(x) == 0L) return(list())
  if (!(inherits(x, "SeqAes") || (is.list(x) && !is.null(names(x))))) {
    # leaf at the root — can't happen unless prefix is non-empty
    if (nzchar(prefix)) {
      out <- list(); out[[prefix]] <- x; return(out)
    }
    return(list())
  }
  out <- list()
  nms <- names(x) %||% rep("", length(x))
  for (i in seq_along(x)) {
    nm <- nms[i]
    if (!nzchar(nm)) next            # drop unnamed entries
    v  <- x[[i]]
    key <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    if (inherits(v, "SeqAes") ||
        (is.list(v) && !is.null(names(v)) && !.is_leaf_list(v))) {
      sub <- .flatten_theme(v, prefix = key)
      for (k in names(sub)) out[[k]] <- sub[[k]]
    } else {
      out[[key]] <- v
    }
  }
  out
}

#' Decide whether a named list should be treated as a theme leaf
#'
#' Some leaf values are themselves named lists (e.g. the `margins` key on
#' the default theme stores `list(top, right, bottom, left)`). Heuristic:
#' if every element of the list is atomic and length 1 and the name set is
#' a subset of known leaf-list names, treat it as a leaf.
#'
#' @param v A named list.
#' @return `TRUE` if `v` should be emitted as a single leaf value.
#' @keywords internal
.is_leaf_list <- function(v) {
  leaf_names <- c("top", "right", "bottom", "left")
  !is.null(names(v)) &&
    length(v) > 0L &&
    all(names(v) %in% leaf_names) &&
    all(vapply(v, function(e) is.atomic(e) && length(e) == 1L, logical(1)))
}

#' Build the inheritance chain for an axis theme key
#'
#' Given a dotted key like `axis.x1.line.col`, returns the sequence of
#' candidate keys to try in order of most to least specific:
#'
#'   axis.x1.line.col → axis.x.line.col → axis.line.col
#'
#' For non-`axis.*` keys, returns just the key itself (no inheritance).
#'
#' @param key_path Dotted key string.
#' @return Character vector of candidate keys.
#' @keywords internal
.axis_inheritance_chain <- function(key_path) {
  parts <- strsplit(key_path, ".", fixed = TRUE)[[1]]
  if (length(parts) < 2L || parts[1] != "axis") return(key_path)
  side <- parts[2]
  # side is one of x, y, x1, x2, y1, y2 — derive the one-letter general form
  gen  <- substr(side, 1L, 1L)
  tail <- if (length(parts) > 2L) parts[-c(1L, 2L)] else character(0)
  candidates <- c(
    paste(c("axis", side, tail), collapse = "."),
    if (!identical(side, gen))
      paste(c("axis", gen,  tail), collapse = "."),
    paste(c("axis",              tail), collapse = ".")
  )
  unique(candidates)
}

#' Resolve a theme key from a flat theme map
#'
#' Walks `.axis_inheritance_chain()` and returns the value stored at the
#' most specific key present. Falls back to `default` if no entry exists.
#'
#' @param flat_theme Named list of dotted keys (from `.flatten_theme()`,
#'   usually merged with `default_theme`).
#' @param key_path Dotted key string (e.g. `"axis.x1.line.col"`).
#' @param default Value returned when no entry is found.
#' @return The resolved value, or `default`.
#' @keywords internal
.resolve_theme <- function(flat_theme, key_path, default = NULL) {
  chain <- .axis_inheritance_chain(key_path)
  for (k in chain) {
    if (k %in% names(flat_theme)) {
      v <- flat_theme[[k]]
      if (!is.null(v)) return(v)
    }
  }
  default
}

#' Merge a track's flat theme over the plot's flat theme
#'
#' Last-write-wins semantics. Track keys override plot keys for identical
#' dotted paths.
#'
#' @param plot_theme Flattened plot-level theme.
#' @param track_theme Flattened track-level theme.
#' @return Merged flat theme.
#' @keywords internal
.merge_themes <- function(plot_theme, track_theme) {
  if (length(track_theme) == 0L) return(plot_theme)
  if (length(plot_theme)  == 0L) return(track_theme)
  out <- plot_theme
  for (k in names(track_theme)) out[[k]] <- track_theme[[k]]
  out
}

#' Translate bare NA values on structural axis keys to `visible = FALSE`
#'
#' Mirrors ggplot's `element_blank()`: when a user writes
#' `aes(axis.x.line = NA)` the intent is to hide that piece. After
#' flattening, scan the theme for the known structural sub-keys
#' (`axis.<side>.line`, `.title`, `.ticks`, `.labels`, `.text`,
#' `.gridline`) whose value is bare `NA`, and rewrite as
#' `<key>.visible = FALSE` (dropping the raw NA leaf so it does not
#' confuse downstream resolvers). Also recognises a bare string on
#' `axis.<side>.title` as a shorthand for `axis.<side>.title.label`.
#'
#' @param flat A flat theme map (from `.flatten_theme()`).
#' @return The transformed flat theme.
#' @keywords internal
.normalize_blanks <- function(flat) {
  if (length(flat) == 0L) return(flat)
  structural <- c("line", "title", "ticks", "labels", "text", "gridline")
  sides      <- c("x", "y", "x1", "x2", "y1", "y2")
  # Pre-compute candidate parent keys.
  parents <- as.vector(outer(
    paste0("axis.", sides),
    paste0(".", structural),
    paste0
  ))
  also_parents <- c("axis.line", "axis.title", "axis.ticks",
                    "axis.labels", "axis.text", "axis.gridline")
  all_parents <- c(parents, also_parents)

  drop <- character(0)
  for (k in intersect(all_parents, names(flat))) {
    v <- flat[[k]]
    if (length(v) == 1L && is.logical(v) && is.na(v)) {
      flat[[paste0(k, ".visible")]] <- FALSE
      drop <- c(drop, k)
      next
    }
    if (length(v) == 1L && is.atomic(v) && is.na(v)) {
      flat[[paste0(k, ".visible")]] <- FALSE
      drop <- c(drop, k)
      next
    }
    # Bare string on a `*.title` parent → label shorthand.
    if (endsWith(k, ".title") && is.character(v) && length(v) == 1L) {
      flat[[paste0(k, ".label")]] <- v
      drop <- c(drop, k)
    }
  }
  if (length(drop) > 0L) flat[drop] <- NULL
  flat
}

# ── Default theme ────────────────────────────────────────────────────────────

#' The built-in default theme for SeqPlotR
#'
#' Flat dotted-key list combining layout parameters (`track_gaps`,
#' `window_gaps`, `margins`) with the hierarchical axis and track chrome
#' keys. Users override via `seq_plot(aesthetics = aes(...))` or
#' `seq_track(aesthetics = aes(...))`.
#'
#' Layout-only keys (`track_gaps`, `window_gaps`, `margins`) are not
#' subject to axis inheritance — they're looked up by exact name.
#'
#' @keywords internal
.default_theme <- function() {
  list(
    # Layout parameters (exact match; no inheritance).
    # Window gap default lives only on the deprecated alias `window_gaps`
    # so an explicit user aes("window_gaps" = ...) can still take effect.
    # The canonical key `window.gap.width` is unset by default; the layout
    # chain resolves to it first, then falls back to `window_gaps`.
    track_gaps  = 0.01,
    window_gaps = 0.01,   # deprecated alias for window.gap.width
    # `margins` is intentionally omitted here so layoutGrid() can detect
    # the absence of user-set margins and fall back to `plot_margin`.

    # Axis defaults — declared at the most general level and inherited by
    # axis.x / axis.y / axis.x1 / axis.x2 / axis.y1 / axis.y2.
    "axis.line.col"            = "#1C1B1A",
    "axis.line.lwd"            = 1,
    "axis.line.alpha"          = 1,
    "axis.line.visible"        = TRUE,
    "axis.ticks.col"           = "#1C1B1A",
    "axis.ticks.lwd"           = 1,
    "axis.ticks.length"        = 0.005,
    "axis.ticks.visible"       = TRUE,
    "axis.text.size"           = 0.6,
    "axis.text.col"            = "#1C1B1A",
    "axis.text.angle"          = 0,
    "axis.text.visible"        = TRUE,
    "axis.title.size"          = 0.8,
    "axis.title.col"           = "#1C1B1A",
    "axis.title.visible"       = TRUE,
    "axis.title.text"          = NULL,
    "axis.scale.cap"           = "capped",
    "axis.scale.expand"        = c(0.025, 0),
    "axis.scale.n_breaks"      = 5,
    "axis.scale.breaks"        = NULL,
    "axis.scale.minor_breaks"  = NULL,
    "axis.scale.limits"        = NULL,
    "axis.scale.oob"           = "exclude",
    "axis.scale.pretty"        = NULL,

    # Side-specific defaults (override the more general axis.* keys).
    "axis.x1.position"         = "bottom",
    "axis.x2.position"         = "top",
    "axis.y1.position"         = "left",
    "axis.y2.position"         = "right",
    "axis.x.text.offset"       = 0.015,
    "axis.y.text.offset"       = 0.010,
    "axis.y.per_window"        = FALSE,

    # In-panel placement and label-style options.
    # `position`: "axis" (default — render in margin band) or c(x_npc, y_npc)
    # for in-panel placement.
    # `style` (labels only): "tick" (default per-tick labels) or "range"
    # (single "[lo–hi]" label).
    "axis.y.title.position"    = "axis",
    "axis.x.title.position"    = "axis",
    "axis.y.labels.style"      = "tick",
    "axis.y.labels.position"   = "axis",
    "axis.x.labels.style"      = "tick",
    "axis.x.labels.position"   = "axis",

    # Gridlines — drawn at tick positions, off by default.
    # Inherit: axis.x.gridline.* → axis.gridline.*  (same hierarchy as all axis.* keys)
    # Enable:  axis.x.gridline = TRUE  OR  axis.x.gridline = aes(color, lwd, lty, alpha)
    "axis.gridline.visible" = FALSE,
    "axis.gridline.color"   = "grey85",
    "axis.gridline.lwd"     = 0.5,
    "axis.gridline.lty"     = 1,
    "axis.gridline.alpha"   = 1,

    # Track and per-window chrome.
    "track.background.fill"           = NA,
    "track.background.alpha"          = 1,
    "track.border.col"                = NA,
    "track.border.lwd"                = 0.5,
    "track.border.alpha"              = 1,
    "track.window.background.fill"    = NA,
    "track.window.background.alpha"   = 1,
    "track.window.border.col"         = NA,
    "track.window.border.lwd"         = 1,
    "track.window.border.alpha"       = 1
  )
}

# ── Resolved per-axis spec ───────────────────────────────────────────────────

#' Precompute a per-side axis spec from a flat theme
#'
#' For one of `c("x1","x2","y1","y2")`, walk the hierarchy once per leaf
#' key and package the resolved values into a nested list the axis-draw
#' helpers consume directly. This avoids repeated chain walks at draw
#' time.
#'
#' @param flat_theme Flattened theme (already merged plot+track).
#' @param side One of `"x1"`, `"x2"`, `"y1"`, `"y2"`.
#' @return A list with elements `position`, `line`, `ticks`, `text`,
#'   `title`, `scale`, plus convenience scalars `axis_dim` and
#'   `axis_index`.
#' @keywords internal
.build_axis_spec <- function(flat_theme, side) {
  .get <- function(leaf, default = NULL)
    .resolve_theme(flat_theme, paste0("axis.", side, ".", leaf), default)

  list(
    side       = side,
    axis_dim   = substr(side, 1L, 1L),            # "x" or "y"
    axis_index = as.integer(substr(side, 2L, 2L)),# 1 or 2
    position   = .get("position",
                      if (side == "x1") "bottom"
                      else if (side == "x2") "top"
                      else if (side == "y1") "left"
                      else "right"),
    visible    = .get("visible", NA),             # NA = auto
    line = list(
      col     = .get("line.col", "#1C1B1A"),
      lwd     = .get("line.lwd", 1),
      alpha   = .get("line.alpha", 1),
      visible = .get("line.visible", TRUE)
    ),
    ticks = list(
      col     = .get("ticks.col", .get("line.col", "#1C1B1A")),
      lwd     = .get("ticks.lwd", .get("line.lwd", 1)),
      length  = .get("ticks.length", 0.005),
      alpha   = .get("ticks.alpha", 1),
      visible = .get("ticks.visible", TRUE)
    ),
    text = list(
      size     = .get("text.size", 0.6),
      col      = .get("text.col",  .get("line.col", "#1C1B1A")),
      angle    = .get("text.angle", 0),
      offset   = .get("text.offset", if (substr(side,1,1) == "x") 0.015 else 0.010),
      visible  = .get("text.visible", TRUE),
      style    = .get("labels.style",    "tick"),
      position = .get("labels.position", "axis"),
      hjust    = .get("labels.hjust",    NULL),
      vjust    = .get("labels.vjust",    NULL)
    ),
    title = list(
      size     = .get("title.size", 0.8),
      col      = .get("title.col",  .get("line.col", "#1C1B1A")),
      text     = .get("title.text", NULL),
      label    = .get("title.label", NULL),
      visible  = .get("title.visible", TRUE),
      position = .get("title.position", "axis"),
      hjust    = .get("title.hjust",    NULL),
      vjust    = .get("title.vjust",    NULL)
    ),
    scale = list(
      cap           = .get("scale.cap", "capped"),
      expand        = .get("scale.expand", c(0.025, 0)),
      n_breaks      = .get("scale.n_breaks", 5),
      breaks        = .get("scale.breaks", NULL),
      minor_breaks  = .get("scale.minor_breaks", NULL),
      limits        = .get("scale.limits", NULL),
      labels        = .get("scale.labels", NULL),
      oob           = .get("scale.oob", "exclude"),
      pretty        = .get("scale.pretty", NULL)
    )
  )
}

#' Precompute the resolved track chrome spec
#'
#' @param flat_theme Flattened theme.
#' @return Nested list with `background`, `border`, and `window` sections.
#' @keywords internal
.build_chrome_spec <- function(flat_theme) {
  .g <- function(k, d = NULL)
    if (k %in% names(flat_theme)) flat_theme[[k]] else d
  list(
    background = list(
      fill  = .g("track.background.fill",  NA),
      alpha = .g("track.background.alpha", 1)
    ),
    border = list(
      col   = .g("track.border.col",  NA),
      lwd   = .g("track.border.lwd",  0.5),
      alpha = .g("track.border.alpha", 1)
    ),
    window = list(
      background = list(
        fill  = .g("track.window.background.fill",  "whitesmoke"),
        alpha = .g("track.window.background.alpha", 1)
      ),
      border = list(
        col   = .g("track.window.border.col",  "grey50"),
        lwd   = .g("track.window.border.lwd",  1),
        alpha = .g("track.window.border.alpha", 1)
      )
    )
  )
}

#' Build a complete resolved theme for one track
#'
#' Combines the per-axis specs for all four sides, the chrome spec, and
#' a few convenience scalars. Stored on the track for the draw pipeline.
#'
#' @param flat_theme Merged flat theme (plot + track).
#' @return Nested list.
#' @keywords internal
.build_resolved_theme <- function(flat_theme) {
  list(
    flat   = flat_theme,
    axes   = list(
      x1 = .build_axis_spec(flat_theme, "x1"),
      x2 = .build_axis_spec(flat_theme, "x2"),
      y1 = .build_axis_spec(flat_theme, "y1"),
      y2 = .build_axis_spec(flat_theme, "y2")
    ),
    chrome = .build_chrome_spec(flat_theme),
    y_per_window =
      isTRUE(.resolve_theme(flat_theme, "axis.y.per_window", FALSE))
  )
}
