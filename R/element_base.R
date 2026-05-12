# ── SeqElementR6 ─────────────────────────────────────────────────────────────
#
# Internal R6 base class for every drawable element (points, lines, polygons,
# text, etc.). Concrete element subclasses live in their own files and follow
# the init / prep / draw three-stage contract described in the project spec.

# ── Shared aesthetic-scaling helpers ──────────────────────────────────────────

# Return TRUE if every value in v is already a valid R color.
.looks_like_color <- function(v) {
  if (is.null(v) || length(v) == 0L) return(FALSE)
  v <- as.character(v)
  ok <- grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", v) |
        v %in% grDevices::colors() |
        is.na(v)
  all(ok)
}


# Apply a color scale to a vector of mapped values; return concrete colors and
# an auto-legend spec.
#
# Arguments:
#   vals        Vector from resolved mapping (character, factor, or numeric).
#   col_name    Source column name string — used as legend title.
#   aes_name    "color" or "fill" — which aesthetic to build the legend for.
#   scale       Optional SeqScale object overriding the default palette/limits.
#   position    Legend position. Default "inside".
#   x, y        Anchor within the target area. Default 0.02 / 0.95.
#   hjust       Horizontal justification. Default 0.
#   orientation Layout direction. Default "vertical".
#   side        For margin positions only.
#
# Returns list(colors = chr_vector, legend = SeqLegendSpec or GradientLegendSpec or NULL)
.auto_scale_colors <- function(vals, col_name = "value", aes_name = "color",
                                scale = NULL,
                                position = "track_margin", x = 0.05, y = 0.8,
                                hjust = 0, orientation = "vertical",
                                side = "right") {
  if (is.null(vals)) return(list(colors = NULL, legend = NULL))

  # Already concrete colors — pass through with no legend
  if (.looks_like_color(vals))
    return(list(colors = as.character(vals), legend = NULL))

  is_discrete <- is.factor(vals) || is.character(vals)

  if (is_discrete) {
    vals_f <- as.factor(vals)
    lvls   <- levels(vals_f)
    pal    <- if (!is.null(scale) && !is.null(scale$values)) {
      scale$values
    } else if (!is.null(scale) && !is.null(scale$palette) &&
               is.function(scale$palette)) {
      scale$palette(length(lvls))
    } else {
      flexoki_palette(length(lvls))
    }
    pal_named <- setNames(pal[seq_along(lvls)], lvls)
    colors    <- unname(pal_named[as.character(vals_f)])
    colors[is.na(colors)] <- "grey80"

    keys <- lapply(lvls, function(lv) {
      LegendKey(label = lv, color = pal_named[[lv]], fill = pal_named[[lv]])
    })
    leg <- seq_legend(keys = keys, title = col_name,
                      position = position, x = x, y = y,
                      hjust = hjust, orientation = orientation, side = side)
    return(list(colors = colors, legend = leg))
  }

  # Numeric — continuous gradient
  vals_n  <- as.numeric(vals)
  pal_nm  <- if (!is.null(scale)) scale$palette %||% "viridis" else "viridis"
  lims    <- if (!is.null(scale) && !is.null(scale$limits)) scale$limits
             else range(vals_n, na.rm = TRUE)
  if (diff(lims) == 0) lims <- lims + c(-1, 1)

  t_vals <- pmax(0, pmin(1, (vals_n - lims[1]) / diff(lims)))
  colors <- .palette_color_at(t_vals, pal_nm)
  colors[is.na(vals_n)] <- "grey80"

  leg <- seq_gradient_legend(palette = pal_nm, limits = lims, title = col_name,
                              position = position, x = x, y = y,
                              hjust = hjust, orientation = orientation, side = side)
  list(colors = colors, legend = leg)
}


# Apply a discrete shape scale to a mapped vector.
# Only character/factor values are supported; numeric falls back to NULL.
#
# Returns list(shapes = chr_vector, legend = SeqLegendSpec or NULL)
.auto_scale_shapes <- function(vals, col_name = "value",
                                position = "track_margin", x = 0.05, y = 0.7,
                                hjust = 0, orientation = "vertical",
                                side = "right") {
  if (is.null(vals)) return(list(shapes = NULL, legend = NULL))
  if (!is.factor(vals) && !is.character(vals))
    return(list(shapes = NULL, legend = NULL))

  shape_cycle <- c("circle", "square", "triangle", "diamond")
  vals_f    <- as.factor(vals)
  lvls      <- levels(vals_f)
  shape_map <- setNames(
    shape_cycle[((seq_along(lvls) - 1L) %% length(shape_cycle)) + 1L],
    lvls
  )
  shapes <- unname(shape_map[as.character(vals_f)])

  keys <- lapply(lvls, function(lv) {
    LegendKey(label = lv, shape = shape_map[[lv]], color = "#1C1B1A")
  })
  leg <- seq_legend(keys = keys, title = col_name,
                    position = position, x = x, y = y,
                    hjust = hjust, orientation = orientation, side = side)
  list(shapes = shapes, legend = leg)
}


# Extract the source column name from a SeqMap entry (language object).
# Returns a plain string, e.g. deparse(quote(type)) = "type".
.map_col_name <- function(mapping, field) {
  expr <- mapping[[field]]
  if (is.null(expr)) return(field)
  deparse(expr)
}


#' Normalise `map(axis.x = ...)` / `map(axis.y = ...)` to a scalar integer
#'
#' Axis selectors are per-element — they route the element to the primary
#' (1) or secondary (2) axis. `.resolve_mapping()` broadcasts scalars to
#' the data length; collapse back to a single integer and validate.
#'
#' @param v Raw resolved value (may be `NULL`, vector, or scalar).
#' @param which Either `"x"` or `"y"` — used for the error message.
#' @return Integer scalar in `{1L, 2L}`. Defaults to `1L` when `v` is NULL.
#' @keywords internal
.normalize_axis_selector <- function(v, which = "x") {
  if (is.null(v)) return(1L)
  uv <- unique(as.integer(v))
  if (length(uv) != 1L || !uv %in% c(1L, 2L))
    stop(sprintf("map(axis.%s = ) must be a single 1 or 2 (got %s).",
                 which, paste(unique(v), collapse = ", ")),
         call. = FALSE)
  uv
}

#' SeqElement R6 base class
#'
#' Internal R6 generator for the base class of every drawable element.
#' Subclasses override `prep()` and `draw()`. Not user-facing — concrete
#' elements expose snake_case constructors (e.g. `seq_point()`) defined in
#' their own files.
#'
#' @keywords internal
SeqElementR6 <- R6::R6Class("SeqElement",
  public = list(
    #' @field data GRanges — element's own data (overrides track data when set).
    data        = NULL,
    #' @field mapping SeqMap — element's own mapping (merged field-by-field
    #'   with the parent track's mapping).
    mapping     = NULL,
    #' @field aesthetics SeqAes — constant aesthetics applied to all glyphs.
    aesthetics  = NULL,
    #' @field resolved Named list populated by `resolve()`. Includes resolved
    #'   mapping fields and `.data` / `.mapping` for downstream use in `prep()`.
    resolved    = NULL,
    #' @field coordCanvas List populated by `prep()` with canvas npc
    #'   coordinates ready for `grid` drawing primitives.
    coordCanvas = NULL,
    #' @field legend A single `LegendKey`, a named list of `LegendKey` objects,
    #'   or `NULL`. Defines the legend entry/entries contributed by this element.
    legend      = NULL,
    #' @field show_legend Logical. When `FALSE`, this element contributes no keys
    #'   to any legend regardless of the `legend` field. Default `TRUE`.
    show_legend = TRUE,
    #' @field auto_legend Auto-generated legend spec (set by `prep()`). A
    #'   `SeqLegendSpec`, `GradientLegendSpec`, list of those, or `NULL`.
    #'   Only consulted when `legend` is `NULL`. Users should not set this
    #'   directly — use the `legend` field instead.
    auto_legend = NULL,

    #' @description Construct a new SeqElementR6.
    #' @param data Optional `GRanges` for the element.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`. Defaults to an empty `aes()`.
    #' @param legend Optional `LegendKey` or list of `LegendKey` objects.
    #' @param show_legend Logical. Set to `FALSE` to suppress legend output.
    #' @param ... Unused — accepted so subclasses can pass extra arguments.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          legend = NULL, show_legend = TRUE, ...) {
      # Allow seq_foo(map(...)) shorthand: a SeqMap passed positionally as
      # `data` is interpreted as the mapping.
      if (inherits(data, "SeqMap") && is.null(mapping)) {
        mapping <- data
        data    <- NULL
      }
      self$data        <- data
      self$mapping     <- mapping
      self$aesthetics  <- if (is.null(aesthetics)) aes() else aesthetics
      self$legend      <- legend
      self$show_legend <- show_legend
    },

    #' @description Resolve the effective data + mapping for this element by
    #'   merging in the parent track's defaults. Element fields take priority;
    #'   missing fields are inherited from the track.
    #' @param track_data Optional `GRanges` from the parent track.
    #' @param track_mapping Optional `SeqMap` from the parent track.
    #' @return The element, invisibly.
    resolve = function(track_data = NULL, track_mapping = NULL) {
      eff_data <- self$data %||% track_data

      eff_mapping <- if (is.null(self$mapping) && is.null(track_mapping)) {
        NULL
      } else if (is.null(self$mapping)) {
        track_mapping
      } else if (is.null(track_mapping)) {
        self$mapping
      } else {
        merged <- as.list(track_mapping)
        for (nm in names(self$mapping)) merged[[nm]] <- self$mapping[[nm]]
        structure(merged, class = "SeqMap")
      }

      self$resolved          <- .resolve_mapping(eff_data, eff_mapping)
      self$resolved$.data    <- eff_data
      self$resolved$.mapping <- eff_mapping

      # Extract axis selectors (`axis.x`, `axis.y`) from resolved mapping.
      # These are per-element, not per-observation — reduce to a single
      # integer in {1, 2}. Keys carry dots and are stored verbatim in the
      # SeqMap list, so we look them up by the dotted name.
      ax <- self$resolved[["axis.x"]]
      ay <- self$resolved[["axis.y"]]
      self$resolved$axis_x <- .normalize_axis_selector(ax, which = "x")
      self$resolved$axis_y <- .normalize_axis_selector(ay, which = "y")
      invisible(self)
    },

    #' @description Override in subclasses. Default implementation errors.
    #' @param layout_track The current track's panel metadata list.
    #' @param track_windows The current track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      stop("prep() must be implemented by ", class(self)[1], call. = FALSE)
    },

    #' @description Override in subclasses. Default implementation errors.
    draw = function() {
      stop("draw() must be implemented by ", class(self)[1], call. = FALSE)
    },

    #' @description Infer a `SeqPositionScale` for the y-axis from the
    #'   resolved `y` mapping. Used by `SeqPlotR6$layoutGrid()` before any
    #'   `prep()` call so that `.compute_track_yscale()` has something to work
    #'   with. Returns `NULL` when no usable y vector is available.
    .infer_scale_y = function() {
      if (is.null(self$resolved) || is.null(self$resolved$y)) return(NULL)
      y <- self$resolved$y
      if (!is.numeric(y) || length(y) == 0) return(NULL)
      seq_scale_continuous(limits = range(y, na.rm = TRUE))
    },

    #' @description
    #' Collect legend keys contributed by this element.
    #' Returns a list of entries, each a named list with fields:
    #'   \describe{
    #'     \item{title}{Character or `NULL`. The legend group title, taken from
    #'       the `LegendKey`'s `title` field when present.}
    #'     \item{key}{A `LegendKey` object.}
    #'     \item{element_class}{Character. The R6 class name of the contributing
    #'       element.}
    #'   }
    #' Returns `NULL` when `show_legend` is `FALSE` or `legend` is `NULL`.
    #' @return A list of entries, or `NULL`.
    collect_legend_keys = function() {
      if (!isTRUE(self$show_legend)) return(NULL)
      if (is.null(self$legend))      return(NULL)

      cls <- class(self)[1]
      leg <- self$legend

      if (inherits(leg, "LegendKey")) {
        return(list(list(
          title         = leg$title,
          key           = leg,
          element_class = cls
        )))
      }

      if (inherits(leg, "SeqLegendSpec")) {
        keys <- leg$keys
        if (length(keys) == 0L) return(NULL)
        out <- vector("list", length(keys))
        for (i in seq_along(keys)) {
          k <- keys[[i]]
          if (!inherits(k, "LegendKey"))
            stop(sprintf("legend keys[[%d]] must be a LegendKey.", i))
          out[[i]] <- list(
            title         = k$title %||% leg$title,
            key           = k,
            element_class = cls
          )
        }
        return(out)
      }

      if (is.list(leg)) {
        out <- vector("list", length(leg))
        for (i in seq_along(leg)) {
          k <- leg[[i]]
          if (!inherits(k, "LegendKey")) {
            stop(sprintf(
              "legend[[%d]] must be a LegendKey object (got %s).",
              i, paste(class(k), collapse = "/")
            ))
          }
          out[[i]] <- list(
            title         = k$title,
            key           = k,
            element_class = cls
          )
        }
        return(out)
      }

      stop("legend must be a LegendKey or a list of LegendKey objects.")
    }
  )
)
