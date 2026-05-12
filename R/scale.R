# ── Internal aesthetic resolver ────────────────────────────────────────────────

#' Resolve a SeqAes + SeqScale into concrete per-observation vectors
#'
#' @param data_mcols  data.frame of GRanges metadata columns.
#' @param aes_obj     A `SeqAes` object (from `aes()`), or `NULL`.
#' @param scale_obj   A `SeqScale` object (from `seq_scale_*()`), or `NULL`.
#' @param n           Number of observations.
#' @param default_color Fallback color when no mapping is specified.
#' @return Named list with resolved vectors for `color`, `fill`, `alpha`,
#'   `size`, `shape` (only those present in `aes_obj`).
#'
#' @keywords internal
.resolve_aes <- function(data_mcols, aes_obj, scale_obj, n,
                         default_color = "#1C1B1A") {
  result <- list(color = rep(default_color, n))
  if (is.null(aes_obj)) return(result)

  for (aes_name in c("color", "fill", "alpha", "size", "shape")) {
    col_name <- aes_obj[[aes_name]]
    if (is.null(col_name) || !col_name %in% names(data_mcols)) next
    raw <- data_mcols[[col_name]]
    sc  <- if (!is.null(scale_obj) && scale_obj$aesthetic == aes_name) scale_obj else NULL

    if (aes_name %in% c("color", "fill")) {
      is_discrete <- is.factor(raw) || is.character(raw) ||
                     (!is.null(sc) && sc$type == "discrete")

      if (is_discrete) {
        raw_f  <- as.factor(raw)
        lvls   <- levels(raw_f)
        pal    <- if (!is.null(sc) && !is.null(sc$values)) {
          sc$values
        } else if (!is.null(sc) && !is.null(sc$palette)) {
          sc$palette(length(lvls))
        } else {
          flexoki_palette(length(lvls))
        }
        names(pal)              <- lvls[seq_along(pal)]
        resolved                <- pal[as.character(raw_f)]
        resolved[is.na(resolved)] <- if (!is.null(sc)) sc$na_value %||% "grey80" else "grey80"

      } else {
        raw_n  <- as.numeric(raw)
        lims   <- if (!is.null(sc) && !is.null(sc$limits)) sc$limits
                  else range(raw_n, na.rm = TRUE)
        pal_nm <- if (!is.null(sc)) sc$palette %||% "viridis" else "viridis"
        stops  <- switch(pal_nm,
                    viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
                    plasma  = c("#0d0887", "#cc4778", "#f0f921"),
                    magma   = c("#000004", "#b63679", "#fcfdbf"),
                    blues   = c("#f7fbff", "#2171b5", "#08306b"),
                    reds    = c("#fff5f0", "#ef3b2c", "#67000d"),
                    stop("Unknown palette: ", pal_nm))
        ramp   <- grDevices::colorRamp(stops)
        t_vals <- pmax(0, pmin(1, (raw_n - lims[1]) / diff(lims)))
        cols   <- ramp(t_vals)
        resolved <- grDevices::rgb(cols[, 1], cols[, 2], cols[, 3],
                                   maxColorValue = 255)
        resolved[is.na(raw_n)] <- if (!is.null(sc)) sc$na_value %||% "grey80" else "grey80"
      }

      result[[aes_name]] <- unname(resolved)

    } else {
      # size, alpha, shape — numeric pass-through
      result[[aes_name]] <- as.numeric(raw)
    }
  }
  result
}

# ── Aesthetic scale constructors ──────────────────────────────────────────────

#' Continuous color scale for SeqAes mappings
#'
#' @param palette One of `"viridis"`, `"plasma"`, `"magma"`, `"blues"`, `"reds"`.
#' @param limits  Optional numeric vector of length 2 clamping the scale range.
#' @param na_value Color for NA values (default `"grey80"`).
#' @param midpoint Optional midpoint for diverging scales (not yet implemented).
#' @param ... Arguments forwarded from `seq_scale_fill_continuous()` to
#'   `seq_scale_color_continuous()`.
#' @return A `SeqScaleContinuous` / `SeqScale` object.
#' @export
seq_scale_color_continuous <- function(palette  = "viridis",
                                       limits   = NULL,
                                       na_value = "grey80",
                                       midpoint = NULL) {
  structure(
    list(aesthetic = "color", type = "continuous",
         palette = palette, limits = limits,
         na_value = na_value, midpoint = midpoint),
    class = c("SeqScaleContinuous", "SeqScale")
  )
}

#' Discrete color scale for SeqAes mappings
#'
#' @param values  Optional named character vector mapping level names to colors.
#' @param palette Optional function `function(n)` returning `n` hex colors.
#'   Falls back to `flexoki_palette()` if `NULL`.
#' @param na_value Color for NA / unmatched values (default `"grey80"`).
#' @param ... Arguments forwarded from `seq_scale_fill_discrete()` to
#'   `seq_scale_color_discrete()`.
#' @return A `SeqScaleDiscrete` / `SeqScale` object.
#' @export
seq_scale_color_discrete <- function(values   = NULL,
                                     palette  = NULL,
                                     na_value = "grey80") {
  structure(
    list(aesthetic = "color", type = "discrete",
         values = values, palette = palette, na_value = na_value),
    class = c("SeqScaleDiscrete", "SeqScale")
  )
}

#' @rdname seq_scale_color_continuous
#' @export
seq_scale_fill_continuous <- function(...) {
  s <- seq_scale_color_continuous(...); s$aesthetic <- "fill"; s
}

#' @rdname seq_scale_color_discrete
#' @export
seq_scale_fill_discrete <- function(...) {
  s <- seq_scale_color_discrete(...); s$aesthetic <- "fill"; s
}

# ── Position scale constructors ──────────────────────────────────────────────

#' Genomic position scale
#'
#' Creates a position scale based on genomic coordinates from a `GRanges` object.
#' Used for axes that represent physical genome positions (bp, kb, Mb).
#'
#' @param windows A `GRanges` object defining the genomic windows.
#' @param scale_factor Numeric vector of per-window scale factors controlling
#'   the unit label (1e-6 = Mb, 1e-3 = kb, 1 = bp). If `NULL`, reads from
#'   `mcols(windows)$scale` or defaults to `1e-6`.
#' @param breaks Optional numeric vector of explicit break positions. When
#'   supplied, these are used instead of `pretty()`-generated breaks.
#' @param minor_breaks Optional scalar (number of sub-divisions between
#'   major breaks) or numeric vector of explicit minor-break positions.
#' @param expand Length-2 `c(mul, add)` specifying multiplicative and
#'   additive padding around the data range. Default `c(0, 0)` for
#'   genomic (windows already set the visible range).
#' @param cap One of `"capped"` (axis line spans the break range),
#'   `"full"` (spans the expanded plot range), `"exact"` (spans the
#'   unexpanded data range), or `"ticks"` (no axis line, only ticks).
#' @param labels Optional character vector of tick labels (same length as
#'   `breaks`). If `NULL`, breaks are formatted as decimal numbers.
#' @return A `SeqScaleGenomic` / `SeqPositionScale` object.
#' @export
seq_scale_genomic <- function(windows, scale_factor = NULL,
                              breaks = NULL, minor_breaks = NULL,
                              expand = c(0, 0),
                              cap = c("capped", "full", "exact", "ticks"),
                              labels = NULL) {
  stopifnot(inherits(windows, "GRanges"))
  cap <- match.arg(cap)
  if (is.null(scale_factor)) {
    scale_factor <- if ("scale" %in% names(S4Vectors::mcols(windows)))
      S4Vectors::mcols(windows)$scale
    else
      rep(1e-6, length(windows))
  }
  structure(
    list(type = "genomic", windows = windows, scale_factor = scale_factor,
         breaks = breaks, minor_breaks = minor_breaks,
         expand = expand, cap = cap, labels = labels),
    class = c("SeqScaleGenomic", "SeqPositionScale")
  )
}

#' Continuous position scale
#'
#' Creates a numeric position scale for axes displaying continuous data.
#'
#' @param limits Optional numeric vector of length 2 clamping the axis range.
#'   If `NULL`, the range is auto-computed from element data.
#' @param n_breaks Target number of pretty breaks (default 5). Ignored when
#'   `breaks` is supplied.
#' @param breaks Optional numeric vector of explicit break positions.
#' @param minor_breaks Optional scalar (number of sub-divisions between
#'   major breaks) or numeric vector of explicit minor-break positions.
#' @param expand Length-2 `c(mul, add)` specifying padding around the
#'   data range. `mul` is multiplied by `diff(data_range)`; `add` is
#'   added in data units. Default `c(0.025, 0)` — a 2.5% breath.
#' @param cap One of `"capped"` (default — axis line spans the break
#'   range, so tick labels don't look stranded), `"full"` (spans the
#'   expanded plot range), `"exact"` (spans the unexpanded data range),
#'   or `"ticks"` (suppress the axis line entirely).
#' @param labels Optional character vector of tick labels (same length
#'   as `breaks`). When `NULL`, breaks are formatted automatically.
#' @return A `SeqScaleContinuous_Pos` / `SeqPositionScale` object.
#' @export
seq_scale_continuous <- function(limits = NULL, n_breaks = 5,
                                 breaks = NULL, minor_breaks = NULL,
                                 expand = c(0.025, 0),
                                 cap = c("capped", "full", "exact", "ticks"),
                                 labels = NULL) {
  cap <- match.arg(cap)
  structure(
    list(type = "continuous", limits = limits, n_breaks = n_breaks,
         breaks = breaks, minor_breaks = minor_breaks,
         expand = expand, cap = cap, labels = labels),
    class = c("SeqScaleContinuous_Pos", "SeqPositionScale")
  )
}

#' Discrete position scale
#'
#' Creates a categorical position scale for axes with discrete levels
#' (e.g., cell types, sample names).
#'
#' @param levels Character vector of category levels, in display order.
#'   If `NULL`, levels are auto-detected from element data.
#' @param labels Optional display labels (same length as `levels`).
#'   If `NULL`, level names are used as labels.
#' @param expand Length-2 `c(mul, add)` expansion. Default `c(0, 0.5)`
#'   (half a category of padding on each side).
#' @param cap Axis-line cap mode. See [seq_scale_continuous()].
#' @return A `SeqScaleDiscrete_Pos` / `SeqPositionScale` object.
#' @export
seq_scale_discrete <- function(levels = NULL, labels = NULL,
                               expand = c(0, 0.5),
                               cap = c("capped", "full", "exact", "ticks")) {
  cap <- match.arg(cap)
  structure(
    list(type = "discrete", levels = levels, labels = labels,
         expand = expand, cap = cap),
    class = c("SeqScaleDiscrete_Pos", "SeqPositionScale")
  )
}

# ── Scale break / expansion helpers ──────────────────────────────────────────

#' Expand a data range by a `c(mul, add)` specification
#'
#' Returns `c(lo - mul*span - add, hi + mul*span + add)` where
#' `span = diff(r)`. Used to introduce a small "breath" around the data
#' so points are not crammed against the panel edge.
#'
#' @param r Length-2 numeric vector (the data range).
#' @param expand Length-1 or length-2 numeric. Scalars become `c(mul, 0)`.
#' @return Length-2 numeric vector.
#' @keywords internal
.expand_limits <- function(r, expand = c(0, 0)) {
  if (is.null(r) || length(r) < 2L || any(!is.finite(r)))
    return(r)
  if (length(expand) == 1L) expand <- c(expand, 0)
  mul <- expand[1] %||% 0
  add <- expand[2] %||% 0
  span <- diff(r)
  c(r[1] - mul * span - add, r[2] + mul * span + add)
}

#' Compute minor breaks between majors
#'
#' If `user_val` is a numeric vector of length > 1, it is returned as-is
#' (filtered to `plot_range`). If it is a single integer `N`, `N`
#' equally-spaced interior subdivisions are placed between adjacent major
#' breaks. If `NULL`, returns `NULL` (no minor breaks).
#'
#' @param user_val User-supplied minor_breaks value, or `NULL`.
#' @param major_breaks Numeric vector of major-break positions.
#' @param plot_range Length-2 numeric vector of the visible plot range.
#' @return Numeric vector of minor-break positions, or `NULL`.
#' @keywords internal
.compute_minor_breaks <- function(user_val, major_breaks, plot_range) {
  if (is.null(user_val)) return(NULL)
  if (length(user_val) > 1L) {
    mb <- as.numeric(user_val)
    return(mb[mb >= plot_range[1] & mb <= plot_range[2]])
  }
  n_sub <- as.integer(user_val)
  if (!is.finite(n_sub) || n_sub <= 0L || length(major_breaks) < 2L)
    return(NULL)
  mb <- unlist(lapply(seq_len(length(major_breaks) - 1L), function(i) {
    a <- major_breaks[i]; b <- major_breaks[i + 1L]
    seq.int(a, b, length.out = n_sub + 2L)[-c(1L, n_sub + 2L)]
  }))
  mb[mb >= plot_range[1] & mb <= plot_range[2]]
}

#' Compute breaks, labels, minor breaks, and axis_range for a scale
#'
#' Given a scale (continuous / genomic / discrete) and the raw data
#' range, returns a list with the full axis-drawing metadata: `breaks`,
#' `labels`, `minor_breaks`, `axis_range` (where the axis line should
#' span, governed by `cap`), `plot_range` (the expanded visible range),
#' and `data_range`.
#'
#' Algorithm ported from THEfunc's `make_axis_meta()`; uses base
#' `pretty()` instead of `scales::pretty_breaks()` so no new dependency
#' is introduced.
#'
#' @param scale A `SeqPositionScale`.
#' @param data_range Length-2 numeric vector of the unexpanded data
#'   range. For genomic scales, this is the window range.
#' @param plot_range Optional pre-expanded plot range. When `NULL`
#'   (default), computed from `data_range` via `.expand_limits()` using
#'   `scale$expand`.
#' @return A list with `breaks`, `labels`, `minor_breaks`, `axis_range`,
#'   `plot_range`, `data_range`.
#' @keywords internal
.compute_scale_breaks <- function(scale, data_range, plot_range = NULL) {
  if (is.null(scale)) scale <- list(type = "continuous")

  expand <- scale$expand %||% c(0, 0)
  if (is.null(plot_range)) plot_range <- .expand_limits(data_range, expand)
  cap <- scale$cap %||% "capped"
  n_breaks <- scale$n_breaks %||% 5

  if (identical(scale$type, "discrete")) {
    lvls <- scale$levels
    br   <- seq_along(lvls)
    lbl  <- scale$labels %||% lvls
    ar <- switch(cap,
      full   = plot_range,
      capped = if (length(br)) range(br) else plot_range,
      exact  = data_range,
      ticks  = NULL)
    return(list(breaks = br, labels = as.character(lbl),
                minor_breaks = NULL,
                axis_range = ar, plot_range = plot_range,
                data_range = data_range))
  }

  # Continuous / genomic path.
  br <- if (!is.null(scale$breaks)) as.numeric(scale$breaks)
        else pretty(if (!is.null(data_range)) data_range else plot_range,
                    n = n_breaks)

  # When 0 appears in the breaks and the data range is strictly positive,
  # snap the lower plot_range to 0 so it is not filtered out — 0 is the
  # natural axis origin for genomic coordinates.
  if (0 %in% br && isTRUE(plot_range[1] > 0))
    plot_range[1] <- 0

  # Trim to the expanded plot range with a small tolerance.
  tol <- max(diff(plot_range) * 1e-9, .Machine$double.eps)
  br <- br[is.finite(br) & br >= (plot_range[1] - tol) &
                          br <= (plot_range[2] + tol)]
  br <- sort(unique(br))

  lbl <- if (!is.null(scale$labels)) {
    if (length(scale$labels) != length(br))
      rep(scale$labels, length.out = length(br))
    else scale$labels
  } else {
    br
  }

  ar <- switch(cap,
    full   = plot_range,
    capped = if (length(br)) range(br) else plot_range,
    exact  = data_range,
    ticks  = NULL)

  mbr <- .compute_minor_breaks(scale$minor_breaks, br, plot_range)

  list(breaks = br, labels = lbl, minor_breaks = mbr,
       axis_range = ar, plot_range = plot_range, data_range = data_range)
}

#' Merge theme-shortcut scale fields into an (optional) scale object
#'
#' The theme can express scale settings via keys like
#' `axis.x1.scale.limits` / `.breaks` / `.expand` / `.cap`. This helper
#' builds or augments the canonical scale object:
#'
#' * If `scale` is `NULL` and the theme has at least one relevant entry,
#'   a `seq_scale_continuous()` is constructed from the theme values.
#' * If `scale` is non-`NULL`, any `NULL` fields on `scale` are filled
#'   from the theme. User-supplied fields always win.
#'
#' @param scale An existing scale or `NULL`.
#' @param axis_spec A per-side axis spec from `.build_axis_spec()`.
#' @return An updated scale, or `NULL` when neither source supplies one.
#' @keywords internal
.merge_scale_with_theme <- function(scale, axis_spec) {
  sc <- axis_spec$scale
  if (is.null(scale)) {
    # Build a continuous scale only if the theme provides limits; without
    # limits the scale isn't useful and we prefer to let inference fill it.
    if (is.null(sc$limits) && is.null(sc$breaks)) return(NULL)
    return(seq_scale_continuous(
      limits       = sc$limits,
      n_breaks     = sc$n_breaks %||% 5,
      breaks       = sc$breaks,
      minor_breaks = sc$minor_breaks,
      expand       = sc$expand %||% c(0.025, 0),
      cap          = sc$cap %||% "capped",
      labels       = sc$labels
    ))
  }

  # Augment: copy in fields only where the scale has a NULL.
  if (is.null(scale$limits)       && !is.null(sc$limits))       scale$limits       <- sc$limits
  if (is.null(scale$breaks)       && !is.null(sc$breaks))       scale$breaks       <- sc$breaks
  if (is.null(scale$minor_breaks) && !is.null(sc$minor_breaks)) scale$minor_breaks <- sc$minor_breaks
  if (is.null(scale$labels)       && !is.null(sc$labels))       scale$labels       <- sc$labels
  if (is.null(scale$expand)       && !is.null(sc$expand))       scale$expand       <- sc$expand
  if (is.null(scale$cap)          && !is.null(sc$cap))          scale$cap          <- sc$cap
  if (is.null(scale$n_breaks)     && !is.null(sc$n_breaks))     scale$n_breaks     <- sc$n_breaks
  scale
}
