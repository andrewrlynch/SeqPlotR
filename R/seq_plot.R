# ── SeqPlotR6 ─────────────────────────────────────────────────────────────────

#' SeqPlot R6 class
#'
#' Internal R6 class backing [seq_plot()]. Public users should go through the
#' snake_case constructor.
#'
#' @keywords internal
SeqPlotR6 <- R6::R6Class("SeqPlot",
  public = list(
    #' @field rows List of lists of SeqTrackR6 — positional layout. `NULL` in
    #'   patchwork mode.
    rows             = NULL,
    #' @field layout_str Raw layout string — patchwork mode (authoritative).
    layout_str       = NULL,
    #' @field aesthetics SeqAes for plot-wide aesthetics.
    aesthetics       = NULL,
    #' @field tracks Flat list of SeqTrackR6 — populated in patchwork mode.
    tracks           = NULL,
    #' @field plot_links Plot-level SeqLink objects — deferred, drawn last.
    plot_links       = NULL,
    #' @field plot_annotations Plot-level SeqAnnotation objects — deferred.
    plot_annotations = NULL,
    #' @field layout Layout metadata produced by `$layoutGrid()` — list with
    #'   `panelBounds` and `trackBounds`.
    layout           = NULL,

    #' @field flat_theme The merged flat theme map (plot-level default theme
    #'   + user overrides). Populated at `layoutGrid()` time.
    flat_theme        = NULL,
    #' @field show_legend Logical. Global legend switch. When `FALSE`, no legend
    #'   output is produced for any track. Default `TRUE`.
    show_legend       = TRUE,

    #' @description Construct a SeqPlotR6.
    #' @param layout Either NULL (positional layout, default) or a multiline
    #'   character string defining a patchwork layout.
    #' @param aesthetics A SeqAes object.
    #' @param tracks Optional list of `SeqTrackR6` objects to pre-populate the
    #'   plot. Elements can also be added with `%+%`. When supplied with no
    #'   `layout` string, tracks are placed in one positional row.
    #' @param show_legend Logical. Global legend switch. Default `TRUE`.
    #' @param legend Convenience argument. `legend = FALSE` is sugar for
    #'   `show_legend = FALSE`. All other values are ignored.
    #' @param ... Additional arguments (currently ignored).
    initialize = function(layout = NULL, aesthetics = aes(),
                          tracks = list(), show_legend = TRUE,
                          legend = NULL, ...) {
      self$aesthetics       <- aesthetics
      self$plot_links       <- list()
      self$plot_annotations <- list()
      self$show_legend      <- if (identical(legend, FALSE)) FALSE
                               else isTRUE(show_legend)
      if (is.null(layout)) {
        # Pre-populate positional rows from the tracks argument if supplied.
        if (length(tracks) > 0L) {
          self$rows       <- list(tracks)
        } else {
          self$rows       <- list(list())
        }
        self$layout_str <- NULL
        self$tracks     <- NULL
      } else {
        self$layout_str <- layout
        self$tracks     <- if (length(tracks) > 0L) tracks else list()
        self$rows       <- NULL
      }
    },

    #' @description Add a track. Uses `direction` in positional mode; appends
    #'   to the flat tracks list in patchwork mode.
    #' @param track A SeqTrackR6 instance.
    addTrack = function(track) {
      if (!is.null(self$layout_str)) {
        self$tracks <- append(self$tracks, list(track))
        return(invisible(self))
      }
      dir <- track$direction %||% "right"
      if (length(self$rows) == 0 || length(self$rows[[1]]) == 0) {
        if (length(self$rows) == 0) self$rows <- list(list())
        self$rows[[1]] <- append(self$rows[[1]], list(track))
      } else if (dir == "right") {
        n <- length(self$rows)
        self$rows[[n]] <- append(self$rows[[n]], list(track))
      } else {
        self$rows <- append(self$rows, list(list(track)))
      }
      invisible(self)
    },

    #' @description Return all tracks as a flat list regardless of mode.
    allTracks = function() {
      if (!is.null(self$layout_str)) return(self$tracks)
      unlist(self$rows, recursive = FALSE)
    },

    #' @description Return a character vector of every track's `track_id` in
    #'   the order tracks were added. Tracks without a `track_id` contribute
    #'   `NA_character_`.
    trackIds = function() {
      vapply(self$allTracks(),
             function(t) t$track_id %||% NA_character_,
             character(1))
    },

    #' @description Compute the layout grid for all tracks and windows.
    #'   Errors if any track has `windows = NULL` (every track must have
    #'   windows defined in SeqPlotR — there are no global plot-level windows).
    #'   Builds `scale_x` from `seq_scale_genomic(windows)` when missing, and
    #'   auto-infers `scale_y` from element data when missing. Then dispatches
    #'   to either `.build_positional_layout()` or
    #'   `.build_patchwork_layout()` and opens a `grid` viewport.
    #' @return The plot, invisibly.
    layoutGrid = function() {
      plot_theme_user <- .flatten_theme(self$aesthetics)
      plot_flat <- .merge_themes(.default_theme(), plot_theme_user)
      self$flat_theme <- plot_flat

      all_trks <- self$allTracks()
      for (trk in all_trks) {
        if (is.null(trk$windows))
          stop("seq_track() is missing 'windows'. ",
               "Every track must have windows defined.", call. = FALSE)
      }

      # Per-track resolved theme (plot + track-level overrides).
      for (trk in all_trks) {
        trk_flat <- .merge_themes(plot_flat,
                                  .flatten_theme(trk$aesthetics))
        trk$resolved_theme <- .build_resolved_theme(trk_flat)
      }

      # Auto-detect axis orientation from the track's mapping. A bare
      # `y = start / end / mid / width` mapping signals a genomic y-axis;
      # a non-genomic `x` mapping (with no explicit scale_x) signals a
      # scalar x-axis whose range is inferred from element data.
      for (trk in all_trks) {
        mx <- trk$mapping$x
        my <- trk$mapping$y

        y_is_genomic <- .is_genomic_special(my)
        x_is_genomic <- is.null(mx) || .is_genomic_special(mx)

        if (y_is_genomic && is.null(trk$scale_y) &&
            !isTRUE(trk$uses_genomic_y)) {
          y_win <- trk$y_windows %||% trk$windows
          trk$scale_y       <- seq_scale_genomic(y_win)
          trk$uses_genomic_y <- TRUE
          trk$y_windows     <- y_win
        }

        if (!x_is_genomic && is.null(trk$scale_x)) {
          x_rng <- .infer_x_range(trk)
          if (!is.null(x_rng))
            trk$scale_x <- seq_scale_continuous(limits = x_rng)
        } else if (inherits(trk$scale_x, "SeqScaleContinuous_Pos") &&
                   is.null(trk$scale_x$limits)) {
          # User passed an explicit continuous scale_x but left limits NULL;
          # infer the range from element data so the axis doesn't collapse
          # to the c(0, 1) placeholder.
          x_rng <- .infer_x_range(trk)
          if (!is.null(x_rng)) trk$scale_x$limits <- x_rng
        }
      }

      # Any track that still has no scale_x falls back to genomic x.
      for (trk in all_trks) {
        if (is.null(trk$scale_x))
          trk$scale_x <- seq_scale_genomic(trk$windows)
      }

      # Resolve elements once so axis_y / axis_x selectors are known, then
      # infer primary and secondary scales from the elements assigned to
      # each. Primary scales take existing inference logic; secondary is
      # inferred from elements whose `map(axis.y=2)` / `map(axis.x=2)`
      # routes them to the secondary axis.
      for (trk in all_trks) {
        for (elem in trk$elements) {
          tryCatch(
            elem$resolve(track_data    = trk$data,
                         track_mapping = trk$mapping),
            error = function(e) NULL
          )
        }
        prim_y_elems <- Filter(function(e)
          isTRUE(e$resolved$axis_y %||% 1L) == 1L ||
          identical(e$resolved$axis_y, 1L), trk$elements)
        sec_y_elems  <- Filter(function(e)
          identical(e$resolved$axis_y, 2L), trk$elements)
        sec_x_elems  <- Filter(function(e)
          identical(e$resolved$axis_x, 2L), trk$elements)

        # Helper: pick the first element-inferred y limits, or NULL.
        .first_inferred_y_limits <- function(elems) {
          for (elem in elems) {
            inf <- tryCatch(elem$.infer_scale_y(), error = function(e) NULL)
            if (!is.null(inf) && !is.null(inf$limits)) return(inf$limits)
          }
          NULL
        }

        if (is.null(trk$scale_y) && length(prim_y_elems) > 0L) {
          for (elem in prim_y_elems) {
            inf <- tryCatch(elem$.infer_scale_y(),
                            error = function(e) NULL)
            if (!is.null(inf)) { trk$scale_y <- inf; break }
          }
        } else if (inherits(trk$scale_y, "SeqScaleContinuous_Pos") &&
                   is.null(trk$scale_y$limits) &&
                   length(prim_y_elems) > 0L) {
          y_rng <- .first_inferred_y_limits(prim_y_elems)
          if (!is.null(y_rng)) trk$scale_y$limits <- y_rng
        }
        if (is.null(trk$scale_y2) && length(sec_y_elems) > 0L) {
          for (elem in sec_y_elems) {
            inf <- tryCatch(elem$.infer_scale_y(),
                            error = function(e) NULL)
            if (!is.null(inf)) { trk$scale_y2 <- inf; break }
          }
        } else if (inherits(trk$scale_y2, "SeqScaleContinuous_Pos") &&
                   is.null(trk$scale_y2$limits) &&
                   length(sec_y_elems) > 0L) {
          y_rng <- .first_inferred_y_limits(sec_y_elems)
          if (!is.null(y_rng)) trk$scale_y2$limits <- y_rng
        }
        sec_x_range <- function() {
          xs <- unlist(lapply(sec_x_elems, function(e) {
            x <- e$resolved$x
            if (is.numeric(x) && length(x)) range(x, na.rm = TRUE) else NULL
          }))
          if (length(xs) >= 2L && all(is.finite(xs))) range(xs) else NULL
        }
        if (is.null(trk$scale_x2) && length(sec_x_elems) > 0L) {
          x_rng <- sec_x_range()
          if (!is.null(x_rng))
            trk$scale_x2 <- seq_scale_continuous(limits = x_rng)
        } else if (inherits(trk$scale_x2, "SeqScaleContinuous_Pos") &&
                   is.null(trk$scale_x2$limits) &&
                   length(sec_x_elems) > 0L) {
          x_rng <- sec_x_range()
          if (!is.null(x_rng)) trk$scale_x2$limits <- x_rng
        }

        # Merge theme shortcuts (axis.<side>.scale.*) into each scale.
        trk$scale_x  <- .merge_scale_with_theme(trk$scale_x,
                           trk$resolved_theme$axes$x1)
        trk$scale_y  <- .merge_scale_with_theme(trk$scale_y,
                           trk$resolved_theme$axes$y1)
        trk$scale_x2 <- .merge_scale_with_theme(trk$scale_x2,
                           trk$resolved_theme$axes$x2)
        trk$scale_y2 <- .merge_scale_with_theme(trk$scale_y2,
                           trk$resolved_theme$axes$y2)

        # Flag which secondary axes are "active" for draw decisions.
        visible_x2 <- trk$resolved_theme$axes$x2$visible
        visible_y2 <- trk$resolved_theme$axes$y2$visible
        trk$has_axis_x2 <- length(sec_x_elems) > 0L ||
                           !is.null(trk$scale_x2) ||
                           isTRUE(visible_x2)
        trk$has_axis_y2 <- length(sec_y_elems) > 0L ||
                           !is.null(trk$scale_y2) ||
                           isTRUE(visible_y2)
      }

      if (!is.null(self$layout_str)) {
        self$layout <- .build_patchwork_layout(
          self$tracks, self$layout_str, plot_flat
        )
      } else {
        self$layout <- .build_positional_layout(self$rows, plot_flat)
      }

      margins_m <- self$flat_theme$margins %||%
                   list(top = 0, right = 0, bottom = 0, left = 0)
      self$layout$trackMarginBounds <- lapply(
        self$layout$trackBounds, function(tb) {
          list(
            top    = list(x0 = tb$x0, x1 = tb$x1,
                          y0 = 1 - margins_m$top, y1 = 1),
            bottom = list(x0 = tb$x0, x1 = tb$x1,
                          y0 = 0, y1 = margins_m$bottom),
            left   = list(x0 = 0, x1 = margins_m$left,
                          y0 = tb$y0, y1 = tb$y1),
            right  = list(x0 = 1 - margins_m$right, x1 = 1,
                          y0 = tb$y0, y1 = tb$y1)
          )
        }
      )
      self$layout$canvasMarginBounds <- list(
        top    = list(x0 = margins_m$left,       x1 = 1 - margins_m$right,
                      y0 = 1 - margins_m$top,    y1 = 1),
        bottom = list(x0 = margins_m$left,       x1 = 1 - margins_m$right,
                      y0 = 0,                    y1 = margins_m$bottom),
        left   = list(x0 = 0,                    x1 = margins_m$left,
                      y0 = margins_m$bottom,     y1 = 1 - margins_m$top),
        right  = list(x0 = 1 - margins_m$right,  x1 = 1,
                      y0 = margins_m$bottom,     y1 = 1 - margins_m$top)
      )

      grid::grid.newpage()
      grid::pushViewport(grid::viewport(name = "seqplotr_root"))
      invisible(self)
    },

    #' @description Draw track backgrounds / borders and per-window panel
    #'   chrome. Delegates to `.draw_track_chrome()` for each track using
    #'   the track's resolved theme.
    drawGrid = function() {
      stopifnot(is.list(self$layout),
                !is.null(self$layout$panelBounds),
                !is.null(self$layout$trackBounds))

      all_trks <- self$allTracks()
      panelBounds <- self$layout$panelBounds
      trackBounds <- self$layout$trackBounds

      for (i in seq_along(all_trks)) {
        trk <- all_trks[[i]]
        key <- if (!is.null(self$layout_str))
          (trk$track_id %||% NA_character_)
        else
          i
        panels <- panelBounds[[key]]
        if (is.null(panels)) next
        tb <- trackBounds[[key]]
        .draw_track_chrome(trk, panels, tb)
      }
      invisible(self)
    },

    #' @description Draw x and y gridlines at axis break positions for all
    #'   track windows. Gridlines sit after window backgrounds and before
    #'   elements. Enable per axis via `axis.x.gridline = TRUE` (or `=
    #'   aes(color, lwd, lty, alpha)`) and `axis.y.gridline = TRUE` in
    #'   `seq_plot()` or `seq_track()` aesthetics. Styling inherits from
    #'   `axis.gridline.*` defaults in the theme hierarchy.
    #' @return Renders gridlines to the graphics device; returns invisibly.
    drawGridlines = function() {
      stopifnot(is.list(self$layout),
                !is.null(self$layout$panelBounds))

      # Decide whether gridlines are enabled for a given dimension ("x" or "y").
      # Three-way check covers:
      #   axis.x.gridline = TRUE            (leaf boolean — most common)
      #   axis.x.gridline = aes(color, ...) (flattens to axis.x.gridline.* sub-keys)
      #   axis.x.gridline.visible = TRUE    (explicit visible sub-key)
      # The first matching form wins; default is FALSE (off).
      .gridline_enabled <- function(flat, dim) {
        root_key <- paste0("axis.", dim, ".gridline")
        root_val <- .resolve_theme(flat, root_key, NULL)
        if (!is.null(root_val) && is.logical(root_val)) return(root_val)
        sub_prefix <- paste0(root_key, ".")
        if (any(startsWith(names(flat), sub_prefix))) return(TRUE)
        isTRUE(.resolve_theme(flat, paste0(root_key, ".visible"), FALSE))
      }

      all_trks    <- self$allTracks()
      panelBounds <- self$layout$panelBounds

      for (i in seq_along(all_trks)) {
        trk <- all_trks[[i]]
        key <- if (!is.null(self$layout_str))
          (trk$track_id %||% NA_character_)
        else
          i
        panels <- panelBounds[[key]]
        if (is.null(panels) || length(panels) == 0L) next

        flat <- trk$resolved_theme$flat %||% list()
        draw_x <- .gridline_enabled(flat, "x")
        draw_y <- .gridline_enabled(flat, "y")
        if (!draw_x && !draw_y) next

        for (win in panels) {
          p <- win$inner

          if (draw_x && !is.null(win$xscale) && !is.null(trk$scale_x)) {
            meta <- .compute_scale_breaks(trk$scale_x, win$xscale,
                                          plot_range = win$xscale)
            gp_x <- grid::gpar(
              col   = .resolve_theme(flat, "axis.x.gridline.color", "grey85"),
              lwd   = .resolve_theme(flat, "axis.x.gridline.lwd", 0.5),
              lty   = .resolve_theme(flat, "axis.x.gridline.lty", 1),
              alpha = .resolve_theme(flat, "axis.x.gridline.alpha", 1)
            )
            for (b in meta$breaks) {
              xpos <- .axis_map_npc(b, win$xscale, p$x0, p$x1)
              grid::grid.lines(x = grid::unit(c(xpos, xpos), "npc"),
                               y = grid::unit(c(p$y0, p$y1), "npc"),
                               gp = gp_x)
            }
          }

          if (draw_y && !is.null(win$yscale) && !is.null(trk$scale_y)) {
            meta <- .compute_scale_breaks(trk$scale_y, win$yscale,
                                          plot_range = win$yscale)
            gp_y <- grid::gpar(
              col   = .resolve_theme(flat, "axis.y.gridline.color", "grey85"),
              lwd   = .resolve_theme(flat, "axis.y.gridline.lwd", 0.5),
              lty   = .resolve_theme(flat, "axis.y.gridline.lty", 1),
              alpha = .resolve_theme(flat, "axis.y.gridline.alpha", 1)
            )
            for (b in meta$breaks) {
              ypos <- .axis_map_npc(b, win$yscale, p$y0, p$y1)
              grid::grid.lines(x = grid::unit(c(p$x0, p$x1), "npc"),
                               y = grid::unit(c(ypos, ypos), "npc"),
                               gp = gp_y)
            }
          }
        }
      }

      invisible()
    },

    #' @description Draw x and y axes for every track. Delegates to
    #'   `.draw_track_axes()`, which reads the track's `resolved_theme`
    #'   and renders up to four axes (x1, x2, y1, y2) with hierarchical
    #'   aesthetic control, break computation, and cap modes.
    drawAxes = function() {
      if (is.null(self$layout)) return(invisible(self))

      all_trks <- self$allTracks()
      panelBounds <- self$layout$panelBounds

      for (i in seq_along(all_trks)) {
        trk <- all_trks[[i]]
        key <- if (!is.null(self$layout_str))
          (trk$track_id %||% NA_character_)
        else
          i
        panels <- panelBounds[[key]]
        if (is.null(panels) || length(panels) == 0L) next
        .draw_track_axes(trk, panels)
      }
      invisible(self)
    },

    #' @description Draw all elements in drawing order: within-track
    #'   non-links, then within-track links, then plot-level deferred links,
    #'   then plot-level annotations.
    drawElements = function() {
      all_trks <- self$allTracks()
      if (length(all_trks) == 0) return(invisible(self))

      keys <- names(self$layout$panelBounds) %||% seq_along(all_trks)

      layout_all_tracks <- setNames(
        self$layout$panelBounds,
        vapply(all_trks, function(t) t$track_id %||% NA_character_,
               character(1))
      )
      track_windows_list <- setNames(
        lapply(all_trks, function(t) t$windows),
        vapply(all_trks, function(t) t$track_id %||% NA_character_,
               character(1))
      )

      # When a track was laid out with combine_windows = TRUE, swap in
      # the virtualized data and the combined window so element prep
      # logic sees a single-window panel with concatenated coordinates.
      .effective_track_io <- function(track, layout_panels) {
        vmap <- if (length(layout_panels)) layout_panels[[1]]$virtual_map_x
                else NULL
        if (is.null(vmap) || !inherits(track$data, "GRanges"))
          return(list(data = track$data, windows = track$windows))
        list(data    = .virtualize_granges(track$data, vmap),
             windows = vmap$combined_window)
      }

      # Inject parent-track data + mapping into every panel of
      # `layout_all_tracks` so that link prep methods can resolve their
      # mappings against the referenced track's data when the link itself
      # carries no `data` of its own.
      for (i in seq_along(all_trks)) {
        track <- all_trks[[i]]
        key <- if (!is.null(self$layout_str))
          (track$track_id %||% NA_character_)
        else
          i
        if (is.null(layout_all_tracks[[key]])) next
        eff <- .effective_track_io(track, layout_all_tracks[[key]])
        for (p in seq_along(layout_all_tracks[[key]])) {
          layout_all_tracks[[key]][[p]]$track_data    <- eff$data
          layout_all_tracks[[key]][[p]]$track_mapping <- track$mapping
        }
        track_windows_list[[key]] <- eff$windows
      }

      # 1. within-track non-link elements
      for (i in seq_along(all_trks)) {
        track <- all_trks[[i]]
        layout_track_key <- if (!is.null(self$layout_str))
          (track$track_id %||% NA_character_)
        else
          i
        layout_track <- self$layout$panelBounds[[layout_track_key]]
        if (is.null(layout_track)) next  # patchwork: skipped track
        eff <- .effective_track_io(track, layout_track)
        for (p in seq_along(layout_track)) {
          layout_track[[p]]$track_data    <- eff$data
          layout_track[[p]]$track_mapping <- track$mapping
        }
        for (elem in track$elements) {
          if (inherits(elem, "SeqLink")) next
          # Ensure axis selectors are set before we swap scales.
          tryCatch(
            elem$resolve(track_data    = eff$data,
                         track_mapping = track$mapping),
            error = function(e) NULL
          )
          lt_elem <- .panels_for_element(layout_track, elem$resolved)
          elem$prep(lt_elem, eff$windows)
          elem$draw()
        }
      }

      # Build an axis-aware view of `layout_all_tracks` for a given link.
      .link_layout_view <- function(link) {
        tryCatch(
          link$resolve(track_data    = NULL,
                       track_mapping = NULL),
          error = function(e) NULL
        )
        lapply(layout_all_tracks, function(panels)
          .panels_for_element(panels, link$resolved))
      }

      # 2. within-track links
      for (i in seq_along(all_trks)) {
        track <- all_trks[[i]]
        for (elem in track$elements) {
          if (!inherits(elem, "SeqLink")) next
          la <- .link_layout_view(elem)
          elem$prep(layout_all_tracks  = la,
                    track_windows_list = track_windows_list,
                    plot_track_index   = i)
          elem$draw()
        }
      }

      # 3. plot-level deferred links
      for (lnk in self$plot_links) {
        la <- .link_layout_view(lnk)
        lnk$prep(layout_all_tracks  = la,
                 track_windows_list = track_windows_list,
                 plot_track_index   = NULL)
        lnk$draw()
      }

      # 4. plot-level annotations (placeholder — wrap defensively)
      for (ann in self$plot_annotations) {
        tryCatch(ann$draw(), error = function(e) NULL)
      }

      invisible(self)
    },

    #' @description Draw legends for all tracks. Dispatches on `position` in
    #'   each `SeqLegendSpec` found on element `legend` fields. Bare
    #'   `LegendKey` or list-of-`LegendKey` on an element are automatically
    #'   wrapped in a default `"inside"` spec. Phase 1 handles `"inside"` and
    #'   `"track_margin"` per element; Phase 2 aggregates all
    #'   `"canvas_margin"` specs and draws once per side.
    #'
    #'   Call after `drawElements()`.
    #' @return Renders legends to the graphics device; returns invisibly.
    drawLegends = function() {
      if (!isTRUE(self$show_legend)) return(invisible())
      if (is.null(self$layout))      return(invisible())

      all_trks <- self$allTracks()

      # --- Phase 1: inside and track_margin (per element) ---
      for (i in seq_along(all_trks)) {
        track <- all_trks[[i]]
        if (!isTRUE(track$show_legend)) next

        key <- if (!is.null(self$layout_str))
          (track$track_id %||% NA_character_)
        else
          i
        layout_track <- self$layout$panelBounds[[key]]
        if (is.null(layout_track) || length(layout_track) == 0L) next

        for (elem in track$elements) {
          if (!isTRUE(elem$show_legend)) next

          # Determine which legend(s) to draw: explicit `legend` field first,
          # then the auto-generated `auto_legend` when `legend` is NULL.
          specs_to_draw <- if (!is.null(elem$legend)) {
            list(elem$legend)
          } else if (!is.null(elem$auto_legend)) {
            # auto_legend may be a single spec or a list of specs
            if (is.list(elem$auto_legend) &&
                !inherits(elem$auto_legend, "SeqLegendSpec") &&
                !inherits(elem$auto_legend, "GradientLegendSpec")) {
              elem$auto_legend
            } else {
              list(elem$auto_legend)
            }
          } else {
            list()
          }

          for (spec in specs_to_draw) {
            # Normalise bare LegendKey / list-of-LegendKey to SeqLegendSpec
            if (inherits(spec, "LegendKey")) {
              spec <- seq_legend(spec)
            } else if (is.list(spec) && !inherits(spec, "SeqLegendSpec") &&
                       !inherits(spec, "GradientLegendSpec")) {
              spec <- seq_legend(spec)
            }

            if (inherits(spec, "GradientLegendSpec")) {
              if (spec$position == "inside") {
                .draw_gradient_legend(spec, layout_track[[1]])
              } else if (spec$position == "track_margin") {
                side <- spec$side %||% "right"
                margin_rect <- self$layout$trackMarginBounds[[key]][[side]]
                if (!is.null(margin_rect)) {
                  fake_pm <- list(full = margin_rect, inner = margin_rect)
                  .draw_gradient_legend(spec, fake_pm)
                }
              }
              next
            }

            if (!inherits(spec, "SeqLegendSpec")) next

            if (spec$position == "inside") {
              .draw_legend_inside(spec, layout_track[[1]])

            } else if (spec$position == "track_margin") {
              side <- spec$side %||% "top"
              margin_rect <- self$layout$trackMarginBounds[[key]][[side]]
              if (!is.null(margin_rect)) {
                .draw_legend_track_margin(spec, margin_rect)
              }
            }
            # canvas_margin entries are handled in Phase 2 below
          }
        }
      }

      # --- Phase 2: canvas_margin (aggregated, one draw per side) ---
      canvas_entries <- .collect_canvas_legend_specs(all_trks)

      if (length(canvas_entries) > 0L) {
        for (side in c("top", "bottom", "left", "right")) {
          side_entries <- Filter(function(e) {
            (e$spec$side %||% "top") == side
          }, canvas_entries)

          if (length(side_entries) == 0L) next

          merged_spec <- .merge_canvas_specs(side_entries, side)
          if (is.null(merged_spec)) next

          canvas_rect <- self$layout$canvasMarginBounds[[side]]
          if (!is.null(canvas_rect)) {
            .draw_legend_canvas_margin(merged_spec, canvas_rect)
          }
        }
      }

      invisible()
    },

    #' @description Run the full plot pipeline:
    #'   `layoutGrid()` -> `drawGrid()` -> `drawAxes()` -> `drawElements()`
    #'   -> `drawLegends()`.
    #'   Resets base-graphics `par(mar, oma, mai, omi)` to zero so that any
    #'   residual margins inherited from the active device (knitr's PNG
    #'   device, RStudio's plot pane, etc.) do not leave whitespace around
    #'   the grid viewport.
    plot = function() {
      old_par <- graphics::par(
        mar = c(0, 0, 0, 0),
        oma = c(0, 0, 0, 0),
        mai = c(0, 0, 0, 0),
        omi = c(0, 0, 0, 0)
      )
      on.exit(graphics::par(old_par), add = TRUE)
      self$layoutGrid()
      self$drawGrid()
      self$drawGridlines()
      self$drawAxes()
      self$drawElements()
      self$drawLegends()
      invisible(self)
    }
  )
)

# ── seq_plot() constructor ───────────────────────────────────────────────────

#' Create a new seq_plot
#'
#' @param layout Either NULL (positional layout, default) or a multiline character
#'   string defining a patchwork layout. When a layout string is given, track
#'   positions are determined entirely by `track_id` matching — `direction` on
#'   [seq_track()] is ignored.
#' @param aesthetics A SeqAes object from [aes()] for plot-wide aesthetics.
#' @param show_legend Logical. Global legend switch — when `FALSE` no legend is
#'   drawn regardless of element or track settings. Default `TRUE`.
#' @param legend Convenience alias: `legend = FALSE` is sugar for
#'   `show_legend = FALSE`. All other values are ignored.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqPlotR6` instance (S3 class `"SeqPlot"`).
#' @examples
#' seq_plot()
#' seq_plot(layout = "AB\nCC")
#' @export
seq_plot <- function(layout = NULL, aesthetics = aes(), show_legend = TRUE,
                     legend = NULL, ...) {
  SeqPlotR6$new(layout = layout, aesthetics = aesthetics,
                show_legend = show_legend, legend = legend, ...)
}

#' Auto-print a SeqPlot
#'
#' Renders the plot to the current graphics device. Matches ggplot2's
#' convention that a bare `seq_plot() %+% ...` expression in the console
#' or a knitr chunk draws itself.
#'
#' @param x A `SeqPlot` (i.e. a `SeqPlotR6` instance).
#' @param ... Ignored.
#' @return The plot, invisibly.
#' @export
print.SeqPlot <- function(x, ...) {
  x$plot()
  invisible(x)
}
