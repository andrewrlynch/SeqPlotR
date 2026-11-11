# NODE-LINK ELEMENT

# Node shape name -> grid pch and generic legend shape.
.node_pch <- c(circle = 21L, triangle_up = 24L, triangle_down = 25L,
               diamond = 23L, square = 22L)
.node_key_shape <- c(circle = "circle", triangle_up = "triangle",
                     triangle_down = "triangle", diamond = "diamond",
                     square = "square")

# Build the legend spec of an element from the user argument and auto keys.
# NULL uses the auto keys, FALSE gives no legend, a spec or key list is used as given.
.element_legend_spec <- function(arg, auto_keys, extra, position, side) {
  if (isFALSE(arg)) return(NULL)
  if (inherits(extra, "LegendKey")) extra <- list(extra)
  mk <- function(keys) {
    if (length(keys) == 0L) return(NULL)
    seq_legend(keys, position = position, side = side, x = 0.05, y = 0.95)
  }
  if (inherits(arg, "SeqLegendSpec")) return(arg)
  if (inherits(arg, "LegendKey")) return(mk(c(list(arg), extra)))
  if (is.list(arg)) return(mk(c(arg, extra)))
  mk(c(auto_keys, extra))
}

# Vectorised colour opacity (adjustcolor needs a scalar alpha).
.with_alpha <- function(col, alpha) {
  m <- grDevices::col2rgb(col, alpha = TRUE)
  grDevices::rgb(m[1, ], m[2, ], m[3, ], m[4, ] * alpha, maxColorValue = 255)
}

# Check that a table has the required columns.
.need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop("seq_node_link: ", what, " is missing column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
}

# Panel with nodes, edges and labelled group boxes at given positions.
SeqNodeLinkR6 <- R6::R6Class("SeqNodeLink",
  inherit = SeqElementR6,
  public = list(
    nodes = NULL,
    edges = NULL,
    groups = NULL,
    legend_arg = NULL,
    legend_labels = NULL,
    legend_extra = NULL,
    legend_position = "canvas_margin",
    legend_side = "right",

    initialize = function(nodes = NULL, edges = NULL, groups = NULL,
                          aesthetics = aes(), legend = NULL,
                          legend_labels = NULL, legend_extra = NULL,
                          show_legend = TRUE,
                          legend_position = "canvas_margin",
                          legend_side = "right", ...) {
      nodes  <- if (is.null(nodes))  NULL else as.data.frame(nodes,  stringsAsFactors = FALSE)
      edges  <- if (is.null(edges))  NULL else as.data.frame(edges,  stringsAsFactors = FALSE)
      groups <- if (is.null(groups)) NULL else as.data.frame(groups, stringsAsFactors = FALSE)
      if (!is.null(nodes))  .need_cols(nodes,  c("id", "x", "y"), "nodes")
      if (!is.null(edges))  .need_cols(edges,  c("from", "to"), "edges")
      if (!is.null(groups)) .need_cols(groups, c("x", "y", "label"), "groups")
      if (!is.null(edges) && is.null(nodes))
        stop("seq_node_link: edges need nodes.", call. = FALSE)
      if (is.null(nodes) && is.null(groups))
        stop("seq_node_link: supply nodes or groups.", call. = FALSE)
      if (!is.null(nodes$shape)) {
        bad <- setdiff(unique(nodes$shape[!is.na(nodes$shape)]), names(.node_pch))
        if (length(bad))
          stop("seq_node_link: unknown node shape(s) ", paste(bad, collapse = ", "),
               ".", call. = FALSE)
      }
      super$initialize(data = nodes %||% groups, mapping = map(x = x, y = y),
                       aesthetics = aesthetics, legend = NULL,
                       show_legend = show_legend)
      self$nodes <- nodes
      self$edges <- edges
      self$groups <- groups
      self$legend_arg <- legend
      self$legend_labels <- legend_labels
      self$legend_extra <- legend_extra
      self$legend_position <- legend_position
      self$legend_side <- legend_side
    },

    prep = function(layout_track, track_windows) {
      pm  <- layout_track[[1]]
      xpr <- pm$xplot_range %||% pm$xscale
      ypr <- pm$yplot_range %||% pm$yscale
      xw  <- pm$inner$x1 - pm$inner$x0
      yw  <- pm$inner$y1 - pm$inner$y0
      cx <- function(v) pm$inner$x0 + pmax(pmin((v - xpr[1]) / diff(xpr), 1), 0) * xw
      cy <- function(v) pm$inner$y0 + pmax(pmin((v - ypr[1]) / diff(ypr), 1), 0) * yw
      col_or <- function(df, nm, default) {
        if (is.null(df[[nm]])) rep(default, nrow(df)) else df[[nm]]
      }

      nd <- self$nodes
      nodes_out <- NULL
      if (!is.null(nd) && nrow(nd) > 0L) {
        nodes_out <- data.frame(
          id = as.character(nd$id), x = cx(nd$x), y = cy(nd$y),
          lo = if (is.null(nd$bar_lo)) rep(NA_real_, nrow(nd)) else cy(nd$bar_lo),
          hi = if (is.null(nd$bar_hi)) rep(NA_real_, nrow(nd)) else cy(nd$bar_hi),
          shape = as.character(col_or(nd, "shape", "circle")),
          fill = as.character(col_or(nd, "fill", "grey60")),
          alpha = as.numeric(col_or(nd, "alpha", 1)),
          label = as.character(col_or(nd, "label", NA_character_)),
          stringsAsFactors = FALSE
        )
        nodes_out$shape[is.na(nodes_out$shape)] <- "circle"
      }

      ed <- self$edges
      edges_out <- NULL
      if (!is.null(ed) && nrow(ed) > 0L) {
        i0 <- match(as.character(ed$from), nodes_out$id)
        i1 <- match(as.character(ed$to), nodes_out$id)
        bad <- which(is.na(i0) | is.na(i1))
        if (length(bad))
          stop(sprintf("seq_node_link: edge row %d references an unknown node id (%s -> %s).",
                       bad[1], ed$from[bad[1]], ed$to[bad[1]]), call. = FALSE)
        edges_out <- data.frame(
          from = as.character(ed$from), to = as.character(ed$to),
          x0 = nodes_out$x[i0], y0 = nodes_out$y[i0],
          x1 = nodes_out$x[i1], y1 = nodes_out$y[i1],
          arrow = as.logical(col_or(ed, "arrow", TRUE)),
          color = as.character(col_or(ed, "color", "grey10")),
          linetype = col_or(ed, "linetype", 1),
          linewidth = as.numeric(col_or(ed, "linewidth", 1.4)),
          label = as.character(col_or(ed, "label", NA_character_)),
          flag = as.logical(col_or(ed, "flag", FALSE)),
          stringsAsFactors = FALSE
        )
        edges_out$arrow[is.na(edges_out$arrow)] <- TRUE
        edges_out$flag[is.na(edges_out$flag)] <- FALSE
      }

      gr <- self$groups
      groups_out <- NULL
      if (!is.null(gr) && nrow(gr) > 0L) {
        groups_out <- data.frame(
          x = cx(gr$x), y = cy(gr$y), label = as.character(gr$label),
          fill = as.character(col_or(gr, "fill", "white")),
          color = as.character(col_or(gr, "color", "grey40")),
          stringsAsFactors = FALSE
        )
      }
      self$coordCanvas <- list(nodes = nodes_out, edges = edges_out,
                               groups = groups_out)

      # legend: one key per node shape present, extra keys appended
      auto <- list()
      if (!is.null(nodes_out)) {
        lab <- self$legend_labels
        for (s in intersect(names(.node_pch), unique(nodes_out$shape)))
          auto <- c(auto, list(LegendKey(
            label = if (!is.null(lab[s]) && !is.na(lab[s])) unname(lab[s]) else s,
            shape = unname(.node_key_shape[s]), color = "black", fill = "grey80",
            pch = unname(.node_pch[s])
          )))
      }
      self$legend <- if (isTRUE(self$show_legend))
        .element_legend_spec(self$legend_arg, auto, self$legend_extra,
                             self$legend_position, self$legend_side)
      else NULL
      invisible()
    },

    draw = function() {
      cc <- self$coordCanvas
      if (is.null(cc)) return(invisible())
      a <- self$aesthetics
      node_mm  <- a$node_mm %||% 3
      bar_lw   <- a$bar_linewidth %||% 3
      bar_a    <- a$bar_alpha %||% 0.35
      arrow_mm <- a$arrow_mm %||% 2
      lab_cex  <- a$label_cex %||% 0.6
      elab_cex <- a$edge_label_cex %||% 0.55
      box_r    <- a$box_r_mm %||% 1
      flag_col <- a$flag_color %||% "#D62728"
      nd <- cc$nodes
      ed <- cc$edges

      # node bars
      if (!is.null(nd)) {
        b <- nd[!is.na(nd$lo) & !is.na(nd$hi), , drop = FALSE]
        if (nrow(b) > 0L)
          grid::grid.segments(
            x0 = grid::unit(b$x, "npc"), x1 = grid::unit(b$x, "npc"),
            y0 = grid::unit(b$lo, "npc"), y1 = grid::unit(b$hi, "npc"),
            gp = grid::gpar(col = .with_alpha(b$fill, b$alpha * bar_a),
                            lwd = bar_lw, lineend = "butt"),
            name = "seqplotr.node_link.bars"
          )
      }

      # edges: plain first, then arrowed with the tip pulled back to the node rim
      if (!is.null(ed)) {
        plain <- ed[!ed$arrow, , drop = FALSE]
        if (nrow(plain) > 0L)
          grid::grid.segments(
            x0 = grid::unit(plain$x0, "npc"), x1 = grid::unit(plain$x1, "npc"),
            y0 = grid::unit(plain$y0, "npc"), y1 = grid::unit(plain$y1, "npc"),
            gp = grid::gpar(col = plain$color, lwd = plain$linewidth,
                            lty = plain$linetype),
            name = "seqplotr.node_link.edges"
          )
        arr <- ed[ed$arrow, , drop = FALSE]
        if (nrow(arr) > 0L) {
          ix0 <- grid::convertX(grid::unit(arr$x0, "npc"), "inches", valueOnly = TRUE)
          ix1 <- grid::convertX(grid::unit(arr$x1, "npc"), "inches", valueOnly = TRUE)
          iy0 <- grid::convertY(grid::unit(arr$y0, "npc"), "inches", valueOnly = TRUE)
          iy1 <- grid::convertY(grid::unit(arr$y1, "npc"), "inches", valueOnly = TRUE)
          len <- sqrt((ix1 - ix0)^2 + (iy1 - iy0)^2)
          pull <- pmin((node_mm / 2 + 0.3) / 25.4, len / 2)
          f <- ifelse(len > 0, pull / len, 0)
          grid::grid.segments(
            x0 = grid::unit(ix0, "inches"), y0 = grid::unit(iy0, "inches"),
            x1 = grid::unit(ix1 - f * (ix1 - ix0), "inches"),
            y1 = grid::unit(iy1 - f * (iy1 - iy0), "inches"),
            arrow = grid::arrow(type = "closed", length = grid::unit(arrow_mm, "mm")),
            gp = grid::gpar(col = arr$color, fill = arr$color, lwd = arr$linewidth,
                            lty = arr$linetype),
            name = "seqplotr.node_link.edges_arrow"
          )
        }

        # edge labels and flags at the midpoints
        mx <- (ed$x0 + ed$x1) / 2
        my <- (ed$y0 + ed$y1) / 2
        lab <- !is.na(ed$label)
        if (any(lab))
          grid::grid.text(
            ed$label[lab], x = grid::unit(mx[lab], "npc"),
            y = grid::unit(my[lab], "npc"), just = c("left", "bottom"),
            gp = grid::gpar(cex = elab_cex, col = "grey20"),
            name = "seqplotr.node_link.edge_labels"
          )
        if (any(ed$flag))
          grid::grid.points(
            x = grid::unit(mx[ed$flag], "npc"), y = grid::unit(my[ed$flag], "npc"),
            pch = 23, size = grid::unit(2.5, "mm"),
            gp = grid::gpar(col = flag_col, fill = flag_col),
            name = "seqplotr.node_link.flags"
          )
      }

      # nodes and labels
      if (!is.null(nd)) {
        for (s in intersect(names(.node_pch), unique(nd$shape))) {
          k <- nd[nd$shape == s, , drop = FALSE]
          grid::grid.points(
            x = grid::unit(k$x, "npc"), y = grid::unit(k$y, "npc"),
            pch = .node_pch[[s]], size = grid::unit(node_mm, "mm"),
            gp = grid::gpar(col = "black",
                            fill = .with_alpha(k$fill, k$alpha)),
            name = paste0("seqplotr.node_link.nodes.", s)
          )
        }
        lab <- !is.na(nd$label)
        if (any(lab))
          grid::grid.text(
            nd$label[lab], x = grid::unit(nd$x[lab], "npc") + grid::unit(2.5, "mm"),
            y = grid::unit(nd$y[lab], "npc"), just = c("left", "center"),
            gp = grid::gpar(cex = lab_cex, col = "#1C1B1A"),
            name = "seqplotr.node_link.node_labels"
          )
      }

      # group boxes sized to their text
      gr <- cc$groups
      if (!is.null(gr)) {
        for (i in seq_len(nrow(gr))) {
          tg <- grid::textGrob(gr$label[i], gp = grid::gpar(cex = lab_cex))
          grid::grid.roundrect(
            x = grid::unit(gr$x[i], "npc"), y = grid::unit(gr$y[i], "npc"),
            width = grid::grobWidth(tg) + grid::unit(3 * 2.5, "mm"),
            height = grid::grobHeight(tg) + grid::unit(2 * 2.5, "mm"),
            r = grid::unit(box_r, "mm"),
            gp = grid::gpar(fill = gr$fill[i], col = gr$color[i]),
            name = sprintf("seqplotr.node_link.group_box.%d", i)
          )
          grid::grid.text(
            gr$label[i], x = grid::unit(gr$x[i], "npc"),
            y = grid::unit(gr$y[i], "npc"),
            gp = grid::gpar(cex = lab_cex, col = "#1C1B1A"),
            name = sprintf("seqplotr.node_link.group_text.%d", i)
          )
        }
      }
      invisible()
    }
  )
)

#' Draw a node-link panel
#'
#' Draws nodes, directed or plain edges and labelled group boxes at given
#' positions. No layout is computed: `x` and `y` are used as supplied, in the
#' track's x and y scales.
#'
#' Columns of `nodes`: `id`, `x`, `y` (required), `bar_lo`, `bar_hi` (a
#' translucent vertical bar), `shape` (`circle`, `triangle_up`,
#' `triangle_down`, `diamond` or `square`; default `circle`), `fill`
#' (default `grey60`), `alpha` (default 1), `label`.
#'
#' Columns of `edges`: `from`, `to` (required node ids), `arrow` (default
#' `TRUE`), `color` (default `grey10`), `linetype` (default 1), `linewidth`
#' (default 1.4), `label`, `flag` (default `FALSE`; draws a marker at the edge
#' midpoint).
#'
#' Columns of `groups`: `x`, `y`, `label` (required), `fill`, `color`.
#'
#' @param nodes,edges,groups Data frames (or data tables) as described above.
#' @param aesthetics A `SeqAes`. Recognised: `node_mm` (3), `bar_linewidth`
#'   (3), `bar_alpha` (0.35), `arrow_mm` (2), `label_cex` (0.6),
#'   `edge_label_cex` (0.55), `box_r_mm` (1), `flag_color` (`"#D62728"`).
#' @param legend `NULL` for one key per node shape present, `FALSE` for none,
#'   or a `SeqLegendSpec` or list of `LegendKey` used as given.
#' @param legend_labels Named character vector relabelling the shape keys.
#' @param legend_extra List of `LegendKey` appended to the keys.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the generated legend.
#' @param ... Unused.
#' @return A `SeqNodeLinkR6` instance. Add it to a track with a `windows`
#'   range and fixed scales.
#' @examples
#' nodes <- data.frame(id = c("A", "B"), x = c(1, 2), y = c(0.2, 0.6))
#' edges <- data.frame(from = "A", to = "B", label = "fixed")
#' seq_node_link(nodes = nodes, edges = edges)
#' @export
seq_node_link <- function(nodes = NULL, edges = NULL, groups = NULL,
                          aesthetics = aes(), legend = NULL,
                          legend_labels = NULL, legend_extra = NULL,
                          show_legend = TRUE,
                          legend_position = "canvas_margin",
                          legend_side = "right", ...) {
  SeqNodeLinkR6$new(nodes = nodes, edges = edges, groups = groups,
                    aesthetics = aesthetics, legend = legend,
                    legend_labels = legend_labels, legend_extra = legend_extra,
                    show_legend = show_legend, legend_position = legend_position,
                    legend_side = legend_side, ...)
}
