library(SeqPlotR)
library(GenomicRanges)

# ── Positional layout ──────────────────────────────────────────────────────────

test_that("first track always placed top-left regardless of direction", {
  p <- seq_plot() %+% seq_track(direction = "under", track_id = "A")
  expect_length(p$rows, 1)
  expect_length(p$rows[[1]], 1)
})

test_that("%+% with direction='right' appends to current row", {
  p <- seq_plot() %+%
    seq_track(track_id = "A") %+%
    seq_track(track_id = "B", direction = "right")
  expect_length(p$rows, 1)
  expect_length(p$rows[[1]], 2)
})

test_that("%+% with direction='under' starts a new row", {
  p <- seq_plot() %+%
    seq_track(track_id = "A") %+%
    seq_track(track_id = "B", direction = "under")
  expect_length(p$rows, 2)
  expect_length(p$rows[[2]], 1)
})

test_that("three tracks: two right then one under", {
  p <- seq_plot() %+%
    seq_track(track_id = "A") %+%
    seq_track(track_id = "B", direction = "right") %+%
    seq_track(track_id = "C", direction = "under")
  expect_length(p$rows, 2)
  expect_length(p$rows[[1]], 2)
  expect_length(p$rows[[2]], 1)
})

# ── Convenience aliases ────────────────────────────────────────────────────────

test_that("%|% is equivalent to direction='right'", {
  p1 <- seq_plot() %+% seq_track(track_id="A") %+% seq_track(track_id="B", direction="right")
  p2 <- seq_plot() %|% seq_track(track_id="A") %|% seq_track(track_id="B")
  expect_equal(length(p1$rows[[1]]), length(p2$rows[[1]]))
})

test_that("%__% is equivalent to direction='under'", {
  p1 <- seq_plot() %+% seq_track(track_id="A") %+% seq_track(track_id="B", direction="under")
  p2 <- seq_plot() %|% seq_track(track_id="A") %__% seq_track(track_id="B")
  expect_length(p1$rows, 2)
  expect_length(p2$rows, 2)
})

# ── Patchwork mode ─────────────────────────────────────────────────────────────

test_that("patchwork mode appends to flat tracks list, direction ignored", {
  p <- seq_plot(layout = "AB") %+%
    seq_track(track_id = "A", direction = "under") %+%
    seq_track(track_id = "B", direction = "under")
  expect_length(p$tracks, 2)
  expect_null(p$rows)  # rows not used in patchwork mode
})

# ── Plot-level features ────────────────────────────────────────────────────────

test_that("SeqAnnotation goes to plot_annotations", {
  ann <- structure(list(), class = "SeqAnnotation")
  p <- seq_plot() %+% seq_track(track_id="A") %+% ann
  expect_length(p$plot_annotations, 1)
})

test_that("adding element with no tracks errors", {
  elem <- structure(list(), class = c("SeqPoint", "SeqElement"))
  expect_error(seq_plot() %+% elem, "No tracks")
})

# ── Track existence validation for SeqLink ────────────────────────────────────

test_that("plot-level link errors if t1 track not yet defined", {
  lnk <- structure(list(t0 = "A", t1 = "B"), class = c("SeqString", "SeqLink"))
  expect_error(
    seq_plot() %+% seq_track(track_id = "A") %+% lnk,
    "B"
  )
})

test_that("plot-level link succeeds when both tracks are defined", {
  lnk <- structure(list(t0 = "A", t1 = "B"), class = c("SeqString", "SeqLink"))
  p <- seq_plot() %+%
    seq_track(track_id = "A") %+%
    seq_track(track_id = "B", direction = "under") %+%
    lnk
  expect_length(p$plot_links, 1)
})

# ── Within-track link locking ─────────────────────────────────────────────────

test_that("%+% on SeqTrack locks SeqLink t0/t1 to parent track_id", {
  lnk <- structure(list(t0 = "X", t1 = "Y"), class = c("SeqArc", "SeqLink"))
  trk <- seq_track(track_id = "A")
  trk %+% lnk
  # S3 stub is copy-semantics; inspect the copy stored on the track instead
  # of the outer `lnk`. When SeqLink becomes an R6 class, both will agree.
  expect_equal(trk$elements[[1]]$t0, "A")
  expect_equal(trk$elements[[1]]$t1, "A")
})
