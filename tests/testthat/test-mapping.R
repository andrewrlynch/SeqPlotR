library(SeqPlotR)
library(GenomicRanges)
library(S4Vectors)

gr <- GRanges("chr1", IRanges(c(100, 200, 300), width = 50),
              score = c(1.0, 2.5, 0.8),
              AF    = c(0.1, 0.5, 0.9),
              label = c("a", "b", "c"))

test_that("map() captures unevaluated expressions", {
  m <- map(x = start, y = score)
  expect_s3_class(m, "SeqMap")
  expect_named(m, c("x", "y"))
})

test_that(".resolve_mapping resolves 'start' special", {
  m <- map(x = start)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$x, BiocGenerics::start(gr))
})

test_that(".resolve_mapping resolves 'end' special", {
  m <- map(x = end)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$x, BiocGenerics::end(gr))
})

test_that(".resolve_mapping resolves 'width' special", {
  m <- map(x = width)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$x, BiocGenerics::width(gr))
})

test_that(".resolve_mapping resolves 'mid' special", {
  m <- map(x = mid)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expected <- (BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2
  expect_equal(res$x, expected)
})

test_that(".resolve_mapping resolves bare column name", {
  m <- map(y = score)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$y, gr$score)
})

test_that(".resolve_mapping supports negation", {
  m <- map(y = -AF)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$y, -gr$AF)
})

test_that(".resolve_mapping supports function calls", {
  m <- map(y = log2(score + 1))
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expect_equal(res$y, log2(gr$score + 1))
})

test_that(".resolve_mapping supports arithmetic expressions", {
  m <- map(x = (start + end) / 2)
  res <- SeqPlotR:::.resolve_mapping(gr, m)
  expected <- (BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2
  expect_equal(res$x, expected)
})

test_that(".resolve_mapping returns list() when mapping is NULL", {
  res <- SeqPlotR:::.resolve_mapping(gr, NULL)
  expect_equal(res, list())
})

test_that(".resolve_mapping returns list() when data is NULL", {
  m <- map(x = start)
  res <- SeqPlotR:::.resolve_mapping(NULL, m)
  expect_equal(res, list())
})

test_that("aes() captures evaluated constants", {
  a <- aes(color = "blue", linewidth = 1.5)
  expect_s3_class(a, "SeqAes")
  expect_equal(a$color, "blue")
  expect_equal(a$linewidth, 1.5)
})
