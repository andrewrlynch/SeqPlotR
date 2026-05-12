library(testthat)
library(GenomicRanges)

gr_genes <- GRanges("chr1",
  IRanges(start = c(100,  400, 1000, 1500),
          end   = c(300,  600, 1200, 1700)),
  gene_id   = c("geneA", "geneA", "geneB", "geneB"),
  gene_name = c("GeneA", "GeneA", "GeneB", "GeneB"),
  strand_col = c("+", "+", "-", "-"),
  feat_type  = c("exon", "exon", "exon", "exon")
)
win <- GRanges("chr1", IRanges(1, 2000))

make_layout <- function() {
  trk <- seq_track(data = gr_genes, windows = win)
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  sp$layout$panelBounds[[1]]
}

test_that("default label (no gene.label aes) renders without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  el <- seq_gene(map(group = gene_id, strand = strand_col, type = feat_type))
  el$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("gene.label position = 'top' renders without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  el <- seq_gene(map(group = gene_id, strand = strand_col, type = feat_type),
                 aesthetics = aes(gene.label = aes(position = "top")))
  el$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("gene.label position = c('start','top') renders without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  el <- seq_gene(map(group = gene_id, strand = strand_col, type = feat_type),
                 aesthetics = aes(gene.label = aes(position = c("start","top"),
                                                    color = "navy",
                                                    size = 0.5)))
  el$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("gene.label position = c('end','bottom') renders without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  el <- seq_gene(map(group = gene_id, strand = strand_col, type = feat_type),
                 aesthetics = aes(gene.label = aes(position = c("end","bottom"))))
  el$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})
