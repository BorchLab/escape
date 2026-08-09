# test script for SpatialExperiment / SingleCellExperiment input -
# testcases are NOT comprehensive!
#
# Regression coverage for the bug reported in #180: escape.matrix(normalize =
# TRUE) threw an S4 dispatch error ("assay" on a NULL) for every
# SummarizedExperiment-derived input, not just SpatialExperiment.

gs <- toy_spe_sets()

test_that("escape.matrix(normalize = TRUE) works on a SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe()

  res <- escape.matrix(spe, gene.sets = gs, method = "UCell",
                       normalize = TRUE, min.size = NULL)

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(ncol(spe), length(gs)))
  expect_equal(colnames(res), names(gs))
  expect_equal(rownames(res), colnames(spe))
  expect_true(all(is.finite(res)))
})

test_that("escape.matrix(normalize = TRUE) works on a plain SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()

  res <- escape.matrix(sce, gene.sets = gs, method = "UCell",
                       normalize = TRUE, min.size = NULL)

  expect_equal(dim(res), c(ncol(sce), length(gs)))
  expect_equal(colnames(res), names(gs))
  expect_true(all(is.finite(res)))
})

test_that("normalized scores agree across SPE, SCE and raw matrix input", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe()
  sce <- make_toy_sce()
  mat <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  f <- function(x) escape.matrix(x, gene.sets = gs, method = "UCell",
                                 normalize = TRUE, min.size = NULL)

  expect_equal(f(spe), f(sce), tolerance = 1e-10)
  expect_equal(f(spe), f(mat), tolerance = 1e-10)
})

test_that("runEscape() preserves spatial metadata", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe()
  out <- runEscape(spe, gene.sets = gs, method = "UCell", min.size = NULL)

  expect_s4_class(out, "SpatialExperiment")
  expect_true("escape" %in% SingleCellExperiment::altExpNames(out))
  expect_equal(SpatialExperiment::spatialCoords(out),
               SpatialExperiment::spatialCoords(spe))
  expect_equal(out$sample_id, spe$sample_id)
  expect_equal(nrow(SpatialExperiment::imgData(out)),
               nrow(SpatialExperiment::imgData(spe)))

  # colnames must stay aligned between the altExp and its parent
  alt <- SingleCellExperiment::altExp(out, "escape")
  expect_equal(colnames(alt), colnames(spe))
  expect_equal(rownames(alt), names(gs))
})

test_that("multi-sample SpatialExperiment completes with chunking", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe(n.cells = 60, n.samples = 3)
  expect_gt(length(unique(spe$sample_id)), 1L)

  res <- escape.matrix(spe, gene.sets = gs, method = "UCell",
                       normalize = TRUE, min.size = NULL, groups = 25)
  expect_equal(dim(res), c(ncol(spe), length(gs)))
})

test_that("performNormalization() round-trips through an SCE altExp", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  obj <- runEscape(sce, gene.sets = gs, method = "UCell", min.size = NULL)

  out <- performNormalization(obj, assay = "escape", gene.sets = gs)

  expect_true("escape_normalized" %in% SingleCellExperiment::altExpNames(out))
  expect_equal(
    dim(SingleCellExperiment::altExp(out, "escape_normalized")),
    c(length(gs), ncol(sce))
  )
})

test_that("performPCA() works on a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  obj <- runEscape(sce, gene.sets = gs, method = "UCell", min.size = NULL)

  out <- performPCA(obj, assay = "escape")
  expect_s4_class(out, "SingleCellExperiment")
})

test_that("input.assay selects the expression matrix and errors informatively", {
  skip_if_not_installed("SpatialExperiment")
  spe <- make_toy_spe()

  res <- escape.matrix(spe, gene.sets = gs, method = "UCell",
                       min.size = NULL, input.assay = "logcounts")
  expect_equal(dim(res), c(ncol(spe), length(gs)))
  expect_equal(colnames(res), names(gs))

  expect_error(
    escape.matrix(spe, gene.sets = gs, method = "UCell",
                  min.size = NULL, input.assay = "nope"),
    "needs assay 'nope'"
  )
})

test_that(".resolve_input_assay maps names per object class", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()

  expect_equal(escape:::.resolve_input_assay(sce, "auto", "native"), "counts")
  expect_equal(escape:::.resolve_input_assay(sce, "auto", "plaid"), "logcounts")
  expect_equal(escape:::.resolve_input_assay(SeuratObject::pbmc_small,
                                             "logcounts", "native"), "data")
  expect_equal(escape:::.resolve_input_assay(SeuratObject::pbmc_small,
                                             "auto", "native"), "counts")

  # an SCE without logcounts must say how to get them
  bare <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = SummarizedExperiment::assay(sce, "counts"))
  )
  expect_error(escape:::.resolve_input_assay(bare, "auto", "plaid"),
               "logNormCounts")
})
