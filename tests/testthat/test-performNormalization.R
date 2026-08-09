# test script for performNormalization.R - testcases are NOT comprehensive!


# --------------------------------------------------------------------------
# helper: tiny toy dataset --------------------------------------------------
toy_counts <- Matrix::sparseMatrix(
  i = c(1, 3, 2, 1, 3),          # g1 g3 g2 g1 g3
  j = c(1, 1, 2, 3, 4),          # c1 c1 c2 c3 c4
  x = c(5, 2, 3, 4, 1),
  dims = c(3, 4),
  dimnames = list(c("g1", "g2", "g3"), paste0("c", 1:4))
)

toy_enrich <- matrix(
  c(3, 6, 4, 8,   # Set1
    2, 4, 3, 6),  # Set2
  nrow = 4,
  dimnames = list(paste0("c", 1:4), c("Set1", "Set2"))
)

toy_sets <- list(
  Set1 = c("g1", "g2"),
  Set2 = c("g2", "g3")
)

# --------------------------------------------------------------------------
test_that("matrix input: internal scale factors + log transform", {
  norm <- performNormalization(
    input.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets
  )
  
  # dimensions and finite values
  expect_equal(dim(norm), dim(toy_enrich))
  expect_true(all(is.finite(norm)))
  expect_false(anyNA(norm))
  
  # manual check on first cell / gene-set
  gs_counts_c1 <- c(
    Set1 = sum(toy_counts[c("g1", "g2"), "c1"] != 0),
    Set2 = sum(toy_counts[c("g2", "g3"), "c1"] != 0)
  )
  manual <- log1p(toy_enrich["c1", ] / gs_counts_c1 + 1e-6)
  expect_equal(unname(norm["c1", ]), unname(manual))
})

# --------------------------------------------------------------------------
test_that("matrix input: external scale.factor bypasses log step", {
  ext_sf <- c(2, 2, 2, 2)                       # one per cell
  norm <- performNormalization(
    input.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets,
    scale.factor    = ext_sf
  )
  expect_equal(norm, toy_enrich / ext_sf)       # exact division only
})

# --------------------------------------------------------------------------
test_that("chunked processing (groups) reproduces full result", {
  full <- performNormalization(
    input.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets,
    scale.factor    = rep(1, 4)
  )
  chunked <- performNormalization(
    input.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets,
    scale.factor    = rep(1, 4),
    groups          = 2                         # split into two chunks
  )
  expect_equal(full, chunked)
})

# --------------------------------------------------------------------------
test_that("error handling works", {
  # scale.factor length mismatch
  expect_error(
    performNormalization(
      input.data         = toy_counts,
      enrichment.data = toy_enrich,
      gene.sets       = toy_sets,
      scale.factor    = c(1, 2)                 # wrong length
    ),
    "Length of 'scale.factor'"
  )
  
  # missing enrichment matrix
  expect_error(
    performNormalization(
      input.data   = toy_counts,
      gene.sets = toy_sets
    ),
    "obtain enrichment matrix"
  )
  
  # gene-set names do not match enrichment cols
  bad_sets <- list(Other = c("g1", "g2"))
  expect_error(
    performNormalization(
      input.data         = toy_counts,
      enrichment.data = toy_enrich,
      gene.sets       = bad_sets
    ),
    "None of the supplied gene sets match"
  )

  # a gene set for only some of the columns must not silently misalign
  expect_error(
    performNormalization(
      input.data      = toy_counts,
      enrichment.data = toy_enrich,
      gene.sets       = list(Set1 = c("g1", "g2"))
    ),
    "No gene set supplied for enrichment column"
  )
})

# --------------------------------------------------------------------------
# Underscored set names (HALLMARK_*, GO_*, REACTOME_*) used to be mangled to
# hyphens unconditionally, which dropped every such set for non-Seurat input.
under_sets <- list(HALLMARK_SET_ONE = c("g1", "g2"),
                   Set2             = c("g2", "g3"))
under_enrich <- toy_enrich
colnames(under_enrich) <- names(under_sets)

test_that("underscored gene-set names normalize on matrix input", {
  norm <- performNormalization(
    input.data      = toy_counts,
    enrichment.data = under_enrich,
    gene.sets       = under_sets
  )
  expect_equal(dim(norm), dim(under_enrich))
  expect_equal(colnames(norm), names(under_sets))
  expect_true(all(is.finite(norm)))

  gs_counts_c1 <- c(
    sum(toy_counts[c("g1", "g2"), "c1"] != 0),
    sum(toy_counts[c("g2", "g3"), "c1"] != 0)
  )
  manual <- log1p(under_enrich["c1", ] / gs_counts_c1 + 1e-6)
  expect_equal(unname(norm["c1", ]), unname(manual))
})

test_that("underscored names give identical scores across input classes", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  gs  <- toy_spe_sets()   # HALLMARK_SET_A + SetB
  mat <- as.matrix(SummarizedExperiment::assay(sce, "counts"))

  f <- function(x) escape.matrix(x, gene.sets = gs, method = "UCell",
                                 normalize = TRUE, min.size = NULL)

  from_sce <- f(sce)
  from_mat <- f(mat)
  expect_equal(from_sce, from_mat, tolerance = 1e-10)
  expect_equal(colnames(from_sce), names(gs))
})

# --------------------------------------------------------------------------
test_that("supplied enrichment.data wins over scores held on the object", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  gs  <- toy_spe_sets()
  obj <- runEscape(sce, gene.sets = gs, method = "UCell", min.size = NULL)

  # a matrix that is deliberately nothing like the stored scores
  fake <- matrix(1, nrow = ncol(sce), ncol = length(gs),
                 dimnames = list(colnames(sce), names(gs)))

  expect_warning(
    out <- performNormalization(obj, enrichment.data = fake,
                                assay = "escape", gene.sets = gs),
    "using `enrichment.data`"
  )

  from_fake <- performNormalization(SummarizedExperiment::assay(sce, "counts"),
                                    enrichment.data = fake, gene.sets = gs)
  expect_equal(
    Matrix::t(SummarizedExperiment::assay(
      SingleCellExperiment::altExp(out, "escape_normalized"))),
    from_fake,
    tolerance = 1e-10, ignore_attr = TRUE
  )
})

test_that("a nonexistent enrichment assay names what is available", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  gs  <- toy_spe_sets()
  obj <- runEscape(sce, gene.sets = gs, method = "UCell", min.size = NULL)

  expect_error(
    performNormalization(obj, assay = "not_there", gene.sets = gs),
    "Could not find enrichment assay 'not_there'"
  )
  expect_error(
    performNormalization(obj, assay = "not_there", gene.sets = gs),
    "Available: escape"
  )
})


