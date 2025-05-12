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
    sc.data         = toy_counts,
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
    sc.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets,
    scale.factor    = ext_sf
  )
  expect_equal(norm, toy_enrich / ext_sf)       # exact division only
})

# --------------------------------------------------------------------------
test_that("chunked processing (groups) reproduces full result", {
  full <- performNormalization(
    sc.data         = toy_counts,
    enrichment.data = toy_enrich,
    gene.sets       = toy_sets,
    scale.factor    = rep(1, 4)
  )
  chunked <- performNormalization(
    sc.data         = toy_counts,
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
      sc.data         = toy_counts,
      enrichment.data = toy_enrich,
      gene.sets       = toy_sets,
      scale.factor    = c(1, 2)                 # wrong length
    ),
    "Length of 'scale.factor'"
  )
  
  # missing enrichment matrix
  expect_error(
    performNormalization(
      sc.data   = toy_counts,
      gene.sets = toy_sets
    ),
    "obtain enrichment matrix"
  )
  
  # gene-set names do not match enrichment cols
  bad_sets <- list(Other = c("g1", "g2"))
  expect_error(
    performNormalization(
      sc.data         = toy_counts,
      enrichment.data = toy_enrich,
      gene.sets       = bad_sets
    ),
    "None of the supplied gene sets match"
  )
})


