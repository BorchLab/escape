# test script for plaid.R - testcases are NOT comprehensive!
#
# The suite runs at testthat edition 2 (no Config/testthat/edition in
# DESCRIPTION), so local_mocked_bindings() is unavailable and assignInNamespace()
# would fail R CMD check. .plaid_dispatch() therefore takes an injectable `fn`,
# which lets the dispatch logic be tested on CI whether or not plaid is
# installed. Everything needing the real package is skipped explicitly.

# genes chosen to exist in pbmc_small so UCell does not warn about imputation
mini_gs <- list(B = c("MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A"),
                T = c("CD3E", "CD3D", "CD7", "CD8A", "IL7R"))

# A stand-in for plaid's replaid.* family: honors min.genes/max.genes, returns
# gene sets x cells like the real thing, and records what it was called with.
fake_plaid <- function(X, matG, min.genes = 5, max.genes = 500, alpha = 0, ...) {
  keep <- names(matG)[lengths(matG) >= min.genes & lengths(matG) <= max.genes]
  m <- matrix(seq_len(length(keep) * ncol(X)),
              nrow = length(keep), ncol = ncol(X),
              dimnames = list(keep, colnames(X)))
  attr(m, "call.args") <- list(min.genes = min.genes, max.genes = max.genes,
                               alpha = alpha)
  m
}

big_sets <- list(Big   = paste0("G", 1:900),
                 Small = paste0("G", 1:10))
big_X <- matrix(1, nrow = 1000, ncol = 4,
                dimnames = list(paste0("G", 1:1000), paste0("c", 1:4)))

# --------------------------------------------------------------------------
# Dispatch logic - no plaid needed
# --------------------------------------------------------------------------
test_that("the max.genes cap is lifted above the largest gene set", {
  out <- escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                                  min.size = 5, fn = fake_plaid)
  # plaid's own default (500) would have silently dropped the 900-gene set
  expect_setequal(colnames(out), names(big_sets))
  expect_equal(dim(out), c(4L, 2L))
})

test_that("output is cells x gene sets with the input dimnames", {
  out <- escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                                  min.size = 5, fn = fake_plaid)
  expect_equal(rownames(out), colnames(big_X))
  expect_equal(colnames(out), names(big_sets))
})

test_that("min.size is translated to plaid's min.genes", {
  captured <- NULL
  spy <- function(X, matG, min.genes = 5, max.genes = 500, ...) {
    captured <<- list(min.genes = min.genes, max.genes = max.genes)
    fake_plaid(X, matG, min.genes = min.genes, max.genes = max.genes)
  }
  escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA", min.size = 8, fn = spy)
  expect_equal(captured$min.genes, 8L)
  expect_gte(captured$max.genes, 900L)
})

test_that("backend.args reach the plaid function and typos are rejected", {
  captured <- NULL
  spy <- function(X, matG, min.genes = 5, max.genes = 500, alpha = 0, ...) {
    captured <<- alpha
    fake_plaid(X, matG, min.genes = min.genes, max.genes = max.genes)
  }
  escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                           backend.args = list(alpha = 0.25), fn = spy)
  expect_equal(captured, 0.25)

  expect_error(
    escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                             backend.args = list(alhpa = 1), fn = fake_plaid),
    "Unknown `backend.args`"
  )
})

test_that("dropped gene sets warn rather than silently shrinking the output", {
  # a user-supplied cap that excludes the 900-gene set
  expect_warning(
    out <- escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                                    backend.args = list(max.genes = 100),
                                    fn = fake_plaid),
    "plaid dropped 1 of 2 gene set"
  )
  expect_equal(colnames(out), "Small")

  # nothing scored at all is an error, not an empty matrix
  expect_error(
    escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                             backend.args = list(max.genes = 1),
                             fn = fake_plaid),
    "scored none of the"
  )
})

test_that("a plaid result with the wrong number of cells is rejected", {
  wrong <- function(X, matG, ...) {
    m <- matrix(1, nrow = length(matG), ncol = 2,
                dimnames = list(names(matG), c("c1", "c2")))
    m
  }
  expect_error(
    escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA", fn = wrong),
    "returned 2 cell"
  )
})

# --------------------------------------------------------------------------
# Method / backend resolution
# --------------------------------------------------------------------------
test_that(".resolve_backend normalizes method names and backends", {
  expect_equal(escape:::.resolve_backend("ssgsea", "plaid"),
               list(method = "SSGSEA", backend = "plaid"))
  expect_equal(escape:::.resolve_backend("AUCell", "native"),
               list(method = "AUCELL", backend = "native"))
  expect_error(escape:::.resolve_backend("nope", "native"), "Unknown `method`")
})

test_that("plaid-only methods route to plaid regardless of backend", {
  for (m in c("PLAID", "singscore", "scSE")) {
    expect_message(res <- escape:::.resolve_backend(m, "native"),
                   "provided by the plaid backend")
    expect_equal(res$backend, "plaid")
  }
  # asking for plaid explicitly is silent
  expect_silent(escape:::.resolve_backend("scSE", "plaid"))
})

test_that("escape.matrix rejects a plaid-only method under the native engine only via routing", {
  expect_error(escape.matrix(SeuratObject::pbmc_small, mini_gs, method = "nope"),
               "Unknown `method`")
})

# --------------------------------------------------------------------------
# Guards
# --------------------------------------------------------------------------
test_that("the plaid backend errors informatively when plaid is missing", {
  skip_if(requireNamespace("plaid", quietly = TRUE))
  expect_error(
    escape.matrix(SeuratObject::pbmc_small, mini_gs, min.size = 0,
                  backend = "plaid"),
    "plaid not installed"
  )
})

test_that("dots are refused on the plaid path", {
  expect_error(
    escape.matrix(SeuratObject::pbmc_small, mini_gs, min.size = 0,
                  backend = "plaid", maxRank = 100),
    "not forwarded when backend"
  )
})

test_that("backend.args must be a fully named list", {
  expect_error(
    escape:::.plaid_dispatch(big_X, big_sets, "ssGSEA",
                             backend.args = list(1), fn = fake_plaid),
    "fully named list"
  )
})

# --------------------------------------------------------------------------
# Provenance and native-path regression
# --------------------------------------------------------------------------
test_that("scores are stamped with the engine that produced them", {
  res <- escape.matrix(SeuratObject::pbmc_small, mini_gs, method = "UCell",
                       min.size = 0)
  prov <- attr(res, "escape.backend")
  expect_equal(prov$backend, "native")
  expect_equal(prov$method, "UCell")
})

test_that("runEscape records provenance on the object", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_toy_sce()
  out <- runEscape(sce, gene.sets = toy_spe_sets(), method = "UCell",
                   min.size = NULL)
  prov <- S4Vectors::metadata(out)[["escape_backend"]]
  expect_equal(prov$backend, "native")

  obj <- runEscape(SeuratObject::pbmc_small, gene.sets = mini_gs,
                   method = "UCell", min.size = 0)
  expect_equal(SeuratObject::Misc(obj, "escape_backend")$backend, "native")
})

test_that("the native path is unchanged by the new arguments", {
  a <- escape.matrix(SeuratObject::pbmc_small, mini_gs, min.size = 0)
  b <- escape.matrix(SeuratObject::pbmc_small, mini_gs, min.size = 0,
                     backend = "native", input.assay = "auto")
  expect_equal(a, b)
})

# --------------------------------------------------------------------------
# Integration - needs the real package
# --------------------------------------------------------------------------
test_that("plaid backend reproduces escape's orientation and dimnames", {
  skip_if_not_installed("plaid")
  pbmc <- SeuratObject::pbmc_small

  for (m in c("ssGSEA", "GSVA", "UCell", "AUCell")) {
    res <- escape.matrix(pbmc, mini_gs, method = m, min.size = 0,
                         backend = "plaid")
    expect_equal(dim(res), c(ncol(pbmc), length(mini_gs)),
                 info = paste("method:", m))
    expect_equal(colnames(res), names(mini_gs), info = paste("method:", m))
    expect_equal(rownames(res), colnames(pbmc), info = paste("method:", m))
    expect_true(all(is.finite(res)), info = paste("method:", m))
  }
})

test_that("plaid-only methods return sane matrices", {
  skip_if_not_installed("plaid")
  pbmc <- SeuratObject::pbmc_small

  for (m in c("PLAID", "singscore", "scSE")) {
    res <- escape.matrix(pbmc, mini_gs, method = m, min.size = 0)
    expect_equal(dim(res), c(ncol(pbmc), length(mini_gs)),
                 info = paste("method:", m))
    expect_true(all(is.finite(res)), info = paste("method:", m))
  }
})

test_that("no gene sets are dropped at Hallmark scale under defaults", {
  skip_if_not_installed("plaid")
  pbmc <- SeuratObject::pbmc_small
  feats <- rownames(pbmc)
  # a set larger than plaid's default max.genes = 500
  wide <- list(Wide = rep(feats, length.out = 600), Narrow = feats[1:10])

  # the drop warning must not fire - that is what "nothing dropped" means here
  expect_no_warning(
    res <- suppressMessages(
      escape.matrix(pbmc, wide, method = "ssGSEA", min.size = 0,
                    backend = "plaid")
    )
  )
  expect_equal(colnames(res), names(wide))
})

test_that("plaid backend is broadly concordant with the native engines", {
  skip_if_not_installed("plaid")
  skip_if_not_installed("UCell")
  pbmc <- SeuratObject::pbmc_small

  nat <- escape.matrix(pbmc, mini_gs, method = "UCell", min.size = 0)
  pla <- escape.matrix(pbmc, mini_gs, method = "UCell", min.size = 0,
                       backend = "plaid")

  # Deliberately loose. plaid documents replaid.ucell as near-identical to
  # UCell, but measured pooled Pearson is ~0.79 on pbmc_small (230 genes) and
  # ~0.83 on a 2000-gene simulation - the scores are correlated but on a
  # different scale. This guards against the backend breaking outright, not
  # against drift, and must not be tightened into an equivalence claim the
  # package cannot make.
  expect_gt(stats::cor(as.vector(nat), as.vector(pla)), 0.6)
})

test_that("backend.args change the scores", {
  skip_if_not_installed("plaid")
  pbmc <- SeuratObject::pbmc_small
  a <- escape.matrix(pbmc, mini_gs, method = "ssGSEA", min.size = 0,
                     backend = "plaid")
  b <- escape.matrix(pbmc, mini_gs, method = "ssGSEA", min.size = 0,
                     backend = "plaid", backend.args = list(alpha = 0.25))
  expect_false(isTRUE(all.equal(unname(a), unname(b))))
})

test_that("plaid scores feed the downstream workflow", {
  skip_if_not_installed("plaid")
  pbmc <- SeuratObject::pbmc_small
  obj <- runEscape(pbmc, gene.sets = mini_gs, method = "UCell", min.size = 0,
                   backend = "plaid")

  expect_true("escape" %in% SeuratObject::Assays(obj))
  expect_equal(SeuratObject::Misc(obj, "escape_backend")$backend, "plaid")
  expect_no_error(performNormalization(obj, assay = "escape",
                                       gene.sets = mini_gs))
  expect_no_error(performPCA(obj, assay = "escape"))
})
