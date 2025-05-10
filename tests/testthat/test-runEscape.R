# test script for runEscape.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------- helpers --
mini_gs <- list(
  B = c("MS4A1", "CD79B", "CD79A"),
  T = c("CD3E", "CD3D", "CD3G")
)

get_score <- function(method = "ssGSEA", ...) {
  escape.matrix(pbmc_small,
                gene.sets      = mini_gs,
                method         = method,
                groups         = 200,  # small chunk for speed
                min.size       = 3,
                normalize      = FALSE,
                make.positive  = FALSE,
                min.expr.cells = 0,
                min.filter.by  = NULL,
                BPPARAM        = BiocParallel::SerialParam(),
                ...)
}

# ------------------------------------------------------------- interface -----
test_that("escape.matrix() accepts Seurat, SCE and matrix", {
  sce <- as.SingleCellExperiment(pbmc_small)
  mtx <- pbmc_small[["RNA"]]@counts
  
  expect_silent(get_score(method = "ssGSEA"))
  expect_silent(escape.matrix(sce,   mini_gs))
  expect_silent(escape.matrix(mtx,   mini_gs))
})

test_that("invalid method triggers error", {
  expect_error(get_score(method = "foobar"),
               "must be one of")
})

# ---------------------------------------------------------- output shape -----
test_that("output matrix has cells × gene-sets and ordered columns", {
  sc <- get_score()
  expect_equal(dim(sc), c(ncol(pbmc_small), length(mini_gs)))
  expect_equal(colnames(sc), names(mini_gs))
  expect_true(setequal(rownames(sc), colnames(pbmc_small)))
})

# ------------------------------------------------------- min.size filter -----
test_that("gene-sets failing min.size are dropped with message", {
  gs_bad <- c(mini_gs, Junk = "ZZZ_UNKNOWN_GENE")
  expect_message(
    sc <- escape.matrix(pbmc_small, gs_bad, min.size = 3),
    "No.*ZZZ_UNKNOWN_GENE"
  )
  expect_false("Junk" %in% colnames(sc))
})

# --------------------------------------------------- min.expr.cells (global) -
test_that("min.expr.cells filters genes globally", {
  sc0 <- get_score(min.expr.cells = 0)
  sc5 <- get_score(min.expr.cells = 0.5)   # keep genes in ≥50% of cells
  expect_true(is.matrix(sc5) && is.matrix(sc0))
  # dimension equality (gene filter should not affect cell × set shape)
  expect_equal(dim(sc0), dim(sc5))
})

# ------------------------------------------ min.expr.cells with min.filter.by -
test_that("per-group gene filter behaves and is cluster-specific", {
  # Use seurat_clusters as grouping; expect same shape but different values
  sc_global <- get_score(min.expr.cells = 0.2)
  sc_group  <- get_score(min.expr.cells = 0.2,
                         min.filter.by  = "seurat_clusters")
  expect_equal(dim(sc_global), dim(sc_group))
  expect_false(isTRUE(all.equal(sc_global, sc_group)))
})

# --------------------------------------------------------- chunk invariance --
test_that("different 'groups' chunking gives identical results", {
  sc_small <- get_score(groups = ncol(pbmc_small))  # one chunk
  sc_many  <- get_score(groups = 20)                # many chunks
  expect_equal(sc_small, sc_many, tolerance = 1e-10)
})

# ---------------------------------------------------- normalise / positive ---
test_that("normalisation and make.positive shift range correctly", {
  norm <- get_score(normalize = TRUE, make.positive = TRUE)
  expect_true(all(norm >= 0))
})

# ---------------------------------------------------------- back-end tests ---
backends <- c("ssGSEA", "GSVA", "UCell", "AUCell")
for (m in backends) {
  test_that(paste0("method = '", m, "' runs if backend present"), {
    pkg <- switch(m,
                  GSVA  = "GSVA",
                  UCell = "UCell",
                  AUCell= "AUCell",
                  ssGSEA= NA)
    skip_if(!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE),
            paste("skip:", pkg, "not installed"))
    expect_silent(get_score(method = m))
  })
}

# ----------------------------------------------------- runEscape integration --
test_that("runEscape adds assay (default & custom names)", {
  gs <- mini_gs
  obj1 <- runEscape(pbmc_small, gene.sets = gs, groups = 200)
  expect_true("escape" %in% Assays(obj1))
  
  obj2 <- runEscape(pbmc_small, gene.sets = gs,
                    groups = 200, new.assay.name = "myESCAPE")
  expect_true("myESCAPE" %in% Assays(obj2))
})

# -------------------------------------------------------- error pathways -----
test_that("runEscape propagates escape.matrix errors", {
  gs_bad <- list(bad = "NOT_A_GENE")
  expect_error(runEscape(pbmc_small, gs_bad, min.size = 3),
               "No gene-sets meet")
})
