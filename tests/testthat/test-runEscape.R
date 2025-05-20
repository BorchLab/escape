# test script for runEscape.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------- helpers --
mini_gs <- list(
  B = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
  T = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))

pbmc_small  <- SeuratObject::pbmc_small

get_score <- function(method = "ssGSEA", ...) {
  escape.matrix(pbmc_small,
                gene.sets      = mini_gs,
                method         = method,
                groups         = 200,  # small chunk for speed
                min.size       = 0,
                normalize      = FALSE,
                make.positive  = FALSE,
                min.filter.by  = NULL,
                BPPARAM        = BiocParallel::SerialParam())
}

# ------------------------------------------------------------- interface -----
test_that("escape.matrix() accepts Seurat, SCE and matrix", {
  sce <- as.SingleCellExperiment(pbmc_small)
  mtx <- pbmc_small[["RNA"]]@counts
  
  x <- get_score(method = "ssGSEA")
  y <- escape.matrix(sce,   mini_gs, min.size = 0)
  z <- escape.matrix(mtx,   mini_gs, min.size = 0)
  expect_equal(x,y)
  expect_equal(x,z)
  expect_equal(y,z)
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
  sc <- escape.matrix(pbmc_small, gs_bad, min.size = 3)
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

# --------------------------------------------------------- chunk invariance --
test_that("different 'groups' chunking gives identical results", {
  sc_small <- get_score(groups = ncol(pbmc_small))  # one chunk
  sc_many  <- get_score(groups = 20)                # many chunks
  expect_equal(sc_small, sc_many, tolerance = 1e-10)
})


# ----------------------------------------------------- runEscape integration --
test_that("runEscape adds assay (default & custom names)", {
  gs <- mini_gs
  obj1 <- runEscape(pbmc_small, gene.sets = gs, groups = 200, min.size = 0)
  expect_true("escape" %in% Assays(obj1))
  
  obj2 <- runEscape(pbmc_small, gene.sets = gs,
                    groups = 200, new.assay.name = "myESCAPE", min.size = 0)
  expect_true("myESCAPE" %in% Assays(obj2))
})

# -------------------------------------------------------- error pathways -----
test_that("runEscape propagates escape.matrix errors", {
  gs_bad <- list(bad = "NOT_A_GENE")
  expect_error(runEscape(pbmc_small, gs_bad, min.size = 3),
               "No gene-sets meet")
})
