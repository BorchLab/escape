# test script for performPCA.R - testcases are NOT comprehensive!

# -------------------------------------------------------------------------
# 1.  Matrix input utilities ------------------------------------------------
# -------------------------------------------------------------------------
set.seed(123)
mat_small  <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20,
                     dimnames = list(paste0("cell", 1:100),
                                     paste0("set",  1:20)))

test_that("Matrix input returns well-formed list", {
  pca_res <- performPCA(mat_small, scale = FALSE, n.dim = 5)
  
  expect_type(pca_res, "list")
  expect_named(pca_res,
               c("PCA", "eigen_values", "contribution", "rotation"),
               ignore.order = TRUE)
  expect_equal(dim(pca_res$PCA), c(100, 5))          # 100 cells × 5 PCs
  expect_length(pca_res$eigen_values, 20)
  expect_length(pca_res$contribution, 20)
  expect_equal(dim(pca_res$rotation), c(20, 5))       # gene sets × loadings
})

test_that("Scaling alters the embeddings", {
  pca_unscaled <- performPCA(mat_small, scale = FALSE, n.dim = 5)$PCA
  pca_scaled   <- performPCA(mat_small, scale = TRUE,  n.dim = 5)$PCA
  expect_false(isTRUE(all.equal(pca_unscaled, pca_scaled)))
})

test_that("n.dim supplied as a vector is honoured", {
  pca_res <- performPCA(mat_small, n.dim = 1:7)
  expect_equal(ncol(pca_res$PCA), 7)
})

# -------------------------------------------------------------------------
# 2.  Seurat workflow -------------------------------------------------------
# -------------------------------------------------------------------------
if (requireNamespace("SeuratObject", quietly = TRUE) &&
    requireNamespace("Seurat",        quietly = TRUE)) {
  
  test_that("Seurat object gains a DimReduc slot", {
    pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")  # helper fixture
    pbmc_small <- performPCA(pbmc_small, assay = "escape",
                            n.dim = 6, reduction.name = "escPCA")
    
    expect_s4_class(pbmc_small[["escPCA"]], "DimReduc")
    emb <- SeuratObject::Embeddings(pbmc_small[["escPCA"]])
    expect_equal(dim(emb)[2], 2)                           
  })
}

# -------------------------------------------------------------------------
# 3.  Error handling --------------------------------------------------------
# -------------------------------------------------------------------------
test_that("performPCA() fails on invalid input types", {
  expect_error(performPCA("not a matrix"),
               "must be a matrix/data.frame or a Seurat/SCE object")
})

test_that("performPCA() fails on non-numeric matrix", {
  bad_mat <- matrix(letters[1:20], nrow = 4)
  expect_error(performPCA(bad_mat), "Enrichment matrix must be numeric")
})
