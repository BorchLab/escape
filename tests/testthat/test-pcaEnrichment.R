# test script for pcaEnrichment.R - testcases are NOT comprehensive!

pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")

# PCA (small data → very fast)
pbmc_small <- escape::performPCA(pbmc_small, assay = "escape")

# Convenience: pull the raw list returned by .grabDimRed()
pca_list <- escape:::.grabDimRed(pbmc_small, "escape.PCA")


## -----------------------------------------------------------------
## 1.  Basic behaviour ---------------------------------------------
## -----------------------------------------------------------------
test_that("returns a ggplot object for Seurat input", {
  g <- escape::pcaEnrichment(pbmc_small,
                             dimRed = "escape.PCA",
                             x.axis = "PC1",
                             y.axis = "PC2")
  expect_s3_class(g, "gg")
  expect_true(ggplot2::is_ggplot(g))
})

test_that("returns a ggplot object when supplied the raw PCA list", {
  g <- escape::pcaEnrichment(pca_list,
                             x.axis = "PC1",
                             y.axis = "PC2")
  expect_s3_class(g, "gg")
})

## -----------------------------------------------------------------
## 2.  Axis-label handling -----------------------------------------
## -----------------------------------------------------------------
test_that("percentage labels are appended when requested", {
  g <- escape::pcaEnrichment(pbmc_small,
                             dimRed  = "escape.PCA",
                             x.axis  = "PC1",
                             y.axis  = "PC2",
                             add.percent.contribution = TRUE)
  expect_match(g$labels$x, "PC1.*%")
  expect_match(g$labels$y, "PC2.*%")
})

## -----------------------------------------------------------------
## 3.  Faceting -----------------------------------------------------
## -----------------------------------------------------------------
test_that("faceting works and errors appropriately", {
  g <- escape::pcaEnrichment(pbmc_small,
                             dimRed   = "escape.PCA",
                             facet.by = "groups")
  expect_true("FacetGrid" %in% class(g$facet))
  
  # facet.by with raw list → error
  expect_error(
    escape::pcaEnrichment(pca_list, facet.by = "groups"),
    "input.data' must be a Seurat / SCE object or the list from performPCA().",
    fixed = TRUE
  )
  
  # invalid facet.by column
  expect_error(
    escape::pcaEnrichment(pbmc_small,
                          dimRed   = "escape.PCA",
                          facet.by = "not_a_col"),
    "Please select a variable in your meta data to use for facet.by.",
    fixed = TRUE
  )
})

## -----------------------------------------------------------------
## 4.  Plot styles --------------------------------------------------
## -----------------------------------------------------------------
test_that("`style = 'hex'` produces a `GeomHex` layer (when hexbin present)", {
  skip_if_not_installed("hexbin")
  g <- escape::pcaEnrichment(pbmc_small,
                             dimRed = "escape.PCA",
                             style  = "hex")
  geoms <- vapply(g$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomHex" %in% geoms)
})

## -----------------------------------------------------------------
## 5.  Biplot overlay ----------------------------------------------
## -----------------------------------------------------------------
test_that("display.factors adds segment & text layers", {
  g <- escape::pcaEnrichment(pbmc_small,
                             dimRed            = "escape.PCA",
                             display.factors   = TRUE,
                             number.of.factors = 5)
  geoms <- vapply(g$layers, function(x) class(x$geom)[1], character(1))
  expect_true(any(c("GeomSegment", "GeomLabel") %in% geoms))
})

## -----------------------------------------------------------------
## 6.  Error handling for bad inputs -------------------------------
## -----------------------------------------------------------------
test_that("bad inputs are rejected with informative errors", {
  expect_error(
    escape::pcaEnrichment(mtcars),
    "input.data' must be a Seurat / SCE object or the list from performPCA().",
    fixed = TRUE
  )
})
