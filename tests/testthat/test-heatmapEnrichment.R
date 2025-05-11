# test script for heatmapEnrichment.R - testcases are NOT comprehensive!

test_that("setup: example dataset is available", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("escape")        # runEscape & helpers
  expect_silent(
    seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  )
  expect_s3_class(seuratObj, "Seurat")
})

# ----------------------------------------------------------------
#  1. Basic functionality & return type
# ----------------------------------------------------------------
test_that("default call returns a ggplot object", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  p <- heatmapEnrichment(seuratObj, assay = "escape")
  expect_s3_class(p, "ggplot")
  expect_true(c("group", "variable", "value") %in% names(p$data))
  # default summary = mean; check at least one numeric value present
  expect_true(is.numeric(p$data$value))
})

# ----------------------------------------------------------------
#  2. Gene-set sub-selection
# ----------------------------------------------------------------
test_that("gene.set.use filters rows correctly", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  chosen <- c("Bcells", "Tcells")
  p <- heatmapEnrichment(seuratObj,
                         assay        = "escape",
                         gene.set.use = chosen)
  expect_setequal(unique(p$data$variable), chosen)
})

# ----------------------------------------------------------------
#  3. Scaling (Z-transform)
# ----------------------------------------------------------------
test_that("scale = TRUE centres each gene set to mean ≈ 0", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  p <- heatmapEnrichment(seuratObj,
                         assay = "escape",
                         scale = TRUE)
  z_by_gene <- split(p$data$value, p$data$variable)
  # Mean of each scaled column should be 0 (tolerance for FP error)
  z_means <- vapply(z_by_gene, mean, numeric(1))
  expect_true(all(abs(z_means) < 1e-6))
})

# ----------------------------------------------------------------
#  4. Summary statistics (median, custom, error handling)
# ----------------------------------------------------------------
test_that("summary.stat = 'median' gives expected result", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  gs <- "Bcells"
  # Manual median for reference
  x <- FetchData(seuratObj, vars = gs, slot = "data", assay = "escape")[, 1]
  grp <- Idents(seuratObj)
  ref_median <- tapply(x, grp, median)
  p <- heatmapEnrichment(seuratObj,
                         assay         = "escape",
                         gene.set.use  = gs,
                         summary.stat  = "median")
  # Extract tile corresponding to first group
  med_calc <- subset(p$data,
                     variable == gs & group == names(ref_median)[1])$value
  expect_equal(med_calc, unname(ref_median[1]), tolerance = 1e-8)
})

test_that("custom summary function is accepted", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  range_fun <- function(x) max(x) - min(x)
  p <- heatmapEnrichment(seuratObj,
                         assay        = "escape",
                         summary.stat = range_fun)
  expect_s3_class(p, "ggplot")
})

test_that("invalid summary keyword errors cleanly", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  expect_error(
    heatmapEnrichment(seuratObj,
                      assay        = "escape",
                      summary.stat = "foobar"),
    "Unsupported summary keyword"
  )
})

# ----------------------------------------------------------------
#  5. Clustering options
# ----------------------------------------------------------------
test_that("row/column clustering re-orders factors", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  p <- heatmapEnrichment(seuratObj,
                         assay            = "escape",
                         cluster.rows     = TRUE,
                         cluster.columns  = TRUE)
  # After clustering, factors keep their specified order
  expect_true(is.factor(p$data$variable))
  expect_true(is.factor(p$data$group))
})

# ----------------------------------------------------------------
#  6. Faceting
# ----------------------------------------------------------------
test_that("facet.by adds facetting column to output", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  p <- heatmapEnrichment(seuratObj,
                         assay    = "escape",
                         facet.by = "letter.idents")
  expect_true("letter.idents" %in% names(p$data))
  # ggplot2 stores facet mapping in the plot's Facets object
  expect_true(inherits(p$facet, "Facet"))
})

# ----------------------------------------------------------------
#  7. Argument validation
# ----------------------------------------------------------------
test_that("unknown gene set triggers informative error", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  expect_error(
    heatmapEnrichment(seuratObj,
                      assay        = "escape",
                      gene.set.use = "NonExistentGS"),
    "No gene-set columns found"
  )
})