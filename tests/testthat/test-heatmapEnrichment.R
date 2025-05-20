# test script for heatmapEnrichment.R - testcases are NOT comprehensive!


pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")
# ----------------------------------------------------------------
#  1. Basic functionality & return type
# ----------------------------------------------------------------
test_that("default call returns a ggplot object", {
  p <- heatmapEnrichment(pbmc_small, assay = "escape")
  expect_s3_class(p, "ggplot")
  expect_true(all(c("group", "variable", "value") %in% names(p$data)))
  # default summary = mean; check at least one numeric value present
  expect_true(is.numeric(p$data$value))
})

# ----------------------------------------------------------------
#  2. Gene-set sub-selection
# ----------------------------------------------------------------
test_that("gene.set.use filters rows correctly", {
  chosen <- c("Bcells", "Tcells")
  p <- heatmapEnrichment(pbmc_small,
                         assay        = "escape",
                         gene.set.use = chosen)
  expect_setequal(unique(p$data$variable), chosen)
})

# ----------------------------------------------------------------
#  3. Scaling (Z-transform)
# ----------------------------------------------------------------
test_that("scale = TRUE centres each gene set to mean ≈ 0", {
  p <- heatmapEnrichment(pbmc_small,
                         assay = "escape",
                         scale = TRUE)
  z_by_gene <- split(p$data$value, p$data$variable)
  # Mean of each scaled column should be 0 (tolerance for FP error)
  z_means <- vapply(z_by_gene, mean, numeric(1))
  expect_true(all(abs(z_means) < 0.1))
})

# ----------------------------------------------------------------
#  4. Summary statistics (median, custom, error handling)
# ----------------------------------------------------------------
test_that("summary.stat = 'median' gives expected result", {
  gs <- "Bcells"
  # Manual median for reference
  x <- pbmc_small[["escape"]]@data[gs,]
  grp <- Idents(pbmc_small)
  ref_median <- tapply(x, grp, median)
  p <- heatmapEnrichment(pbmc_small,
                         assay         = "escape",
                         gene.set.use  = gs,
                         summary.stat  = "median")
  # Extract tile corresponding to first group
  med_calc <- subset(p$data,
                     variable == gs & group == names(ref_median)[1])$value
  expect_equal(med_calc, unname(ref_median[1]), tolerance = 1e-8)
})

test_that("invalid summary keyword errors cleanly", {
  expect_error(
    heatmapEnrichment(pbmc_small,
                      assay        = "escape",
                      summary.stat = "foobar"),
    "Unsupported summary keyword"
  )
})

# ----------------------------------------------------------------
#  5. Clustering options
# ----------------------------------------------------------------
test_that("row/column clustering re-orders factors", {
  p <- heatmapEnrichment(pbmc_small,
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
  p <- heatmapEnrichment(pbmc_small,
                         assay    = "escape",
                         facet.by = "letter.idents")
  expect_true("letter.idents" %in% names(p$data))
  # ggplot2 stores facet mapping in the plot's Facets object
  expect_true(inherits(p$facet, "Facet"))
})

