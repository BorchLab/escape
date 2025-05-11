# test script for geyserEnrichment.R - testcases are NOT comprehensive!

# ────────────────────────────────────────────────────────────────────────────────
#  Test-data set-up -------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(SeuratObject)
})

pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")

# helper to make repeated plotting calls tidy
plot_fun <- function(...) {
  geyserEnrichment(pbmc_small, assay = "escape", ...)
}

# ────────────────────────────────────────────────────────────────────────────────
#  Core object / mapping checks --------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("default call returns a ggplot object with expected mappings", {
  p <- plot_fun(gene.set = "Tcells")
  
  expect_s3_class(p, "ggplot")
  
  # x-axis should map to ident (default group.by)
  expect_identical(
    rlang::get_expr(p$mapping$x),
    rlang::expr(.data[["ident"]])
  )
  
  # y-axis should map to the chosen gene-set
  expect_identical(
    rlang::get_expr(p$mapping$y),
    rlang::expr(.data[["Tcells"]])
  )
})

# ────────────────────────────────────────────────────────────────────────────────
#  order.by logic ----------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("order.by = 'mean' sorts x-axis levels by group mean", {
  p <- plot_fun(gene.set = "Tcells", order.by = "mean")
  
  d <- p$data
  means <- tapply(d$Tcells, d$ident, mean, na.rm = TRUE)
  expect_identical(levels(d$ident), names(sort(means)))
})

test_that("invalid order.by triggers an informative error", {
  expect_error(
    plot_fun(gene.set = "Tcells", order.by = "bogus"),
    "order.by must be 'mean' or 'group.by'"
  )
})

# ────────────────────────────────────────────────────────────────────────────────
#  scale = TRUE (z-transformation) ----------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("scale = TRUE centres and scales the enrichment distribution", {
  p <- plot_fun(gene.set = "Tcells", scale = TRUE)
  z <- p$data$Tcells
  
  expect_lt(abs(mean(z, na.rm = TRUE)), 1e-6)       # ~0
  expect_lt(abs(sd(z,   na.rm = TRUE) - 1), 1e-6)   # ~1
})

# ────────────────────────────────────────────────────────────────────────────────
#  colour handling --------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("colour.by = 'group' creates a discrete colour scale", {
  p <- plot_fun(gene.set = "Tcells", color.by = "group")  # maps to ident
  scale_classes <- vapply(p$scales$scales, class, character(1))
  expect_true(any(grepl("ScaleDiscrete", scale_classes)))
})

test_that("colouring by the gene-set itself yields a continuous scale", {
  p <- plot_fun(gene.set = "Tcells", color.by = "Tcells")
  scale_classes <- vapply(p$scales$scales, class, character(1))
  expect_true(any(grepl("ScaleContinuous", scale_classes)))
})

# ────────────────────────────────────────────────────────────────────────────────
#  facetting --------------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("facet.by adds a FacetGrid object", {
  p <- plot_fun(gene.set = "Tcells", facet.by = "groups")
  expect_s3_class(p$facet, "FacetGrid")
})

# ────────────────────────────────────────────────────────────────────────────────
#  edge-case & robustness checks ------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
test_that("missing gene-set triggers an error", {
  expect_error(
    plot_fun(gene.set = "NonExistentGeneSet"),
    "not found|missing|must be present"
  )
})

test_that("missing group.by column triggers an error", {
  expect_error(
    plot_fun(gene.set = "Tcells", group.by = "unknown_column"),
    "Expecting a Seurat or SummarizedExperiment object|column"
  )
})
