# test script for ridgeEnrichment.R - testcases are NOT comprehensive!

pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")

# -------------------------------------------------------------------------
test_that("returns a proper ggplot object", {
  
  
  p <- ridgeEnrichment(
    pbmc_small,
    assay     = "escape",
    gene.set  = "Tcells",
    group.by  = "groups"
  )
  
  expect_s3_class(p, "ggplot")
  # at least one ridge geom layer (gradient or non-gradient)
  ridge_layers <- vapply(
    p$layers,
    \(ly) inherits(ly$geom,
                   c("GeomDensityRidges", "GeomDensityRidgesGradient")),
    logical(1)
  )
  expect_true(any(ridge_layers))
})

# -------------------------------------------------------------------------
test_that("gradient colour mode when colour.by == gene.set", {
  p <- ridgeEnrichment(
    pbmc_small, assay = "escape",
    gene.set  = "Tcells",
    color.by  = "Tcells"   # triggers numeric gradient
  )
  # mapping$fill should be after_stat(x)
  expect_equal(rlang::quo_text(p$mapping$fill), "if (gradient.mode) ggplot2::after_stat(x) else .data[[\"Tcells\"]]")
})

# -------------------------------------------------------------------------
test_that("categorical colour mode when colour.by == group", {
  p <- ridgeEnrichment(
    pbmc_small, assay = "escape",
    gene.set  = "Tcells",
    color.by  = "group",        # will internally map to group.by "groups"
    group.by  = "groups"
  )
  expect_equal(rlang::quo_text(p$mapping$fill), "if (gradient.mode) ggplot2::after_stat(x) else .data[[\"groups\"]]")
})

# -------------------------------------------------------------------------
test_that("scale = TRUE centres distribution at zero", {
  p <- ridgeEnrichment(
    pbmc_small, assay = "escape",
    gene.set = "Tcells",
    scale    = TRUE
  )
  m <- mean(p$data$Tcells, na.rm = TRUE)
  expect_lt(abs(m), 1e-8)
})

# -------------------------------------------------------------------------
test_that("order.by = 'mean' re-orders factor levels by mean score", {
  p <- ridgeEnrichment(
    pbmc_small, assay = "escape",
    gene.set = "Tcells",
    group.by = "groups",
    order.by = "mean"
  )
  grp      <- p$data$groups
  grp_means <- tapply(p$data$Tcells, grp, mean)
  # levels should be sorted by increasing mean
  expect_equal(levels(grp), names(rev(sort(grp_means))))
})

# -------------------------------------------------------------------------
test_that("add.rug = TRUE switches on jittered points", {
  p <- ridgeEnrichment(
    pbmc_small, assay = "escape",
    gene.set = "Tcells",
    add.rug  = TRUE
  )
  expect_true(any(vapply(
    p$layers,
    \(ly) isTRUE(ly$stat_params$jittered_points),
    logical(1)
  )))
})
