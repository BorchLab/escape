# test script for scatterEnrichment.R - testcases are NOT comprehensive!

# ---------------------------------------------------------------------------
#  Load test data ------------------------------------------------------------
# ---------------------------------------------------------------------------
pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")   # helper provided by escape
x.gene <- "Tcells"
y.gene <- "Bcells"

# ---------------------------------------------------------------------------
# 1. Argument validation -----------------------------------------------------
# ---------------------------------------------------------------------------
test_that("invalid 'style' argument throws error", {
  expect_error(
    scatterEnrichment(pbmc_small,
                      assay = "escape", x.axis = x.gene, y.axis = y.gene,
                      style = "foo"),
    regexp = "point"
  )
})

test_that("invalid 'color.by' argument throws error", {
  expect_error(
    scatterEnrichment(pbmc_small,
                      assay = "escape", x.axis = x.gene, y.axis = y.gene,
                      color.by = "foobar"),
    regexp = "density"
  )
})

# ---------------------------------------------------------------------------
# 2. Object type -------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("function returns a ggplot object", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", 
                         x.axis = x.gene, 
                         y.axis = y.gene, 
                         color.by   = "density",
                         style  = "point")
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 3. Layer composition -------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("style = 'point' adds GeomPointdensity layer", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         style = "point")
  geoms <- vapply(p$layers, \(l) class(l$geom)[1], character(1))
  expect_true("GeomPoint" %in% geoms)
})

test_that("style = 'hex' adds StatBinhex layer", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         style = "hex")
  stats <- vapply(p$layers, \(l) class(l$stat)[1], character(1))
  expect_true("StatBinhex" %in% stats)
})

# ---------------------------------------------------------------------------
# 4. Scaling option ----------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("'scale = TRUE' centres and scales gene-set columns", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         scale = TRUE)
  m1 <- mean(p$data[[x.gene]])
  s1 <-  sd(p$data[[x.gene]])
  m2 <- mean(p$data[[y.gene]])
  s2 <-  sd(p$data[[y.gene]])
  expect_lt(abs(m1), 1e-6)
  expect_lt(abs(m2), 1e-6)
  expect_equal(round(s1, 6), 1)
  expect_equal(round(s2, 6), 1)
})

# ---------------------------------------------------------------------------
# 5. Facetting ---------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("facet.by generates expected facets", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         facet.by = "letter.idents")
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(
    sort(unique(p$data$letter.idents)),
    sort(unique(pbmc_small$letter.idents))
  )
})

# ---------------------------------------------------------------------------
# 6. Coloring strategies ----------------------------------------------------
# ---------------------------------------------------------------------------
test_that("color.by = 'group' maps discrete colour aesthetic", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         color.by = "group", group.by = "groups")
  map_vars <- union(names(p$mapping), names(p$layers[[1]]$mapping))
  expect_true("colour" %in% tolower(map_vars))
})

test_that("color.by = 'x' produces continuous colour scale", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         color.by = "x")
  cont_scale <- any(vapply(
    p$scales$scales,
    \(s) inherits(s, "ScaleColourGradient"),
    logical(1)
  ))
  expect_false(cont_scale)
})

# ---------------------------------------------------------------------------
# 7. Correlation overlay -----------------------------------------------------
# ---------------------------------------------------------------------------
test_that("add.corr inserts a GeomText annotation layer", {
  p <- scatterEnrichment(pbmc_small,
                         assay = "escape", x.axis = x.gene, y.axis = y.gene,
                         add.corr = TRUE)
  has_text <- any(vapply(
    p$layers, \(l) inherits(l$geom, "GeomText"),
    logical(1)
  ))
  expect_true(has_text)
})
