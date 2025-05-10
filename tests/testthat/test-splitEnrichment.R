# test script for splitEnrichment.R - testcases are NOT comprehensive!
## helper ----------------------------------------------------------------
geom_names <- function(p) vapply(p$layers, \(x) class(x$geom)[1], character(1))

## fixture ---------------------------------------------------------------
seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")

# ────────────────────────────────────────────────────────────────────────
test_that("returns a ggplot and uses split violins for two levels", {
  
  p <- splitEnrichment(
    seuratObj,
    assay    = "escape",
    split.by = "groups",      # has exactly 2 levels
    gene.set = "Tcells"
  )
  
  expect_s3_class(p, "ggplot")
  expect_true("GeomSplitViolin" %in% geom_names(p))
  expect_false("GeomViolin"      %in% geom_names(p))
})

# ────────────────────────────────────────────────────────────────────────
test_that("uses dodged violins when split.by has >2 levels", {
  
  # add a 3-level grouping variable
  seuratObj$groups3 <- rep(LETTERS[1:3], length.out = ncol(seuratObj))
  
  p <- splitEnrichment(
    seuratObj,
    assay    = "escape",
    split.by = "groups3",     # 3 levels
    gene.set = "Tcells"
  )
  
  expect_true("GeomViolin"      %in% geom_names(p))
  expect_false("GeomSplitViolin" %in% geom_names(p))
})

# ────────────────────────────────────────────────────────────────────────
test_that("scale = TRUE centres the values (≈ mean 0)", {
  
  p  <- splitEnrichment(
    seuratObj,
    assay    = "escape",
    split.by = "groups",
    gene.set = "Tcells",
    scale    = TRUE
  )
  
  yvals <- ggplot_build(p)$data[[1]]$y
  expect_lt(abs(mean(yvals, na.rm = TRUE)), 1e-6)
})

# ────────────────────────────────────────────────────────────────────────
test_that("order.by = 'mean' reorders x-axis levels by descending mean", {
  
  p <- splitEnrichment(
    seuratObj,
    assay    = "escape",
    split.by = "groups",
    gene.set = "Tcells",
    order.by = "mean"
  )
  
  ## compute expected order
  enr <- escape:::.prepData(
    input.data = seuratObj,
    assay      = "escape",
    gene.set   = "Tcells",
    group.by   = "ident",
    split.by   = "groups",
    facet.by   = NULL
  )
  
  expected <- enr %>%
    group_by(ident) %>%
    summarise(mu = mean(.data$Tcells)) %>%
    arrange(desc(mu)) %>%
    pull(ident) %>%
    as.character()
  
  expect_equal(levels(p$data$ident), expected)
})

# ────────────────────────────────────────────────────────────────────────
test_that("missing split.by argument triggers an error", {
  
  expect_error(
    splitEnrichment(
      seuratObj,
      assay    = "escape",
      gene.set = "Tcells"
    ),
    "split.by"   # error message should mention the missing argument
  )
})