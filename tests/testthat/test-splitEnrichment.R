# test script for splitEnrichment.R - testcases are NOT comprehensive!
## helper ----------------------------------------------------------------
geom_names <- function(p) vapply(p$layers, \(x) class(x$geom)[1], character(1))

## fixture ---------------------------------------------------------------
pbmc_small <- getdata("runEscape", "pbmc_small_ssGSEA")

# ────────────────────────────────────────────────────────────────────────
test_that("returns a ggplot and uses split violins for two levels", {
  
  p <- splitEnrichment(
    pbmc_small,
    assay    = "escape",
    split.by = "groups",      # has exactly 2 levels
    gene.set = "Tcells"
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(any(sapply(p$layers, function(layer) inherits(layer$geom, "GeomSplitViolin"))))
})

# ────────────────────────────────────────────────────────────────────────
test_that("uses dodged violins when split.by has >2 levels", {
  
  # add a 3-level grouping variable
  pbmc_small$groups3 <- rep(LETTERS[1:3], length.out = ncol(pbmc_small))
  
  p <- splitEnrichment(
    pbmc_small,
    assay    = "escape",
    split.by = "groups3",     
    gene.set = "Tcells"
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(!any(sapply(p$layers, function(layer) inherits(layer$geom, "GeomSplitViolin"))))
})

# ────────────────────────────────────────────────────────────────────────
test_that("scale = TRUE centres the values (≈ mean 0)", {
  
  p  <- splitEnrichment(
    pbmc_small,
    assay    = "escape",
    split.by = "groups",
    gene.set = "Tcells",
    scale    = TRUE
  )
  
  yvals <- ggplot_build(p)$data[[1]]$y
  expect_lt(abs(mean(yvals, na.rm = TRUE)), 1e-2)
})

# ────────────────────────────────────────────────────────────────────────
test_that("order.by = 'mean' reorders x-axis levels by descending mean", {
  
  p <- splitEnrichment(
    pbmc_small,
    assay    = "escape",
    split.by = "groups",
    gene.set = "Tcells",
    order.by = "mean"
  )
  
  ## compute expected order
  enr <- escape:::.prepData(
    input.data = pbmc_small,
    assay      = "escape",
    gene.set   = "Tcells",
    group.by   = "ident",
    split.by   = "groups",
    facet.by   = NULL
  )
  
  expected <- enr %>%
    dplyr::group_by(ident) %>%
    dplyr::summarise(mu = mean(.data$Tcells)) %>%
    dplyr::arrange(desc(mu)) %>%
    dplyr::pull(ident) %>%
    as.character()
  
  expect_equal(levels(p$data$ident), expected)
})

# ────────────────────────────────────────────────────────────────────────
test_that("missing split.by argument triggers an error", {
  
  expect_error(
    splitEnrichment(
      pbmc_small,
      assay    = "escape",
      gene.set = "Tcells"
    ),
    "split.by"   
  )
})