# test script for densityEnrichment.R - testcases are NOT comprehensive!

pbmc  <- SeuratObject::pbmc_small
GS    <- list(
  Bcells = c("MS4A1", "CD79B", "CD79A", "IGHG1", "IGHG2"),
  Tcells = c("CD3E", "CD3D", "CD3G", "CD7",  "CD8A")
)

# helper: number of groups in default 'ident' column
n_groups <- length(unique(as.character(Idents(pbmc))))


# ── 1  Core functionality returns patchwork object ──────────────────
test_that("densityEnrichment() returns a patchwork / ggplot object", {
  plt <- densityEnrichment(
    input.data   = pbmc,
    gene.set.use = "Tcells",
    gene.sets    = GS
  )
  
  expect_s3_class(plt, "patchwork")          # overall object
  expect_s3_class(plt[[1]], "ggplot")        # top density panel
  expect_s3_class(plt[[2]], "ggplot")        # bottom segment panel
})

# ── 2  Groups are represented correctly in the density plot ─────────
test_that("all groups appear once in the density layer", {
  plt         <- densityEnrichment(pbmc, "Tcells", GS)
  density_pan <- plt[[1]]
  vars_in_df  <- unique(density_pan$data$variable)
  
  expect_equal(length(na.omit(vars_in_df)), n_groups)
})

# ── 3  Alternative palettes run without error ───────────────────────
test_that("custom palette works", {
  expect_no_error(
    densityEnrichment(
      input.data   = pbmc,
      gene.set.use = "Bcells",
      gene.sets    = GS,
      palette      = "viridis"
    )
  )
})

# ── 4  Input validation – wrong object or gene-set names ────────────
test_that("input validation errors are triggered correctly", {
  mat <- matrix(rpois(1000, 5), nrow = 100)   # not a single-cell object
  
  expect_error(
    densityEnrichment(mat, "Tcells", GS),
    "Expecting a Seurat or SummarizedExperiment object"
  )
})

# ── 5  group.by argument overrides default ---------------------------
test_that("group.by selects an alternative metadata column", {
  pbmc$dummy_group <- sample(c("A", "B"), ncol(pbmc), replace = TRUE)
  plt <- densityEnrichment(
    pbmc, gene.set.use = "Tcells", gene.sets = GS,
    group.by = "dummy_group"
  )
  
  # check that new grouping made it into the plot data
  expect_true(all(grepl("^dummy_group\\.", na.omit(unique(plt[[1]]$data$variable)))))
})
