# test script for enrichItPlot.R - testcases are NOT comprehensive!

skip_if_not_installed("fgsea")
skip_if_not_installed("patchwork")       

# helper data: run a very small fgsea ----------------------------------------
set.seed(42)

## ranked statistic ---------------------------------------------------------
gene_ids <- paste0("G", 1:80)
stat_vec <- setNames(rev(seq_along(gene_ids)), gene_ids)   # 80 .. 1 (descending)

## synthetic multi-library gene sets ----------------------------------------
gene_sets <- list(
  DB1_PathA = paste0("G",  1:15),
  DB1_PathB = paste0("G", 16:30),
  DB2_PathC = paste0("G", 21:35),  # overlaps with both A & B -> ensures cnet links
  DB2_PathD = paste0("G", 46:60)
)

res <- enrichIt(input.data = stat_vec,
                 gene.sets = gene_sets,
                 minSize   = 5           
)


res$Database <- ifelse(grepl("^DB1_", res$pathway), "DB1", "DB2")


# 1.  BAR plot ---------------------------------------------------------------
test_that("bar plot returns a patchwork object with ggplot inside", {
  plt <- enrichItPlot(res, plot.type = "bar", top = 3)
  
  expect_s3_class(plt, "patchwork")
  expect_true(inherits(plt[[1]], "ggplot"))
})

# ---------------------------------------------------------------------------
# 2.  DOT plot ---------------------------------------------------------------
test_that("dot plot returns a patchwork object and respects top argument", {
  plt <- enrichItPlot(res, plot.type = "dot", top = 1)
  
  expect_s3_class(plt, "patchwork")
  # only one term per database should survive top = 1
  build <- ggplot2::ggplot_build(plt[[1]])
  n_terms <- length(unique(build$data[[1]]$y))
  expect_lte(n_terms, 2)          # 2 databases  ⇒ ≤2 rows in panel 1
})

# ---------------------------------------------------------------------------
# 3.  CNET plot --------------------------------------------------------------
test_that("cnet plot returns a ggraph object", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("igraph")
  
  plt <- enrichItPlot(res, plot.type = "cnet", top = 4)
  
  expect_s3_class(plt, "ggraph")
})

# ---------------------------------------------------------------------------
# 4.  Error handling ---------------------------------------------------------
test_that("invalid plot.type triggers an informative error", {
  expect_error(
    enrichItPlot(res, plot.type = "heatmap"),
    regexp = "cnet"
  )
})
