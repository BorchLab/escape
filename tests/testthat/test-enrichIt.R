# test script for enrichIt.R - testcases are NOT comprehensive!

# ---- test fixtures ---------------------------------------------------------
gene_sets <- list(
  PathA = c("G1", "G2", "G10"),
  PathB = c("G3", "G4", "G5"),
  PathC = c("G8", "G9")
)

#  numeric vector (already ranked) -----------
vec <- c(G1 = 2.3, G2 = -1.1, G3 = 0.8, G4 = -2.2, G5 = 1.5)
vec <- sort(vec, decreasing = TRUE)

#  data-frame (Seurat-like) ------------------
de_tbl <- data.frame(
  gene        = paste0("G", 1:6),
  avg_log2FC  = c( 2.3, -1.1, 0.8, -2.2, 1.5, 0),
  p_val_adj   = c( 1e-4, 5e-3, 0.03, 1e-10, 0.2, 1),
  stringsAsFactors = FALSE
)
rownames(de_tbl) <- de_tbl$gene          # default layout

# ---- 1.  BASIC FUNCTIONALITY ----------------------------------------------
test_that("numeric vector input returns a proper fgsea table", {
  res <- enrichIt(vec, gene_sets)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pathway", "NES", "pval", "padj", "leadingEdge") %in%
                    names(res)))
  expect_true(is.character(res$leadingEdge))
  expect_true(all(res$padj >= 0 & res$padj <= 1))
})

test_that("data-frame input (default columns) works", {
  res <- enrichIt(de_tbl, gene_sets)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pathway", "NES", "pval") %in% names(res)))
})

# ---- 2.  ALTERNATIVE OPTIONS ----------------------------------------------
test_that("custom gene_col + explicit ranking_fun = 'logFC' works", {
  tbl <- de_tbl
  names(tbl)[names(tbl) == "avg_log2FC"] <- "logFC"
  res <- enrichIt(tbl,
                     gene.sets   = gene_sets,
                     gene_col    = "gene",
                     logFC_col   = "logFC",
                     pval_col    = "p_val_adj",
                     ranking_fun = "logFC")
  expect_s3_class(res, "data.frame")
})


# ---- 3.  ERROR HANDLING ----------------------------------------------------
test_that("error when no genes left after filtering", {
  expect_error(
    enrichIt(de_tbl,
                gene_sets,
                pval_cutoff  = 1e-10,
                logFC_cutoff = 10),
    "No genes left"
  )
})

test_that("error for unlabeled numeric vector", {
  bad_vec <- unname(vec)
  expect_error(
    enrichIt(bad_vec, gene_sets),
    "named numeric vector"
  )
})

test_that("error when required columns are missing", {
  tmp <- de_tbl
  tmp$avg_log2FC <- NULL
  expect_error(
    enrichIt(tmp, gene_sets),
    "logFC_col"
  )
})
