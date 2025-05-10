# test script for utils.R - testcases are NOT comprehensive!

## --------------------------------------------------------------------- ##
##  1. Fast negation operator                                            ##
## --------------------------------------------------------------------- ##
test_that("%!in% negates %in% correctly", {
  x <- 1:5
  y <- 3:7
  expect_identical(x %!in% y, !(x %in% y))
})

## --------------------------------------------------------------------- ##
##  2. Class helpers & .checkSingleObject                                ##
## --------------------------------------------------------------------- ##
test_that("class helpers recognise Seurat / SCE", {
  # Seurat branch -------------------------------------------------------
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    seurat_obj <- SeuratObject::CreateSeuratObject(
      counts = matrix(rpois(20, 1), nrow = 4)
    )
    expect_true(.is_seurat(seurat_obj))
    expect_false(.is_sce(seurat_obj))
    expect_true(.is_seurat_or_sce(seurat_obj))
  }
  
  # SCE branch ----------------------------------------------------------
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = matrix(rpois(20, 1), nrow = 4))
    )
    expect_true(.is_sce(sce))
    expect_false(.is_seurat(sce))
    expect_true(.is_seurat_or_sce(sce))
  }
  
  # Generic error -------------------------------------------------------
  expect_error(.checkSingleObject(list()), "Expecting a Seurat or")
})

## --------------------------------------------------------------------- ##
##  3. .orderFunction                                                   ##
## --------------------------------------------------------------------- ##
test_that(".orderFunction orders by mean correctly", {
  df <- data.frame(value = c(5, 1, 2, 8, 4, 7),
                   grp   = c("A", "B", "A", "C", "B", "C"))
  out <- .orderFunction(df, order.by = "mean", group.by = "grp")
  expect_equal(levels(out$grp), c("C", "A", "B"))  # means 7.5 > 3.5 > 2.5
})

test_that(".orderFunction gives natural alpha-numeric order", {
  df <- data.frame(value = 1:6,
                   bucket = c("G1", "G2", "G10", "G11", "G3", "G20"))
  out <- .orderFunction(df, order.by = "group.by", group.by = "bucket")
  expect_equal(levels(out$bucket)[1:4], c("G1", "G2", "G3", "G10"))
})

test_that(".orderFunction input validation works", {
  expect_error(.orderFunction(data.frame(x = 1), "foo", "x"),
               "order.by must be")
})

## --------------------------------------------------------------------- ##
##  4. Splitters                                                         ##
## --------------------------------------------------------------------- ##
test_that(".split_cols splits into roughly equal column chunks", {
  mat <- matrix(seq_len(20), nrow = 4)   # 4 × 5
  out <- .split_cols(mat, chunk = 2)
  expect_length(out, 3)                  # 2+2+1 columns
  expect_equal(ncol(out[[1]]), 2)
  expect_equal(ncol(out[[3]]), 1)
})

test_that(".split_rows splits rows and preserves data", {
  mat <- matrix(seq_len(20), nrow = 10, ncol = 2)
  out <- .split_rows(mat, chunk.size = 3)
  expect_length(out, 4)                  # 3+3+3+1 rows
  expect_equal(nrow(out[[4]]), 1)
  expect_equal(rbind(out[[1]], out[[2]], out[[3]], out[[4]]), mat)
})

test_that(".split_vector chunks vectors", {
  v <- letters[1:11]
  out <- .split_vector(v, chunk.size = 4)
  expect_equal(lengths(out), c(4, 4, 3))
  expect_equal(unlist(out), v)
})

## --------------------------------------------------------------------- ##
##  5. .colorizer & .colorby                                             ##
## --------------------------------------------------------------------- ##
test_that(".colorizer returns n distinct colours", {
  pal <- .colorizer("viridis", n = 5)
  expect_length(pal, 5)
  expect_true(all(!is.na(pal)))
})

test_that(".colorby adds gradient scale for numeric colour.by", {
  df <- data.frame(val = rnorm(4), group = letters[1:4])
  p  <- ggplot(df, aes(group, 1, fill = val)) + geom_col()
  p2 <- .colorby(df, p, color.by = "val", palette = "mako", type = "fill")
  expect_s3_class(p2, "ggplot")
  expect_true(any(vapply(p2$scales$scales,
                         inherits, logical(1), "ScaleContinuous")))
})

test_that(".colorby adds manual scale for categorical colour.by", {
  df <- data.frame(val = rnorm(4), group = c("C2", "C10", "C1", "C3"))
  p  <- ggplot(df, aes(group, 1, fill = group)) + geom_col()
  p2 <- .colorby(df, p, color.by = "group", palette = "plasma", type = "fill")
  expect_s3_class(p2, "ggplot")
  expect_true(any(vapply(p2$scales$scales,
                         inherits, logical(1), "ScaleDiscrete")))
})

## --------------------------------------------------------------------- ##
##  6. .cntEval                                                          ##
## --------------------------------------------------------------------- ##
test_that(".cntEval drops all-zero rows for plain matrices", {
  m <- matrix(c(0, 0, 1, 2, 0, 0), nrow = 3, byrow = TRUE,
              dimnames = list(paste0("g", 1:3), NULL))
  out <- .cntEval(m)
  expect_equal(rownames(out), c("g2", "g3"))
})

test_that(".cntEval works for Seurat & SCE (if installed)", {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    s <- SeuratObject::CreateSeuratObject(
      counts = matrix(c(0, 0, 1, 0, 3, 4), nrow = 3,
                      dimnames = list(c("g1", "g2", "g3"), NULL))
    )
    out <- .cntEval(s)
    expect_equal(rownames(out), c("g2", "g3"))
  }
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = matrix(c(0, 2, 0, 0, 0, 4), nrow = 3,
                                    dimnames = list(c("g1", "g2", "g3"), NULL)))
    )
    out <- .cntEval(sce)
    expect_equal(rownames(out), c("g1", "g3"))
  }
})

## --------------------------------------------------------------------- ##
##  7. .GS.check                                                         ##
## --------------------------------------------------------------------- ##
test_that(".GS.check validates input", {
  expect_error(.GS.check(NULL), "Please supply")
  expect_equal(.GS.check(list(A = c("a", "b"))), list(A = c("a", "b")))
  
  if (requireNamespace("GSEABase", quietly = TRUE)) {
    gs <- GSEABase::GeneSetCollection(
      GSEABase::GeneSet(setName = "foo", geneIds = c("x", "y"))
    )
    expect_equal(.GS.check(gs), list(foo = c("x", "y")))
  }
})

## --------------------------------------------------------------------- ##
##  8. .load_backend & _compute_enrichment                               ##
## --------------------------------------------------------------------- ##
test_that(".load_backend errors informatively", {
  expect_error(.load_backend("definitelyNotInstalledPackageXYZ"),
               "not installed")
})

test_that(".compute_enrichment rejects unknown method", {
  m <- matrix(rpois(20, 5), nrow = 5)
  expect_error(.compute_enrichment(m, gene_sets = list(s1 = letters[1:3]),
                                   method = "FOOBAR",
                                   BPPARAM = BiocParallel::SerialParam()),
               "Unknown method")
})

## --------------------------------------------------------------------- ##
##  9. Matrix column splitter (second copy at end of file)               ##
## --------------------------------------------------------------------- ##
test_that(".split_cols duplicate definition behaves consistently", {
  mat <- matrix(seq_len(12), nrow = 3)          # 3 × 4
  expect_identical(.split_cols(mat, 5), list(mat))  # <= chunk size
})
