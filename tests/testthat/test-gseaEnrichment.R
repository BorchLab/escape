# test script for gseaEnrichment.R - testcases are NOT comprehensive!

#####  Helper: tiny toy dataset ------------------------------------------ ###
toy_mat <- matrix(c(
  # Gene1 Gene2 Gene3 Gene4 Gene5
  10,   20,    1,    2,   30,   # group A cell 1
  11,   21,    1,    1,   29,   # group A cell 2
  2,    1,   25,   22,    3,   # group B cell 1
  1,    2,   24,   21,    4    # group B cell 2
), nrow = 5, byrow = FALSE,
dimnames = list(
  paste0("Gene", 1:5),
  paste0("Cell", 1:4)
))

toy_groups <- factor(c("A", "A", "B", "B"))
toy_gs     <- list(Pathway = c("Gene1", "Gene3", "Gene5"))

# Expected ES for group A: leading genes 1 & 5 are in gene-set → positive peak
# Expected ES for group B: leading genes 3 is in gene-set → positive peak
# We just assert sign (+) and non-zero magnitude.

##### 1.  Function runs and returns ggplot / patchwork -------------------- ###
test_that("basic run (Seurat) returns a patchwork plot with ES in legend", {
  seu <- CreateSeuratObject(counts = toy_mat)
  seu$grp <- toy_groups
  
  plt <- gseaEnrichment(seu,
                        gene.set.use = "Pathway",
                        gene.sets    = toy_gs,
                        group.by     = "grp")
  
  expect_s3_class(plt, "patchwork")
  # ggplot object exists inside
  expect_true(inherits(plt[[1]], "ggplot"))
  
  # Legend label contains ES =
  build <- ggplot_build(plt[[1]])
  labs  <- build$plot$scales$scales[[1]]$get_labels()
  expect_true(any(grepl("ES\\s*=\\s*", labs)))
})

##### 2.  Works on SummarizedExperiment ----------------------------------- ###
test_that("basic run (SummarizedExperiment) works", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = toy_mat),
    colData = data.frame(grp = toy_groups))
  
  plt <- gseaEnrichment(se,
                        gene.set.use = "Pathway",
                        gene.sets    = toy_gs,
                        group.by     = "grp",
                        summary.fun  = "median")
  
  expect_s3_class(plt, "patchwork")
})

##### 3.  All built-in summary.fun keywords + custom ---------------------- ###
keys <- c("mean", "median", "max", "sum", "geometric")
for (k in keys) {
  test_that(paste("summary.fun =", k, "runs"), {
    seu <- CreateSeuratObject(counts = toy_mat); seu$grp <- toy_groups
    expect_silent(
      gseaEnrichment(seu,
                     gene.set.use = "Pathway",
                     gene.sets    = toy_gs,
                     group.by     = "grp",
                     summary.fun  = k)
    )
  })
}
test_that("custom summary.fun runs", {
  seu <- CreateSeuratObject(counts = toy_mat); seu$grp <- toy_groups
  expect_silent(
    gseaEnrichment(seu,
                   gene.set.use = "Pathway",
                   gene.sets    = toy_gs,
                   group.by     = "grp",
                   summary.fun  = sd)
  )
})

##### 4.  Numerical sanity: ES sign & non-zero ---------------------------- ###
test_that("enrichment score is positive and non-zero for toy data", {
  seu <- CreateSeuratObject(counts = toy_mat); seu$grp <- toy_groups
  plt <- gseaEnrichment(seu,
                        gene.set.use = "Pathway",
                        gene.sets    = toy_gs,
                        group.by     = "grp",
                        digits       = 4)
  
  labs <- ggplot_build(plt[[1]])$plot$scales$scales[[1]]$get_labels()
  es_vals <- as.numeric(sub(".*ES\\s*=\\s*([0-9.+-]+).*", "\\1", labs))
  expect_true(all(es_vals > 0))
})

##### 5.  Error handling --------------------------------------------------- ###
seu_base <- CreateSeuratObject(counts = toy_mat); seu_base$grp <- toy_groups

test_that("errors for multiple gene-set names", {
  expect_error(
    gseaEnrichment(seu_base,
                   gene.set.use = c("x","y"),
                   gene.sets    = toy_gs,
                   group.by     = "grp"),
    "length 1"
  )
})

test_that("errors for unknown gene-set", {
  expect_error(
    gseaEnrichment(seu_base,
                   gene.set.use = "Unknown",
                   gene.sets    = toy_gs,
                   group.by     = "grp"),
    "Unknown gene-set"
  )
})

test_that("errors when <2 groups", {
  seu1 <- seu_base[,1:2]   # only group A
  expect_error(
    gseaEnrichment(seu1,
                   gene.set.use = "Pathway",
                   gene.sets    = toy_gs,
                   group.by     = "grp"),
    "Need ≥2 groups"
  )
})

test_that("errors for zero overlap gene-set", {
  bad_gs <- list(Bad = c("NotInMatrix"))
  expect_error(
    gseaEnrichment(seu_base,
                   gene.set.use = "Bad",
                   gene.sets    = bad_gs,
                   group.by     = "grp"),
    "overlap"
  )
})

test_that("errors when group.by column missing", {
  expect_error(
    gseaEnrichment(seu_base,
                   gene.set.use = "Pathway",
                   gene.sets    = toy_gs,
                   group.by     = "missing"),
    "not found"
  )
})
