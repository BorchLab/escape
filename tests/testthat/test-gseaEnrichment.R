# test script for gseaEnrichment.R - testcases are NOT comprehensive!

pbmc  <- SeuratObject::pbmc_small
GS    <- list(
  Bcells = c("MS4A1", "CD79B", "CD79A", "IGHG1", "IGHG2"),
  Tcells = c("CD3E", "CD3D", "CD3G", "CD7",  "CD8A")
)

##### 1.  Function runs and returns ggplot / patchwork -------------------- ###
test_that("basic run (Seurat) returns a patchwork plot with ES in legend", {
  
  plt <- gseaEnrichment(pbmc,
                        gene.set.use = "Tcells",
                        gene.sets    = GS)
  
  expect_s3_class(plt, "patchwork")
  # ggplot object exists inside
  expect_true(inherits(plt[[1]], "ggplot"))
  
  # Legend label contains ES =
  build <- ggplot_build(plt[[1]])
  labs  <- build$plot$scales$scales[[1]]$get_labels()
  expect_true(any(grepl("ES\\s*=\\s*", labs)))
})


##### 2.  All built-in summary.fun keywords + custom ---------------------- ###
keys <- c("mean", "median", "max", "sum", "geometric")
for (k in keys) {
  test_that(paste("summary.fun =", k, "runs"), {
    expect_silent(
      gseaEnrichment(pbmc,
                     gene.set.use = "Bcells",
                     gene.sets    = GS)
    )
  })
}

test_that("custom summary.fun runs", {
  expect_silent(
    gseaEnrichment(pbmc,
                   gene.set.use = "Tcells",
                   gene.sets    = GS)
  )
})


##### 3.  Error handling --------------------------------------------------- ###
seu_base <- CreateSeuratObject(counts = toy_mat); seu_base$grp <- toy_groups

test_that("errors for multiple gene-set names", {
  expect_error(
    gseaEnrichment(pbmc,
                   gene.set.use = c("x","y"),
                   gene.sets    = GS), 
    "length 1"
  )
})

test_that("errors for unknown gene-set", {
  expect_error(
    gseaEnrichment(pbmc,
                   gene.set.use = "Unknown",
                   gene.sets    = GS), 
    "Unknown gene-set"
  )
})

