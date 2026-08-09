# Toy SpatialExperiment fixture.
#
# Built on the fly rather than shipped as testdata: the object is tiny, fully
# deterministic, and keeps SpatialExperiment out of the hard test dependencies.
# Everything is fully qualified so sourcing this helper never fails when
# SpatialExperiment is absent - call skip_if_not_installed("SpatialExperiment")
# inside the test_that() block instead.

make_toy_spe <- function(n.genes = 40, n.cells = 60, n.samples = 2,
                         seed = 42) {
  set.seed(seed)
  cnts <- Matrix::rsparsematrix(
    n.genes, n.cells, density = 0.4,
    rand.x = function(n) stats::rpois(n, 5) + 1
  )
  dimnames(cnts) <- list(paste0("gene", seq_len(n.genes)),
                         paste0("cell", seq_len(n.cells)))

  # deterministic log values so tests never depend on scuttle/scater
  logc <- log1p(cnts)

  SpatialExperiment::SpatialExperiment(
    assays        = list(counts = cnts, logcounts = logc),
    colData       = S4Vectors::DataFrame(
      group = rep(c("a", "b"), length.out = n.cells),
      row.names = colnames(cnts)
    ),
    spatialCoords = matrix(
      seq_len(2 * n.cells), ncol = 2,
      dimnames = list(colnames(cnts), c("x", "y"))
    ),
    sample_id     = rep(paste0("s", seq_len(n.samples)), length.out = n.cells)
  )
}

# Underscored name on purpose: HALLMARK_/GO_/REACTOME_ sets are the case that
# used to be silently dropped during normalization for non-Seurat input.
toy_spe_sets <- function() {
  list(HALLMARK_SET_A = paste0("gene", 1:8),
       SetB           = paste0("gene", 9:16))
}

# Matching SingleCellExperiment, for the regression guard that the same bug was
# never spatial-specific.
make_toy_sce <- function(...) {
  spe <- make_toy_spe(...)
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts    = SummarizedExperiment::assay(spe, "counts"),
                   logcounts = SummarizedExperiment::assay(spe, "logcounts")),
    colData = SummarizedExperiment::colData(spe)[, "group", drop = FALSE]
  )
}
