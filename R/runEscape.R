#' Calculate Single-Cell Gene-Set Enrichment Scores
#'
#' `escape.matrix()` computes per-cell enrichment for arbitrary gene-set
#' collections using one of four scoring back-ends and returns a dense numeric
#' matrix (cells × gene-sets).  The expression matrix is processed in
#' user-defined *chunks* (`groups`) so that memory use remains predictable;
#' each chunk is dispatched in parallel via a \pkg{BiocParallel} `BPPARAM`
#' backend.  Heavy engines (\pkg{GSVA}, \pkg{UCell}, \pkg{AUCell}) are loaded
#' lazily, keeping them in the package’s \strong{Suggests} field.
#'
#' @section Supported methods:
#' \describe{
#'   \item{`"GSVA"`}{Gene-set variation analysis (Poisson kernel).}
#'   \item{`"ssGSEA"`}{Single-sample GSEA.}
#'   \item{`"UCell"`}{Rank-based UCell scoring.}
#'   \item{`"AUCell"`}{Area-under-the-curve ranking score.}
#' }
#'
#' @param input.data A raw‐counts matrix (`genes × cells`), a
#' \link[SeuratObject]{Seurat} object, or a
#' \link[SingleCellExperiment]{SingleCellExperiment}. Gene identifiers must
#' match those in `gene.sets`.
#' @param gene.sets A named list of character vectors, the result of
#' [getGeneSets()], or the built-in data object [escape.gene.sets]. 
#' List names become column names in the result.
#' @param method Scoring algorithm (case-insensitive). One of `"GSVA"`, 
#' `"ssGSEA"`, `"UCell"`, or `"AUCell"`. Default **`"ssGSEA"`**.
#' @param groups Integer ≥ 1. Number of cells per processing chunk.
#'   Larger values reduce overhead but increase memory usage.  Default **1000**.
#' @param min.size Minimum number of genes from a set that must be detected
#' in the expression matrix for that set to be scored.  Default **5**.
#' Use `NULL` to disable filtering.
#' @param normalize  Logical. If `TRUE`, the score matrix is passed to
#' [performNormalization()] (drop-out scaling and optional log transform).
#' Default **FALSE**.
#' @param make.positive Logical. If `TRUE` *and* `normalize = TRUE`, shifts
#' every gene-set column so its global minimum is zero, facilitating
#' downstream log-ratio analyses.  Default **FALSE**.
#' @param min.expr.cells Numeric. Gene-expression filter threshold (see
#' details above). Default **0** (no gene filtering).
#' @param min.filter.by Character or `NULL`.  Column name in `meta.data`
#' (Seurat) or `colData` (SCE) defining groups within which the
#' `min.expr.cells` rule is applied.  Default **`NULL`**.
#' @param BPPARAM A \pkg{BiocParallel} parameter object describing the
#' parallel backend. 
#' @param ... Extra arguments passed verbatim to the chosen back-end
#' scoring function (`gsva()`, `ScoreSignatures_UCell()`, or
#' `AUCell_calcAUC()`).
#'
#' @return A numeric matrix with one row per cell and one column per gene set,
#' ordered as in `gene.sets`.
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso [runEscape()] to attach scores to a single-cell object;
#' [getGeneSets()] for MSigDB retrieval; [performNormalization()] for the
#' optional normalization workflow.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' pbmc <- SeuratObject::pbmc_small
#' es   <- escape.matrix(pbmc, 
#'                       gene.sets = gs,
#'                       method = "ssGSEA", 
#'                       groups = 500, 
#'                       min.size = 3)
#'
#' @export
escape.matrix <- function(input.data,
                          gene.sets        = NULL,
                          method           = "ssGSEA",
                          groups           = 1000,
                          min.size         = 5,
                          normalize        = FALSE,
                          make.positive    = FALSE,
                          min.expr.cells   = 0,
                          min.filter.by    = NULL,
                          BPPARAM          = NULL,
                          ...) {
  if(is.null(min.size)) min.size <- 0
  
  # ---- 1) resolve gene-sets & counts ----------------------------------------
  egc  <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")  # dgCMatrix
  
  if (is.null(min.filter.by)) {
    cnts <- .filter_genes(cnts, min.expr.cells)
  } else {
    # get grouping factor from object
    group.vec <- .extract_group_vector(input.data, min.filter.by)
    split.idx <- split(seq_len(ncol(cnts)), group.vec)
    
    cnts <- do.call(cbind, lapply(split.idx, function(cols) {
      sub <- cnts[, cols, drop = FALSE]
      .filter_genes(sub, min.expr.cells)
    }))
  }
  
  # ---- 2) drop undersized gene-sets -----------------------------------------
  keep <- vapply(egc, function(gs) sum(rownames(cnts) %in% gs) >= min.size,
                 logical(1))
  if (!all(keep)) {
    egc <- egc[keep]
    if (!length(egc))
      stop("No gene-sets meet the size threshold (min.size = ", min.size, ")")
  }
  
  # ---- 3) split cells into chunks -------------------------------------------
  chunks <- .split_cols(cnts, groups)
  message("escape.matrix(): processing ", length(chunks), " chunk(s)...")
  
  # ---- 4) compute enrichment in parallel ------------------------------------
  res_list <- .plapply(
    chunks,
    function(mat)
      .compute_enrichment(mat, egc, method, BPPARAM, ...),
    BPPARAM  = BPPARAM
  )
  
  # ---- 5) combine + orient (rows = cells) -----------------------------------
  all_sets <- names(egc)
  res_mat  <- do.call(cbind, lapply(res_list, function(m) {
    m <- as.matrix(m)
    m <- m[match(all_sets, rownames(m)), , drop = FALSE]
    m
  }))
  res_mat <- t(res_mat)
  colnames(res_mat) <- all_sets
  
  # ---- 6) optional dropout scaling ------------------------------------------
  if (normalize) {
    res_mat <- performNormalization(
      input.data      = input.data,
      enrichment.data = res_mat,
      assay           = NULL,
      gene.sets       = gene.sets,
      make.positive   = make.positive,
      groups          = groups
    )
    if (.is_seurat_or_sce(input.data)) {
      res_mat <- .pull.Enrich(res_mat, "escape_normalized")
    }
  }
  
  res_mat
}

#' Calculate Enrichment Scores Using Seurat or SingleCellExperiment Objects
#'
#' `runEscape()` is a convenience wrapper around [escape.matrix()] that
#' computes enrichment scores and inserts them as a new assay (default
#' `"escape"`) in a \pkg{Seurat} or \pkg{SingleCellExperiment} object.  All
#' arguments (except `new.assay.name`) map directly to their counterparts in
#' `escape.matrix()`.
#'
#' @inheritParams escape.matrix
#' @param new.assay.name Character. Name for the assay that will store the
#' enrichment matrix in the returned object. Default **"escape"**.
#'
#' @return The input single-cell object with an additional assay containing the
#' enrichment scores (`cells × gene-sets`). Matrix orientation follows
#' standard single-cell conventions (gene-sets as rows inside the assay).
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso [escape.matrix()] for the underlying computation,
#' [performNormalization()] to add normalized scores, [heatmapEnrichment()], 
#' [ridgeEnrichment()] and related plotting helpers for visualization.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' sce <- SeuratObject::pbmc_small
#' sce <- runEscape(sce, 
#'                  gene.sets = gs, 
#'                  method = "GSVA",
#'                  groups = 1000, 
#'                  min.size = 3,
#'                  new.assay.name = "escape")
#'
#' @export
runEscape <- function(input.data,
                      gene.sets,
                      method = c("ssGSEA", "GSVA", "UCell", "AUCell"),
                      groups = 1e3,
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      min.expr.cells   = 0,
                      min.filter.by    = NULL,
                      BPPARAM = NULL,
                      ...) {
    method <- match.arg(method)
    .checkSingleObject(input.data)
    esc <- escape.matrix(input.data, gene.sets, method, groups, min.size,
                         normalize, make.positive, min.expr.cells, 
                         min.filter.by, BPPARAM, ...)
    
    input.data <- .adding.Enrich(input.data, esc, new.assay.name)
    return(input.data)
}


.filter_genes <- function(m, min.expr.cells) {
  if (is.null(min.expr.cells) || identical(min.expr.cells, 0))
    return(m)                        # nothing to do
  
  ncells <- ncol(m)
  
  thr <- if (min.expr.cells < 1)
    ceiling(min.expr.cells * ncells)  # proportion → absolute
  else
    as.integer(min.expr.cells)
  
  keep <- Matrix::rowSums(m > 0) >= thr
  m[keep, , drop = FALSE]
}

# helper: pull a column from meta.data / colData no matter the object
#' @importFrom SummarizedExperiment colData
.extract_group_vector <- function(obj, col) {
  if (.is_seurat(obj))
    return(obj[[col, drop = TRUE]])
  if (.is_sce(obj))
    return(colData(obj)[[col]])
  stop("min.filter.by requires a Seurat or SingleCellExperiment object")
}

.filter_genes <- function(m, min.expr.cells) {
  if (is.null(min.expr.cells) || identical(min.expr.cells, 0))
    return(m)                        # nothing to do
  
  ncells <- ncol(m)
  
  thr <- if (min.expr.cells < 1)
    ceiling(min.expr.cells * ncells)  # proportion → absolute
  else
    as.integer(min.expr.cells)
  
  keep <- Matrix::rowSums(m > 0) >= thr
  m[keep, , drop = FALSE]
}
