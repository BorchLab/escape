#' Calculate Single-Cell Gene-Set Enrichment Scores
#'
#' \code{escape.matrix()} computes per-cell enrichment for arbitrary gene-set
#' collections using one of four scoring back-ends and returns a dense numeric
#' matrix (cells x gene-sets). The expression matrix is processed in
#' user-defined chunks (\code{groups}) so that memory use remains predictable;
#' each chunk is dispatched in parallel via a \pkg{BiocParallel} \code{BPPARAM}
#' backend. Heavy engines (\pkg{GSVA}, \pkg{UCell}, \pkg{AUCell}) are loaded
#' lazily, keeping them in the package's \strong{Suggests} field.
#'
#' @section Supported methods:
#' \describe{
#'   \item{\code{"GSVA"}}{Gene-set variation analysis (Poisson kernel).}
#'   \item{\code{"ssGSEA"}}{Single-sample GSEA.}
#'   \item{\code{"UCell"}}{Rank-based UCell scoring.}
#'   \item{\code{"AUCell"}}{Area-under-the-curve ranking score.}
#'   \item{\code{"PLAID"}}{Average log-intensity of set members (plaid only).}
#'   \item{\code{"singscore"}}{Rank-based singscore (plaid only).}
#'   \item{\code{"scSE"}}{Single-cell signature explorer score (plaid only).}
#' }
#'
#' @section Backends:
#' The first four methods run on escape's own engines by default
#' (\code{backend = "native"}). Setting \code{backend = "plaid"} routes them to
#' \pkg{plaid}'s \code{replaid.*} family instead, which is substantially faster
#' and lighter on memory for large objects. \code{"PLAID"}, \code{"singscore"}
#' and \code{"scSE"} have no native implementation and always use \pkg{plaid}.
#'
#' \strong{The plaid backend approximates rather than reproduces, and the gap is
#' wider than the plaid documentation suggests.} On a simulated 2000-gene,
#' 120-cell matrix, pooled Pearson correlation between the native and plaid
#' scores for the same method was roughly 0.85 (\code{ssGSEA}), 0.83
#' (\code{UCell}, \code{AUCell}) and 0.73 (\code{GSVA}); on the 230-gene
#' \code{pbmc_small} it was lower still. The scores are also on different
#' scales, not just noisier. Notably \code{replaid.ssgsea} is documented as
#' exact at \code{alpha = 0}, but compared directly against
#' \code{GSVA::gsva()} on identical input it correlated at 0.85, not 1. The
#' input assay is not the cause - these are rank-based scores, and counts
#' versus logcounts correlate at exactly 1.
#'
#' Treat \code{backend = "plaid"} as a fast screen, not as a drop-in
#' replacement. Do not mix backends within an analysis, and do not compare
#' plaid scores against previously published \pkg{escape} results. Tuning that
#' may narrow the gap for individual methods:
#' \describe{
#'   \item{\code{ssGSEA}}{\code{backend.args = list(alpha = 0)}.}
#'   \item{\code{GSVA}}{the empirical CDF row transform is approximated by a
#'     z-transform (\code{rowtf = "z"}); pass
#'     \code{backend.args = list(rowtf = "ecdf")} for the slower exact form.}
#'   \item{\code{UCell}}{\code{backend.args = list(rmax = ...)} shifts the score
#'     scale but did not change rank agreement in testing.}
#'   \item{\code{scSE}}{plaid documents a match to the original with
#'     \code{backend.args = list(removeLog2 = TRUE, scoreMean = FALSE)}.}
#' }
#'
#' Users of the plaid backend should cite Zito \emph{et al.}, \emph{Bioinformatics}
#' 2025, 41(12):btaf621 in addition to the original method paper and
#' \pkg{escape}.
#'
#' @param input.data A raw-counts matrix (genes x cells), a
#'   \link[SeuratObject]{Seurat} object, or a
#'   \link[SingleCellExperiment]{SingleCellExperiment} (including a
#'   \link[SpatialExperiment]{SpatialExperiment}). Gene identifiers must
#'   match those in \code{gene.sets}.
#' @param gene.sets A named list of character vectors, the result of
#'   \code{\link{getGeneSets}}, or the built-in data object
#'   \code{\link{escape.gene.sets}}. List names become column names in the
#'   result.
#' @param method Character. Scoring algorithm (case-insensitive). One of
#'   \code{"GSVA"}, \code{"ssGSEA"}, \code{"UCell"}, \code{"AUCell"},
#'   \code{"PLAID"}, \code{"singscore"}, or \code{"scSE"}. The last three are
#'   available only through \code{backend = "plaid"} and select it
#'   automatically. Default is \code{"ssGSEA"}.
#' @param groups Integer. Number of cells per processing chunk. Larger values
#'   reduce overhead but increase memory usage. Default is \code{1000}.
#'   Meaning depends on the backend: chunk size for the \pkg{BiocParallel} loop
#'   when \code{backend = "native"}, forwarded to \code{plaid::plaid(chunk=)}
#'   for \code{method = "PLAID"}, and ignored for the \code{replaid.*} paths
#'   (use \code{backend.args$chunk} there).
#' @param min.size Integer or \code{NULL}. Minimum number of genes from a set
#'   that must be detected in the expression matrix for that set to be scored.
#'   Default is \code{5}. Use \code{NULL} to disable filtering.
#' @param normalize Logical. If \code{TRUE}, the score matrix is passed to
#'   \code{\link{performNormalization}} (drop-out scaling and optional log
#'   transform). Default is \code{FALSE}.
#' @param make.positive Logical. If \code{TRUE} \emph{and}
#'   \code{normalize = TRUE}, shifts every gene-set column so its global
#'   minimum is zero, facilitating downstream log-ratio analyses. Default is
#'   \code{FALSE}.
#' @param min.expr.cells Numeric. Gene-expression filter threshold. Default is
#'   \code{0} (no gene filtering).
#' @param min.filter.by Character or \code{NULL}. Column name in
#'   \code{meta.data} (Seurat) or \code{colData} (SCE) defining groups within
#'   which the \code{min.expr.cells} rule is applied. Default is \code{NULL}.
#' @param BPPARAM A \pkg{BiocParallel} parameter object describing the
#'   parallel backend. Default is \code{NULL} (serial execution).
#' @param ... Extra arguments passed verbatim to the chosen native scoring
#'   function (\code{gsva()}, \code{ScoreSignatures_UCell()}, or
#'   \code{AUCell_calcAUC()}). Not forwarded when \code{backend = "plaid"} -
#'   use \code{backend.args} there.
#' @param backend Character. Scoring engine, \code{"native"} (default) or
#'   \code{"plaid"}. Ignored for methods that only exist in \pkg{plaid}.
#' @param backend.args Named list of method-specific tuning arguments passed to
#'   the underlying \pkg{plaid} function, e.g. \code{alpha}, \code{tau},
#'   \code{rowtf}, \code{aucMaxRank}, \code{rmax}, \code{nsmooth},
#'   \code{stats}, \code{chunk}, \code{removeLog2}, \code{scoreMean}. Names are
#'   validated against the target function. Note that
#'   \code{backend.args$normalize} is \pkg{plaid}'s median normalization of the
#'   scores and is unrelated to escape's \code{normalize} argument. Default is
#'   \code{list()}.
#' @param input.assay Character. Which expression matrix to score.
#'   \code{"auto"} (default) reads raw counts for the native backend and
#'   log-normalized values for \pkg{plaid}. \code{"counts"} and
#'   \code{"logcounts"} map to the right layer for both \pkg{Seurat}
#'   (\code{counts} / \code{data}) and \pkg{SummarizedExperiment}-derived
#'   objects, including \link[SpatialExperiment]{SpatialExperiment}. Any other
#'   string is taken literally.
#'
#' @return A numeric matrix with one row per cell and one column per gene set,
#'   ordered as in \code{gene.sets}.
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{runEscape}} to attach scores to a single-cell object;
#'   \code{\link{getGeneSets}} for MSigDB retrieval;
#'   \code{\link{performNormalization}} for the optional normalization workflow.
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
                          ...,
                          backend          = c("native", "plaid"),
                          backend.args     = list(),
                          input.assay      = "auto") {
  if(is.null(min.size)) min.size <- 0

  # ---- 0) resolve method / backend ------------------------------------------
  res     <- .resolve_backend(method, backend)
  method  <- res$method
  backend <- res$backend

  # ---- 1) resolve gene-sets & counts ----------------------------------------
  egc  <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA",
                   type = .resolve_input_assay(input.data, input.assay,
                                               backend))  # dgCMatrix

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
  
  # ---- 3-5) score -----------------------------------------------------------
  if (backend == "plaid") {
    ## plaid does its own chunked crossprod and its own forking. Running
    ## escape's chunk loop on top would (a) make rank- and median-based scores
    ## depend on `groups`, (b) fork inside a fork, and (c) throw away the sparse
    ## speed-up that is the whole point of the backend.
    if (...length())
      stop("Extra arguments in `...` are not forwarded when ",
           "backend = \"plaid\". Pass method-specific tuning through ",
           "`backend.args = list(...)`.", call. = FALSE)
    if (!is.null(BPPARAM) && !inherits(BPPARAM, "SerialParam"))
      message("backend = \"plaid\" manages its own parallelism; ",
              "`BPPARAM` is ignored.")
    if (normalize && !identical(backend.args$normalize, FALSE))
      message("normalize = TRUE stacks escape's drop-out scaling on top of ",
              "plaid's median normalization. Set ",
              "backend.args = list(normalize = FALSE) to use escape's only.")

    .require_plaid()
    message("escape.matrix(): scoring ", ncol(cnts),
            " cells with the plaid backend...")
    res_mat <- .plaid_dispatch(cnts, egc, method,
                               min.size     = min.size,
                               backend.args = backend.args,
                               groups       = groups)
  } else {
    chunks <- .split_cols(cnts, groups)
    message("escape.matrix(): processing ", length(chunks), " chunk(s)...")

    res_list <- .plapply(
      chunks,
      function(mat)
        .compute_enrichment(mat, egc, method, BPPARAM, ...),
      BPPARAM  = BPPARAM
    )

    ## combine + orient (rows = cells)
    all_sets <- names(egc)
    res_mat  <- do.call(cbind, lapply(res_list, function(m) {
      m <- as.matrix(m)
      m <- m[match(all_sets, rownames(m)), , drop = FALSE]
      m
    }))
    res_mat <- t(res_mat)
    colnames(res_mat) <- all_sets
  }

  # ---- 6) optional dropout scaling ------------------------------------------
  if (normalize) {
    ## assay = NULL keeps this on the matrix path for every input class - the
    ## previous round trip through .adding.Enrich()/.pull.Enrich() is what broke
    ## SingleCellExperiment and SpatialExperiment input.
    res_mat <- performNormalization(
      input.data      = input.data,
      enrichment.data = res_mat,
      assay           = NULL,
      gene.sets       = egc,
      make.positive   = make.positive,
      groups          = groups
    )
  }

  .stamp_backend(res_mat, method, backend)
}

#' Calculate Enrichment Scores Using Seurat or SingleCellExperiment Objects
#'
#' \code{runEscape()} is a convenience wrapper around \code{\link{escape.matrix}}
#' that computes enrichment scores and inserts them as a new assay (default
#' \code{"escape"}) in a \pkg{Seurat} or \pkg{SingleCellExperiment} object. All
#' arguments (except \code{new.assay.name}) map directly to their counterparts
#' in \code{escape.matrix()}.
#'
#' @inheritParams escape.matrix
#' @param new.assay.name Character. Name for the assay that will store the
#'   enrichment matrix in the returned object. Default is \code{"escape"}.
#'
#' @return The input single-cell object with an additional assay containing the
#'   enrichment scores (cells x gene-sets). Matrix orientation follows standard
#'   single-cell conventions (gene-sets as rows inside the assay).
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{escape.matrix}} for the underlying computation;
#'   \code{\link{performNormalization}} to add normalized scores;
#'   \code{\link{heatmapEnrichment}}, \code{\link{ridgeEnrichment}}, and
#'   related plotting helpers for visualization.
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
                      method = c("ssGSEA", "GSVA", "UCell", "AUCell",
                                 "PLAID", "singscore", "scSE"),
                      groups = 1e3,
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      min.expr.cells   = 0,
                      min.filter.by    = NULL,
                      BPPARAM = NULL,
                      ...,
                      backend      = c("native", "plaid"),
                      backend.args = list(),
                      input.assay  = "auto") {
    method  <- match.arg(method)
    backend <- match.arg(backend)
    .checkSingleObject(input.data)

    ## named, not positional - inserting an argument above must never silently
    ## shift what the callee receives
    esc <- escape.matrix(input.data     = input.data,
                         gene.sets      = gene.sets,
                         method         = method,
                         groups         = groups,
                         min.size       = min.size,
                         normalize      = normalize,
                         make.positive  = make.positive,
                         min.expr.cells = min.expr.cells,
                         min.filter.by  = min.filter.by,
                         BPPARAM        = BPPARAM,
                         ...,
                         backend        = backend,
                         backend.args   = backend.args,
                         input.assay    = input.assay)

    prov       <- attr(esc, "escape.backend")
    input.data <- .adding.Enrich(input.data, esc, new.assay.name)
    input.data <- .record_backend(input.data, new.assay.name, prov)
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
