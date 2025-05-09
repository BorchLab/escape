#' Calculate gene set enrichment scores 
#'
#' The function processes the expression matrix in chunks (size controlled by
#' \code{groups}) so memory usage is predictable. Chunks are distributed across
#' the parallel backend defined by \pkg{BiocParallel}. Heavy scoring engines
#' (\pkg{GSVA}, \pkg{UCell}, \pkg{AUCell}) are loaded lazily, so they can live
#' in the package's \strong{Suggests} field.
#'
#' @section Supported methods:
#' \describe{
#'   \item{\code{"GSVA"}}{Gene‑set variation analysis (Poisson kernel).}
#'   \item{\code{"ssGSEA"}}{Single‑sample gene‑set enrichment.}
#'   \item{\code{"UCell"}}{Rank‑based UCell scoring.}
#'   \item{\code{"AUCell"}}{Area‑under‑the‑curve gene‑ranking scoring.}
#' }
#'
#' @param input.data A count matrix (genes × cells), a \pkg{SeuratObject}, or a
#'   \pkg{SingleCellExperiment}. Gene names must match those used in
#'   \code{gene.sets}.
#' @param gene.sets A named list of character vectors, the output of
#'   \code{\link{getGeneSets}}, or the built‑in \code{\link{escape.gene.sets}}. 
#'   List names become column names in the returned matrix.
#' @param method Scoring algorithm to use. One of \code{"GSVA"}, \code{"ssGSEA"},
#'   \code{"UCell"}, or \code{"AUCell"} (case‑insensitive). Default
#'   \code{"ssGSEA"}.
#' @param groups Integer. Number of cells to process per chunk. Affects memory
#'   use and parallel granularity. Default \code{1000}.
#' @param min.size Minimum number of genes from a set that must be present in
#'   the expression matrix for the set to be scored. Default \code{5}. Set to
#'   \code{NULL} to disable filtering.
#' @param normalize Logical; if \code{TRUE} the score matrix is passed to
#'   \code{\link{performNormalization}} for dropout scaling.
#' @param make.positive Logical; if \code{TRUE} (and \code{normalize = TRUE})
#'   shifts the normalized scores so that the minimum value across all cells is
#'   zero.
#' @param BPPARAM A \pkg{BiocParallel} parameter object describing the parallel
#'   backend. Defaults to \code{BiocParallel::SerialParam()} for serial
#'   execution.
#' @param ... Additional arguments forwarded to the chosen back‑end scoring
#'   function.
#'
#' @return A numeric matrix of enrichment scores with cells in rows and gene
#'   sets in columns (ordered as in \code{gene.sets}).
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{runEscape}} to attach the matrix to a single‑cell
#'   object; \code{\link{getGeneSets}} for convenient gene‑set retrieval.
#'
#' @examples
#' gs <- list(B = c("MS4A1", "CD79B", "CD79A"),
#'            T = c("CD3E", "CD3D", "CD3G"))
#' pbmc <- SeuratObject::pbmc_small
#' es <- escape_matrix(pbmc, gene.sets = gs, min.size = 3, groups = 500)
#'
#' @importFrom BiocParallel bplapply SerialParam
#' @export
escape.matrix <- function(input.data, 
                          gene.sets = NULL, 
                          method = "ssGSEA", 
                          groups = 1000, 
                          min.size = 5,
                          normalize = FALSE,
                          make.positive = FALSE,
                          BPPARAM = SerialParam(),
                          ...) {
  # 1) Resolve gene‑sets & counts 
  egc  <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  
  # 2) Filter gene‑sets shorter than min.size ---------------------------------
  keep <- vapply(egc, function(gs) sum(rownames(cnts) %in% gs) >= min.size,
                 logical(1))
  if (!all(keep)) {
    egc <- egc[keep]
    if (!length(egc))
      stop("No gene‑sets meet the size threshold (min.size=", min.size, ")")
  }
  
  # 3) Split cells 
  chunks <- .split_cols(cnts, groups)
  message("escape_matrix(): processing ", length(chunks), " chunk(s)…")
  
  # 4) Compute enrichment per chunk in parallel 
  res_list <- BiocParallel::bplapply(chunks, function(mat) {
    .compute_enrichment(mat, egc, method, BPPARAM, ...)
  }, BPPARAM = BPPARAM)
  
  # 5) Combine, transpose so rows=cells, cols=gene‑sets ------------------------
  all_sets <- names(egc)
  res_mat <- do.call(cbind, lapply(res_list, function(m) {
    m <- as.matrix(m)
    m <- m[match(all_sets, rownames(m)), , drop = FALSE]
    m
  }))
  res_mat <- t(res_mat)
  colnames(res_mat) <- all_sets
  
  # 6) Optional dropout normalization ----------------------------------------
  if (normalize) {
    res_mat <- performNormalization(sc.data = input.data,
                                    enrichment.data = res_mat,
                                    assay = NULL,
                                    gene.sets = gene.sets,
                                    make.positive = make.positive,
                                    groups = groups)
  }
  
  res_mat
}


#' Enrichment calculation for single-cell workflows
#'
#' Run the escape-based gene-set enrichment calculation with 
#' Seurat or SingleCellExperiment pipelines
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'
#' @param input.data The count matrix, Seurat, or Single-Cell Experiment object.
#' @param gene.sets Gene sets can be a list, output from 
#' \code{\link{getGeneSets}}, or the built-in gene sets 
#' in the escape package \code{\link{escape.gene.sets}}.
#' @param method Select the method to calculate enrichment, \strong{AUCell},
#' \strong{GSVA}, \strong{ssGSEA} or \strong{UCell}.
#' @param groups The number of cells to separate the enrichment calculation.
#' @param min.size Minimum number of gene necessary to perform the enrichment
#' calculation
#' @param normalize Whether to divide the enrichment score by the number 
#' of genes \strong{TRUE} or report unnormalized \strong{FALSE}.
#' @param make.positive During normalization shift enrichment values to a 
#' positive range \strong{TRUE} for downstream analysis or not 
#' \strong{TRUE} (default). Will only be applied if \strong{normalize = TRUE}.
#' @param new.assay.name The new name of the assay to append to 
#' the single-cell object containing the enrichment scores.
#' @param BPPARAM A BiocParallel::bpparam() object that for parallelization. 
#' @param ... pass arguments to AUCell GSVA, ssGSEA or UCell call
#' @export
#' @return Seurat or Single-Cell Experiment object with escape enrichment scores
#' in the assay slot. 

runEscape <- function(input.data,
                      gene.sets,
                      method = c("ssGSEA", "GSVA", "UCell", "AUCell"),
                      groups = 1e3,
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      BPPARAM = BiocParallel::SerialParam(),
                      ...) {
    method <- match.arg(method)
    .checkSingleObject(input.data)
    esc <- escape_matrix(input.data, gene.sets, method, groups, min.size,
                         normalize, make.positive, BPPARAM, ...)
    .adding.Enrich(input.data, esc, new.assay.name)
    
    input.data <- .adding.Enrich(input.data, enrichment, new.assay.name)
    return(input.data)
}

.gsva.setup <- function(data, egc) {
  params.to.use <- gsvaParam(exprData = data,
                             geneSets = egc,
                             kcdf = "Poisson")
  return(params.to.use)
}

.ssGSEA.setup <- function(data, egc) {
  params.to.use <- ssgseaParam(exprData = data,
                               geneSets = egc,
                               normalize = FALSE)
  return(params.to.use)
}

