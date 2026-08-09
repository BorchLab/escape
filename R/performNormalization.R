#' Perform Normalization on Enrichment Data
#'
#' Scales each enrichment value by the \strong{number of genes from the set
#' that are expressed} in that cell (non-zero counts). Optionally shifts
#' results into a positive range and/or applies a natural-log transform for
#' compatibility with log-based differential tests.
#'
#' @param input.data A raw-counts matrix (genes x cells), a
#'   \link[SeuratObject]{Seurat} object, or a
#'   \link[SingleCellExperiment]{SingleCellExperiment} (including a
#'   \link[SpatialExperiment]{SpatialExperiment}). Gene identifiers must
#'   match those in \code{gene.sets}.
#' @param enrichment.data Matrix. Output of \code{\link{escape.matrix}} or
#'   \code{NULL} if enrichment scores are already stored in \code{input.data}.
#'   When supplied it takes precedence over anything stored in
#'   \code{input.data}.
#' @param assay Character. Name of the assay holding enrichment scores when
#'   \code{input.data} is a single-cell object. Default is \code{"escape"}.
#'   Ignored when \code{input.data} is a matrix. Set to \code{NULL} to return
#'   the normalized matrix rather than attaching it to the object.
#' @param gene.sets A named list of character vectors, the result of
#'   \code{\link{getGeneSets}}, or the built-in data object
#'   \code{\link{escape.gene.sets}}. List names must match column names in the
#'   enrichment matrix.
#' @param make.positive Logical. If \code{TRUE}, shifts each column so its
#'   minimum is zero. Default is \code{FALSE}.
#' @param scale.factor Numeric vector or \code{NULL}. Optional per-cell scaling
#'   factors (length = number of cells). Use when you want external per-cell
#'   normalization factors. Default is \code{NULL} (compute from gene counts).
#' @param groups Integer or \code{NULL}. Number of cells per processing chunk.
#'   Larger values reduce overhead but increase memory usage. Default is
#'   \code{NULL} (process all cells at once).
#'
#' @examples 
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs,
#'             min.size = NULL)
#'                         
#' pbmc <- performNormalization(pbmc, 
#'                              assay = "escape", 
#'                              gene.sets = gs)
#'
#' @section Which expression values are used:
#' The per-cell scale factor is the number of genes from each set with a
#' \strong{non-zero raw count}, so this function always reads the \code{counts}
#' assay regardless of any \code{input.assay} used when the scores were
#' computed. Detection is identical in count and log space, so this is
#' deliberate rather than an oversight.
#'
#' @return If `input.data` is an object and `assay` is not `NULL`, the same
#'         object with a new assay "<assay>_normalized". Otherwise a matrix of
#'         normalized scores.
#' @export

performNormalization <- function(input.data,
                                 enrichment.data = NULL,
                                 assay           = "escape",
                                 gene.sets       = NULL,
                                 make.positive   = FALSE,
                                 scale.factor    = NULL,
                                 groups          = NULL) {

  ## 1. Retrieve enrichment matrix ---------------------------------------
  assay.present <- FALSE
  if (!is.null(assay) && .is_seurat_or_sce(input.data)) {
    if (.is_seurat(input.data)) {
      if (requireNamespace("SeuratObject", quietly = TRUE)) {
        assay.present <- assay %in% SeuratObject::Assays(input.data)
      } else {
        warning("SeuratObject package is required but not installed.")
      }
    } else if (.is_sce(input.data)) {
      if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        assay.present <- assay %in% SingleCellExperiment::altExpNames(input.data)
      } else {
        warning("SingleCellExperiment package is required but not installed.")
      }
    }
  }

  ## Supplied scores always win - never re-pull something the caller handed in.
  if (!is.null(enrichment.data)) {
    if (assay.present)
      warning("Both `enrichment.data` and assay '", assay, "' are available; ",
              "using `enrichment.data`.", call. = FALSE)
    enriched <- enrichment.data
  } else if (assay.present) {
    enriched <- .pull.Enrich(input.data, assay)
  } else {
    ## nothing to work with - say what was asked for and what exists
    if (!is.null(assay) && .is_seurat_or_sce(input.data)) {
      avail <- if (.is_seurat(input.data)) SeuratObject::Assays(input.data)
               else SingleCellExperiment::altExpNames(input.data)
      .stop_missing_assay(assay, avail, "enrichment assay",
                          "`input.data`, and `enrichment.data` was not supplied")
    }
    stop("Could not obtain enrichment matrix, please set `assay` or supply `enrichment.data`.")
  }

  ## 2. Validate / derive scale factors ----------------------------------
  if (!is.null(scale.factor) && length(scale.factor) != nrow(enriched))
    stop("Length of 'scale.factor' must match number of cells.")

  if (is.null(scale.factor)) {
    ## one gene set per enrichment column, in column order
    egc <- .match_sets_to_cols(.GS.check(gene.sets), colnames(enriched))

    ## counts matrix (genes x cells) - drop after use to save RAM
    cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
    message("Computing expressed-gene counts per cell...")
    scale.mat <- do.call(cbind, lapply(egc, function(gs) {
      vec <- Matrix::colSums(cnts[rownames(cnts) %in% gs, , drop = FALSE] != 0)
      vec[vec == 0] <- 1L  # avoid /0
      vec
    }))
    colnames(scale.mat) <- names(egc)
    rm(cnts)

    ## alignment is guaranteed by .match_sets_to_cols(); assert it anyway so a
    ## future refactor cannot silently divide the wrong column
    stopifnot(identical(colnames(scale.mat), colnames(enriched)))
    ## optionally split large matrices to spare memory
    chunksize <- if (is.null(groups)) nrow(enriched) else min(groups, nrow(enriched))
    sf.split  <- .split_rows(scale.mat, chunk.size  = chunksize)
  } else {
    sf.split  <- .split_vector(scale.factor, chunk.size = if (is.null(groups)) length(scale.factor) else min(groups, length(scale.factor)))
  }
  
  ## 3. Chunked normalization --------------------------------------------
  message("Normalizing enrichment scores...")
  en.split <- .split_rows(enriched, chunk.size = if (is.null(groups)) nrow(enriched) else min(groups, nrow(enriched)))
  norm.lst <- Map(function(sco, fac) sco / fac, en.split, sf.split)
  normalized <- do.call(rbind, norm.lst)
  
  ## 4. Optional positive shift ------------------------------------------
  if (make.positive) {
    shift <- pmax(0, -apply(normalized, 2L, min))
    normalized <- sweep(normalized, 2L, shift, `+`)
  }
  
  ## 5. Log transform (only when scale.factor derived internally) ---------
  if (is.null(scale.factor)) {
    neg <- normalized < 0
    normalized[!neg] <- log1p(normalized[!neg] + 1e-6)
    normalized[neg]  <- -log1p(abs(normalized[neg]) + 1e-6)
  }
  
  ## 6. Return ------------------------------------------------------------
  ## `assay = NULL` means "hand the matrix back" - that is how escape.matrix()
  ## calls this, and it avoids attaching an altExp only to pull it straight off.
  if (.is_seurat_or_sce(input.data) && !is.null(assay)) {
    return(.adding.Enrich(input.data, normalized, paste0(assay, "_normalized")))
  }
  normalized
}
