#' Perform Normalization on Enrichment Data
#' 
#' @description
#' Scales each enrichment value by the **number of genes from the set that are
#' expressed** in that cell (non‑zero counts). Optionally shifts results into a
#' positive range and/or applies a natural‑log transform for compatibility with
#' log‑based differential tests.
#'
#' @param input.data raw‐counts matrix (`genes × cells`), a 
#' \link[SeuratObject]{Seurat} object, or a 
#' \link[SingleCellExperiment]{SingleCellExperiment}. Gene identifiers must
#' match those in `gene.sets`.
#' @param enrichment.data Output of \code{\link{escape.matrix}} or a single‑cell
#' object previously processed by \code{\link{runEscape}}.
#' @param assay Name of the assay holding enrichment scores when
#' `input.data` is a single‑cell object. Ignored otherwise.
#' @param gene.sets A named list of character vectors, the result of
#' [getGeneSets()], or the built-in data object
#' [escape.gene.sets].  List names become column names in the result.
#' @param make.positive Logical; if `TRUE` shifts each column so its minimum is
#' zero.
#' @param scale.factor Optional numeric vector overriding gene‑count scaling
#' (length = #cells). Use when you want external per‑cell normalization factors.
#' @param groups Integer ≥ 1.  Number of cells per processing chunk.
#' Larger values reduce overhead but increase memory usage.  Default **1000**.
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
#' @return If `input.data` is an object, the same object with a new assay
#'         "<assay>_normalized". Otherwise a matrix of normalized scores.
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
      if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        assay.present <- assay %in% names(SummarizedExperiment::altExps(input.data))
      } else {
        warning("SummarizedExperiment package is required but not installed.")
      }
    }
  }
  

  
  enriched <- if (assay.present) .pull.Enrich(input.data, assay) else enrichment.data
  if (is.null(enriched)) {
    stop("Could not obtain enrichment matrix, please set `assay` or supply `enrichment.data`.")
  }
  
  ## 2. Validate / derive scale factors ----------------------------------
  if (!is.null(scale.factor) && length(scale.factor) != nrow(enriched))
    stop("Length of 'scale.factor' must match number of cells.")
  
  if (is.null(scale.factor)) {
    egc <- .GS.check(gene.sets)
    names(egc) <- gsub("_", "-", names(egc), fixed = TRUE)
    egc <- egc[names(egc) %in% colnames(enriched)]
    if (!length(egc)) stop("None of the supplied gene sets match enrichment columns.")
    
    ## counts matrix (genes × cells) – drop after use to save RAM
    cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
    message("Computing expressed-gene counts per cell...")
    scale.mat <- do.call(cbind, lapply(egc, function(gs) {
      vec <- Matrix::colSums(cnts[rownames(cnts) %in% gs, , drop = FALSE] != 0)
      vec[vec == 0] <- 1L  # avoid /0
      vec
    }))
    rm(cnts)
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
  if (.is_seurat_or_sce(input.data)) {
    .adding.Enrich(input.data, normalized, paste0(assay %||% "escape", "_normalized"))
  } else {
    normalized
  }
}