#' Normalize enrichment scores by expressed‑gene counts per cell
#' 
#' @description
#' Scales each enrichment value by the **number of genes from the set that are
#' expressed** in that cell (non‑zero counts). Optionally shifts results into a
#' positive range and/or applies a natural‑log transform for compatibility with
#' log‑based differential tests.
#'
#' @inheritParams escape_matrix
#' @param sc.data      Single‑cell object used to generate *raw* enrichment, or a
#'                     matrix of counts (cells × genes) when `enrichment.data`
#'                     is supplied.
#' @param enrichment.data Matrix with raw enrichment scores (cells × gene sets).
#'                        Required when `sc.data` is a plain matrix.
#' @param assay        Name of the assay to read/write inside `sc.data` when it
#'                     is a Seurat / SCE object. Default is "escape".
#' @param gene.sets    The gene‑set definitions originally used. Needed to count
#'                     expressed genes per set.
#' @param make.positive Logical; if `TRUE` shifts each column so its minimum is
#'                     zero.
#' @param scale.factor Optional numeric vector overriding gene‑count scaling
#'                     (length = #cells). Use when you want external per‑cell
#'                     normalisation factors.
#' @param groups       Chunk size (cells per block) when memory is limited.
#'
#' @example 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' pbmc_small <- performNormalization(pbmc_small, 
#'                                    assay = "escape", 
#'                                    gene.sets = GS)
#'
#' @return If `sc.data` is an object, the same object with a new assay
#'         "<assay>_normalized". Otherwise a matrix of normalised scores.
#' @export

performNormalization <- function(sc.data,
                                 enrichment.data = NULL,
                                 assay           = "escape",
                                 gene.sets       = NULL,
                                 make.positive   = FALSE,
                                 scale.factor    = NULL,
                                 groups          = NULL) {
  ## ----------------------------------------------------------------------
  ## 1. Retrieve enrichment matrix ---------------------------------------
  assay.present <- FALSE
  if (!is.null(assay) && .is_sc_object(sc.data)) {
    if (.is_seurat(sc.data)) {
      assay.present <- assay %in% SeuratObject::Assays(sc.data)
    } else if (.is_sce(sc.data) || .is_se(sc.data)) {
      assay.present <- assay %in% names(SummarizedExperiment::altExps(sc.data))
    }
  }
  
  enriched <- if (assay.present) .pull.Enrich(sc.data, assay) else enrichment.data
  if (is.null(enriched)) stop("Could not obtain enrichment matrix – please set 'assay' or supply 'enrichment.data'.")
  
  ## ----------------------------------------------------------------------
  ## 2. Validate / derive scale factors ----------------------------------
  if (!is.null(scale.factor) && length(scale.factor) != nrow(enriched))
    stop("Length of 'scale.factor' must match number of cells.")
  
  if (is.null(scale.factor)) {
    egc <- .GS.check(gene.sets)
    names(egc) <- gsub("_", "-", names(egc), fixed = TRUE)
    egc <- egc[names(egc) %in% colnames(enriched)]
    if (!length(egc)) stop("None of the supplied gene sets match enrichment columns.")
    
    ## counts matrix (genes × cells) – drop after use to save RAM
    cnts <- .cntEval(sc.data, assay = "RNA", type = "counts")
    message("Computing expressed‑gene counts per cell …")
    scale.mat <- do.call(cbind, lapply(egc, function(gs) {
      vec <- Matrix::colSums(cnts[rownames(cnts) %in% gs, , drop = FALSE] != 0)
      vec[vec == 0] <- 1L  # avoid /0
      vec
    }))
    rm(cnts)
    ## optionally split large matrices to spare memory
    chunksize <- if (is.null(groups)) nrow(enriched) else min(groups, nrow(enriched))
    sf.split  <- .split_rows(scale.mat, chunk = chunksize)
  } else {
    sf.split  <- .split_vector(scale.factor, chunk = if (is.null(groups)) length(scale.factor) else min(groups, length(scale.factor)))
  }
  
  ## ----------------------------------------------------------------------
  ## 3. Chunked normalisation --------------------------------------------
  message("Normalising enrichment scores …")
  en.split <- .split_rows(enriched, chunk = if (is.null(groups)) nrow(enriched) else min(groups, nrow(enriched)))
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
  
  ## ----------------------------------------------------------------------
  ## 6. Return ------------------------------------------------------------
  if (.is_sc_object(sc.data)) {
    .adding.Enrich(sc.data, normalized, paste0(assay %||% "escape", "_normalized"))
  } else {
    normalized
  }
}