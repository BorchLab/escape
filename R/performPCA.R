#' Perform Principal Component Analysis on Enrichment Data
#' 
#' This function allows users to calculate the principal components 
#' for the gene set enrichment values. For single-cell data, the PCA
#' will be stored with the dimensional reductions. If a matrix is used
#' as input, the output is a list for further plotting. Alternatively,
#' users can use functions for PCA calculations based on their desired
#' workflow in lieu of using \code{\link{performPCA}}, but will not be 
#' compatible with downstream \code{\link{pcaEnrichment}} visualization.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single‑cell
#' object previously processed by \code{\link{runEscape}}.
#' @param assay Name of the assay holding enrichment scores when
#' `input.data` is a single‑cell object. Ignored otherwise.
#' @param scale Logical; if `TRUE` standardises each gene-set column
#' before PCA.
#' @param n.dim Integer >= 1 or vector; the **largest** value sets the
#'   number of principal components to compute / keep.
#' @param reduction.name,reduction.key  Names used when writing back to a
#'   Seurat / SCE object.
#'   
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs,
#'             min.size = NULL)           
#'                         
#' pbmc <- performPCA(pbmc, 
#'                    assay = "escape")
#'
#' @return *If* `input.data` is a single-cell object, the same object with a
#'   new dimensional-reduction slot.  *Otherwise* a list with  
#'   `PCA`, `eigen_values`, `contribution`, and `rotation`.
#' @export
performPCA <- function(input.data,
                       assay           = "escape",
                       scale           = TRUE,
                       n.dim           = 10,
                       reduction.name  = "escape.PCA",
                       reduction.key   = "escPC_") {
  
  ## ------------ 1  Get enrichment matrix ------------------------------------
  if (.is_seurat_or_sce(input.data)) {
    mat <- .pull.Enrich(input.data, assay)
  } else if (is.matrix(input.data) || is.data.frame(input.data)) {
    mat <- as.matrix(input.data)
  } else {
    stop("`input.data` must be a matrix/data.frame or a Seurat/SCE object.")
  }
  if (!is.numeric(mat)) stop("Enrichment matrix must be numeric.")
  
  ## ------------ 2  Choose PCA backend ---------------------------------------
  ndim <- max(as.integer(n.dim))
  use_irlba <- requireNamespace("irlba", quietly = TRUE) &&
    min(dim(mat)) > 50                     # heuristic
  
  pca_obj <- if (use_irlba) {
    irlba::prcomp_irlba(mat, n = ndim, center = TRUE, scale. = scale)
  } else {
    stats::prcomp(mat, rank. = ndim, center = TRUE, scale. = scale)
  }
  
  ## ------------ 3  Post-process ---------------------------------------------
  eig  <- pca_obj$sdev ^ 2
  pct  <- round(eig / sum(eig) * 100, 1)
  colnames(pca_obj$x) <- paste0(reduction.key, seq_len(ncol(pca_obj$x)))
  
  misc <- list(eigen_values = eig,
               contribution = pct,
               rotation      = pca_obj$rotation)
  
  ## ------------ 4  Return / write-back --------------------------------------
  if (.is_seurat_or_sce(input.data)) {
    
    if (.is_seurat(input.data)) {
      if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' is required to write PCA results into a Seurat object.")
      }
      
      input.data[[reduction.name]] <- SeuratObject::CreateDimReducObject(
        embeddings = pca_obj$x,
        loadings   = pca_obj$rotation,
        stdev      = pca_obj$sdev,
        key        = reduction.key,
        misc       = misc, 
        assay      = assay
      )
      
    } else if (.is_sce(input.data)) {
      if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' is required to write PCA results into a SingleCellExperiment object.")
      }
      
      SingleCellExperiment::reducedDim(input.data, reduction.name) <- pca_obj$x
      input.data@metadata <- c(input.data@metadata, misc)
      
    } 
    return(input.data)
    
  } else {
    list(PCA          = pca_obj$x,
         eigen_values = eig,
         contribution = pct,
         rotation     = pca_obj$rotation)
  }
}
