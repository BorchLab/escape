#' Perform Normalization on Enrichment Data
#' 
#' This function allows users to normalize the enrichment calculations 
#' by accounting for single-cell dropout and producing positive 
#' values for downstream differential enrichment analyses. Default calculation 
#' uses will scale the enrichment values by the number of genes present from 
#' the gene set and then use a natural log transformation. A positive range
#' values is useful for several downstream analyses, like differential 
#' evaluation for log2-fold change, but will alter the original 
#' enrichment values.
#' 
#' @param sc.data Single-cell object or matrix used in the gene set enrichment calculation in 
#' \code{\link{escape.matrix}} or \code{\link{runEscape}}.
#' @param enrichment.data The enrichment results from \code{\link{escape.matrix}} 
#' or \code{\link{runEscape}} (optional)
#' @param assay Name of the assay to normalize if using a single-cell object
#' @param gene.sets The gene set library to use to extract 
#' the individual gene set information from
#' @param scale.factor A vector to use for normalizing enrichment scores per cell.
#' @param make.positive Shift enrichment values to a positive range \strong{TRUE}
#' for downstream analysis or not \strong{TRUE} (default).
#' @param groups the number of cells to calculate normalization on at once.
#' chunks matrix into groups sized chunks. Useful in case of memory issues.
#' @importFrom stringr str_replace_all
#' @importFrom SeuratObject Assays
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix colSums

#' @examples
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
#' @export
#' @return Single-cell object or matrix of normalized enrichment scores

performNormalization <- function(sc.data,
                                 enrichment.data = NULL,
                                 assay = "escape",
                                 gene.sets = NULL,
                                 make.positive = FALSE,
                                 scale.factor = NULL,
                                 groups = NULL) {
  if(!is.null(assay)) {
    if(is_seurat_object(sc.data)) {
      assay.present <- assay %in% Assays(sc.data)
    } else if (is_se_object(sc.data)) {
      assay.present <- assay %in% assays(sc.data)
    }
  } else {
    assay.present <- FALSE
  }
  
  
  if(is_seurat_or_se_object(sc.data) & !is.null(assay) & assay.present) {
    enriched <- .pull.Enrich(sc.data, assay)
  } else {
    enriched <- enrichment.data
  }
  
  if(!is.null(scale.factor) & length(scale.factor) != dim(sc.data)[2]) {
    stop("If using a vector as a scale factor, please ensure the length matches the number of cells.")
  }
  
  #Getting the gene sets that passed filters
  egc <- .GS.check(gene.sets)
  names(egc) <- str_replace_all(names(egc), "_", "-")
  egc <- egc[names(egc) %in% colnames(enriched)]
  
  #Isolating the number of genes per cell expressed
  if(is.null(groups)){
    chunks <- dim(enriched)[[1]]
  } else{
    chunks <- min(groups, dim(enriched)[[1]])
  }
  
  if (is.null(scale.factor)) {
    cnts <- .cntEval(sc.data, assay = "RNA", type = "counts")
    print("Calculating features per cell...")
    egc.sizes <- lapply(egc, function(x){
      scales<-unname(Matrix::colSums(cnts[which(rownames(cnts) %in% x),]!=0))
      scales[scales==0] <- 1
      scales
    })
    egc.sizes <- split_rows(do.call(cbind,egc.sizes), chunk.size=chunks)
    rm(cnts)
  } else{
    egc.sizes <- split_vector(scale.factor, chunk.size=chunks)
  }
  enriched <- split_rows(enriched, chunk.size=chunks)
  
  print("Normalizing enrichment scores per cell...")
  #Dividing the enrichment score by number of genes expressed

  enriched<-mapply(function(scores, scales){
    scores/scales
  }, enriched, egc.sizes, SIMPLIFY = FALSE)
  enriched <- do.call(rbind, enriched)
  if(make.positive){
    enriched <- apply(enriched, 2, function(x){
      x+max(0, -min(x))
    })
  }
  #Default Scaling using natural log
  if(is.null(scale.factor)) {
    enriched <- suppressWarnings(ifelse(enriched >= 0, 
                                 log1p(enriched + 1e-6), 
                                 -log1p(abs(enriched) + 1e-6)))
  }
  
  if(is_seurat_or_se_object(sc.data)) {
    if(is.null(assay)) {
      assay <- "escape"
    }
    sc.data <- .adding.Enrich(sc.data, enriched, paste0(assay, "_normalized"))
    return(sc.data)
  } else {
    return(enriched)
  }
}
