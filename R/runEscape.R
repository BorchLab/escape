#' Calculate gene set enrichment scores 
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and output the enrichment scores as a matrix.  
#'
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
#' @param BPPARAM A BiocParallel::bpparam() object that for parallelization. 
#' @param ... pass arguments to AUCell GSVA, ssGSEA, or UCell call
#'
#' @importFrom GSVA gsva gsvaParam ssgseaParam
#' @importFrom GSEABase GeneSetCollection 
#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel SerialParam MulticoreParam BatchtoolsParam SerialParam
#'
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' ES <- escape.matrix(pbmc_small, 
#'                     gene.sets = GS, 
#'                     min.size = NULL)
#' 
#' @export
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{getGeneSets}} to collect gene sets.
#' @return matrix of enrichment scores
escape.matrix <- function(input.data, 
                          gene.sets = NULL, 
                          method = "ssGSEA", 
                          groups = 1000, 
                          min.size = 5,
                          normalize = FALSE,
                          make.positive = FALSE,
                          BPPARAM = SerialParam(),
                          ...) {
    egc <- .GS.check(gene.sets)
    cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
    egc.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
    
    # Filtering gene sets
    if (!is.null(min.size)) {
      remove <- which(egc.size < min.size | egc.size == 0)
    } else {
      remove <- which(egc.size == 0)
    }
    if (length(remove) > 0) {
      egc <- egc[-remove]
      egc.size <- egc.size[-remove]
      
      if (!is.null(min.size) && length(egc) == 0) {
        stop("No gene sets passed the minimum length - please reconsider the 'min.size' parameter")
      }
    }
    
    scores <- list()
    splits <- seq(1, ncol(cnts), by=groups)
    print(paste('Using sets of', groups, 'cells. Running', 
                length(splits), 'times.'))
    split.data <- .split_data.matrix(matrix=cnts, chunk.size=groups)
    
    all_gene_sets <- names(egc) # Collect all gene set names
    
    # Running enrichment calculation
    for (i in seq_along(splits)) {
      if (method == "GSVA") {
        parameters <- .gsva.setup(split.data[[i]], egc)
      } else if (method == "ssGSEA") {
        parameters <- .ssGSEA.setup(split.data[[i]], egc)
      }
      if (method %in% c("ssGSEA", "GSVA")) {
        a <- suppressWarnings(gsva(param = parameters, 
                                   verbose = FALSE,
                                   BPPARAM = BPPARAM,
                                   ...))
      } else if (method == "UCell") {
        a <- t(suppressWarnings(
          ScoreSignatures_UCell(matrix = split.data[[i]], 
                                features = egc,
                                name = NULL,
                                BPPARAM = BPPARAM,
                                ...)))
      } else if (method == "AUCell") {
        rankings <- AUCell_buildRankings(split.data[[i]],
                                         plotStats = FALSE,
                                         verbose = FALSE)
        a <- assay(AUCell_calcAUC(geneSets = egc,
                                  rankings,
                                  normAUC = TRUE,
                                  aucMaxRank = ceiling(0.2 * nrow(split.data[[i]])),
                                  verbose = FALSE,
                                  ...))
      }
      
      # Ensure consistent row names (all_gene_sets) across splits
      a <- as.data.frame(a)
      a <- a[match(all_gene_sets, rownames(a), nomatch = NA), , drop = FALSE]
      scores[[i]] <- a
    }
    scores <- do.call(cbind, scores)
    output <- t(as.matrix(scores))
    
    # Normalize based on dropout
    if(normalize) {
      output <- performNormalization(sc.data = input.data,
                                     enrichment.data = output,
                                     assay = NULL,
                                     gene.sets = gene.sets,
                                     make.positive = make.positive,
                                     groups = groups)
    }
    return(output)
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
                      gene.sets = NULL, 
                      method = "ssGSEA", 
                      groups = 1000, 
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      BPPARAM = SerialParam(),
                      ...) {
  .checkSingleObject(input.data)
   enrichment <- escape.matrix(input.data = input.data,
                              gene.sets = gene.sets,
                              method = method,
                              groups = groups,
                              min.size = min.size,
                              BPPARAM = BPPARAM)
  
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
