# -----------------------------------------------------------------------------
#  FAST NEGATION OPERATOR ------------------------------------------------------
# -----------------------------------------------------------------------------
`%!in%` <- Negate(`%in%`)     

# -----------------------------------------------------------------------------
#  CLASS HELPERS ---------------------------------------------------------------
# -----------------------------------------------------------------------------
.is_seurat      <- function(x) inherits(x, "Seurat")
.is_sce         <- function(x) inherits(x, "SummarizedExperiment")
.is_seurat_or_sce <- function(x) .is_seurat(x) || .is_sce(x)

.checkSingleObject <- function(obj) {
  if (!.is_seurat_or_sce(obj))
    stop("Expecting a Seurat or SummarizedExperiment object")
}

# -----------------------------------------------------------------------------
#  ORDERING UTILITY (base R implementation) -----------------------------------
# -----------------------------------------------------------------------------
.orderFunction <- function(dat, order.by, group.by) {
  if (!(order.by %in% c("mean", "group.by")))
    stop("order.by must be 'mean' or 'group.by'")
  
  if (order.by == "mean") {
    means <- tapply(dat[[1]], dat[[group.by]], mean, simplify = TRUE)
    lev   <- names(sort(means, decreasing = TRUE))
    dat[[group.by]] <- factor(dat[[group.by]], levels = lev)
  } else {                              # natural sort of group labels
    if (requireNamespace("stringr", quietly = TRUE)) {
      lev <- stringr::str_sort(unique(dat[[group.by]]), numeric = TRUE)
    } else {
      lev <- sort(unique(dat[[group.by]]), method = "radix")
    }
    dat[[group.by]] <- factor(dat[[group.by]], levels = lev)
  }
  dat
}

# -----------------------------------------------------------------------------
#  DATA.frame BUILDERS ---------------------------------------------------------
# -----------------------------------------------------------------------------
.makeDFfromSCO <- function(input.data, assay = "escape", gene.set = NULL,
                           group.by = NULL, split.by = NULL, facet.by = NULL) {
  if (is.null(assay))
    stop("Please provide assay name")
  cols <- unique(c(group.by, split.by, facet.by))
  cnts <- .cntEval(input.data, assay = assay, type = "data")
  
  if (length(gene.set) == 1 && gene.set == "all")
    gene.set <- rownames(cnts)
  
  meta <- .grabMeta(input.data)
  meta <- meta[, cols, drop = FALSE]
  
  if (length(gene.set) == 1) {
    df <- cbind(value = cnts[gene.set, ], meta)
    colnames(df)[1] <- gene.set
  } else {
    df <- cbind(Matrix::t(cnts[gene.set, , drop = FALSE]), meta)
  }
  df
}

.prepData <- function(input.data, assay, gene.set, group.by, split.by, facet.by) {
  if (.is_seurat_or_sce(input.data)) {
    df <- .makeDFfromSCO(input.data, assay, gene.set, group.by, split.by, facet.by)
    if (identical(gene.set, "all")) {
      gene.set <- setdiff(colnames(df), c(group.by, split.by, facet.by))
    }
  } else {                               # assume plain data.frame / matrix
    if (identical(gene.set, "all"))
      gene.set <- setdiff(colnames(input.data), c(group.by, split.by, facet.by))
    df <- input.data[, c(gene.set, group.by, split.by, facet.by), drop = FALSE]
  }
  colnames(df) <- c(gene.set, group.by, split.by, facet.by)
  df
}

# -----------------------------------------------------------------------------
#  COLOUR SCALES (ggplot helper; tidy‑agnostic) --------------------------------
# -----------------------------------------------------------------------------
.colorizer <- function(palette = "inferno", n = NULL) {
  grDevices::hcl.colors(n = n, palette = palette, fixup = TRUE)
}

#' @importFrom stats setNames
.colorby <- function(enriched,
                     plot,
                     color.by,
                     palette,
                     type = c("fill", "color")) {
  
  type <- match.arg(type)
  
  vec    <- enriched[[color.by]]
  is_num <- is.numeric(vec)
  
  ## pick scale constructors --------------------------------------------------
  scale_discrete <- switch(type,
                           fill  = ggplot2::scale_fill_manual,
                           color = ggplot2::scale_color_manual)
  
  scale_gradient <- switch(type,
                           fill  = ggplot2::scale_fill_gradientn,
                           color = ggplot2::scale_color_gradientn)
  
  ## build scale + legend ------------------------------------------------------
  if (is_num) {
    plot <- plot +
      scale_gradient(colors = .colorizer(palette, 11)) +
      do.call(ggplot2::labs, setNames(list(color.by), type))
  } else {
    lev <- if (requireNamespace("stringr", quietly = TRUE)) {
      stringr::str_sort(unique(vec), numeric = TRUE)
    } else {
      unique(vec)
    }
    
    pal        <- .colorizer(palette, length(lev))
    names(pal) <- lev
    
    plot <- plot +
      scale_discrete(values = pal) +
      do.call(ggplot2::labs, setNames(list(color.by), type))
  }
  
  return(plot)
}

# -----------------------------------------------------------------------------
#  MATRIX / VECTOR SPLITTERS ---------------------------------------------------
# -----------------------------------------------------------------------------
.split_cols <- function(mat, chunk) {
  if (ncol(mat) <= chunk) return(list(mat))
  idx <- split(seq_len(ncol(mat)), ceiling(seq_len(ncol(mat)) / chunk))
  lapply(idx, function(i) mat[, i, drop = FALSE])
}

.split_rows <- function(mat, chunk.size = 1000) {
  if (is.vector(mat)) mat <- matrix(mat, ncol = 1)
  idx <- split(seq_len(nrow(mat)), ceiling(seq_len(nrow(mat)) / chunk.size))
  lapply(idx, function(i) mat[i, , drop = FALSE])
}

.split_vector <- function(vec, chunk.size = 1000) {
  split(vec, ceiling(seq_along(vec) / chunk.size))
}

# -----------------------------------------------------------------------------
#  EXPRESSION MATRIX EXTRACTOR -------------------------------------------------
# -----------------------------------------------------------------------------
.cntEval <- function(obj, assay = "RNA", type = "counts") {
  if (.is_seurat(obj)) {
    # Use generic accessor if available
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      suppressWarnings(
        cnts <- SeuratObject::GetAssayData(obj, assay = assay, slot = type)
      )
    } else {
      cnts <- obj@assays[[assay]][[type]]
    }
    
  } else if (.is_sce(obj)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      pos <- if (assay == "RNA") "counts" else assay
      
      cnts <- if (assay == "RNA") {
        SummarizedExperiment::assay(obj, pos)
      } else {
        SummarizedExperiment::assay(SingleCellExperiment::altExp(obj, pos))
      }
    } else {
      stop("SummarizedExperiment and SingleCellExperiment packages are required but not installed.")
    }
  } else {
    cnts <- obj
  }
  cnts[MatrixGenerics::rowSums2(cnts) != 0, , drop = FALSE]
}


# -----------------------------------------------------------------------------
#  ATTACH / PULL ENRICHMENT MATRICES ------------------------------------------
# -----------------------------------------------------------------------------
.adding.Enrich <- function(sc, enrichment, name) {
  if (.is_seurat(sc)) {
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      major <- as.numeric(substr(sc@version, 1, 1))
      fn    <- if (major >= 5) {
        SeuratObject::CreateAssay5Object
      } else {
        SeuratObject::CreateAssayObject
      }
      suppressWarnings(
        sc[[name]] <- fn(data = as.matrix(Matrix::t(enrichment)))
      )
    } else {
      warning("SeuratObject package is required to add enrichment to Seurat object.")
    }
    
  } else if (.is_sce(sc)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      alt <- SummarizedExperiment::SummarizedExperiment(
        assays = list(data = Matrix::t(enrichment))
      )
      SingleCellExperiment::altExp(sc, name) <- alt
    } else {
      warning("SummarizedExperiment and SingleCellExperiment packages are required to add enrichment to SCE object.")
    }
  }
  
  sc
}

.pull.Enrich <- function(sc, name) {
  if (.is_seurat(sc)) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      Matrix::t(sc[[name]]["data"])
    } else {
      stop("Matrix package is required to transpose Seurat assay data.")
    }
    
  } else if (.is_sce(sc)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      Matrix::t(SummarizedExperiment::assay(SingleCellExperiment::altExp(sc)[[name]]))
    } else {
      stop("SummarizedExperiment and SingleCellExperiment packages are required to pull enrichment from SCE object.")
    }
    
  } else {
    stop("Unsupported object type for pulling enrichment.")
  }
}

# -----------------------------------------------------------------------------
#  GENE‑SET / META HELPERS -----------------------------------------------------
# -----------------------------------------------------------------------------
.GS.check <- function(gene.sets) {
  if (is.null(gene.sets))
    stop("Please supply 'gene.sets'")
  if (inherits(gene.sets, "GeneSetCollection"))
    return(GSEABase::geneIds(gene.sets))
  gene.sets
}

.grabMeta <- function(sc) {
  if (.is_seurat(sc)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject package is required to extract metadata from a Seurat object.")
    }
    out <- data.frame(sc[[]], ident = SeuratObject::Idents(sc))
  } else if (.is_sce(sc)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package is required to extract metadata 
           from a SingleCellExperiment object.")
    }
    cd <- SummarizedExperiment::colData(sc)
    out <- data.frame(cd, stringsAsFactors = FALSE)
    # Preserve rownames explicitly
    rownames(out) <- rownames(cd)
    
    # Ensure 'ident' column exists
    if ("ident" %!in% colnames(out)) {
      out$ident <- NA
    }
  } else {
    stop("Unsupported object type; must be Seurat or SingleCellExperiment.")
  }
  return(out)
}


.grabDimRed <- function(sc, dimRed) {
  if (.is_seurat(sc)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject package is required to access dimensional reduction in Seurat objects.")
    }
    
    red <- sc[[dimRed]]
    return(list(
      PCA          = red@cell.embeddings,
      eigen_values = red@misc$eigen_values,
      contribution = red@misc$contribution,
      rotation     = red@misc$rotation
    ))
    
  } else if (.is_sce(sc)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package is required to access dimensional reduction in SCE objects.")
    }
    
    return(list(
      PCA          = SingleCellExperiment::reducedDim(sc, dimRed),
      eigen_values = sc@metadata$eigen_values,
      contribution = sc@metadata$contribution,
      rotation     = sc@metadata$rotation
    ))
  }
}
# -----------------------------------------------------------------------------
#  Underlying Enrichment Calculations
# -----------------------------------------------------------------------------

#─ Ensures a package is present and attaches quietly; 
.load_backend <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " not installed, install or choose a different `method`.",
         call. = FALSE)
  }
}

#─ Build the *Param* object used by GSVA for classic GSVA / ssGSEA -------------
.build_gsva_param <- function(expr, gene_sets, method) {
  .load_backend("GSVA")
  if (method == "GSVA") {
    GSVA::gsvaParam(exprData = expr, geneSets = gene_sets, kcdf = "Poisson")
  } else {                               # ssGSEA
    GSVA::ssgseaParam(exprData = expr, geneSets = gene_sets, normalize = FALSE)
  }
}

#─ Perform enrichment on one cell chunk ---------------------------------------
.compute_enrichment <- function(expr, gene_sets, method, BPPARAM, ...) {
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (is.null(BPPARAM) || !inherits(BPPARAM, "BiocParallelParam")) {
      BPPARAM <- BiocParallel::SerialParam()   # safe default everywhere
    }
  }
  
  switch(toupper(method),
         "GSVA" = {
           param <- .build_gsva_param(expr, gene_sets, "GSVA")
           GSVA::gsva(param = param, BPPARAM = BPPARAM, verbose = FALSE, ...)
         },
         "SSGSEA" = {
           param <- .build_gsva_param(expr, gene_sets, "ssGSEA")
           GSVA::gsva(param = param, BPPARAM = BPPARAM, verbose = FALSE, ...)
         },
         "UCELL" = {
           .load_backend("UCell")
           Matrix::t(UCell::ScoreSignatures_UCell(matrix  = expr,
                                          features = gene_sets,
                                          name     = NULL,
                                          BPPARAM  = BPPARAM,
                                          ...))
         },
         "AUCELL" = {
           .load_backend("AUCell")
           ranks <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
           SummarizedExperiment::assay(
             AUCell::AUCell_calcAUC(geneSets = gene_sets,
                                    rankings  = ranks,
                                    normAUC   = TRUE,
                                    aucMaxRank = ceiling(0.2 * nrow(expr)),
                                    verbose  = FALSE,
                                    ...))
         },
         stop("Unknown method: ", method, call. = FALSE)
  )
}

#─ Split a matrix into equal‑sized column chunks ------------------------------
.split_cols <- function(mat, chunk) {
  if (ncol(mat) <= chunk) return(list(mat))
  idx <- split(seq_len(ncol(mat)), ceiling(seq_len(ncol(mat)) / chunk))
  lapply(idx, function(i) mat[, i, drop = FALSE])
}

.match_summary_fun <- function(fun) {
  if (is.function(fun)) return(fun)
  
  if (!is.character(fun) || length(fun) != 1L)
    stop("'summary.fun' must be a single character or a function")
  
  kw <- tolower(fun)
  fn <- switch(kw,
               mean      = base::mean,
               median    = stats::median,
               max       = base::max,
               sum       = base::sum,
               geometric = function(x) exp(mean(log(x + 1e-6))),
               stop("Unsupported summary keyword: ", fun))
  attr(fn, "keyword") <- kw               # tag for fast matrixStats branch
  fn
}

.computeRunningES <- function(gene.order, hits, weight = NULL) {
  N   <- length(gene.order)
  hit <- gene.order %in% hits
  Nh  <- sum(hit)
  Nm  <- N - Nh
  if (is.null(weight)) weight <- rep(1, Nh)
  
  Phit          <- rep(0, N)
  Phit[hit]     <- weight / sum(weight)
  Pmiss         <- rep(-1 / Nm, N)
  cumsum(Phit + Pmiss)
}


# Modified from GSVA
#' @importFrom MatrixGenerics rowSds
.filterFeatures <- function(expr) {
  sdGenes <- rowSds(expr)
  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
  }
  
  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")
  
  if(is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")
  expr
}

# Parallel-aware lapply
.plapply <- function(X, FUN, ..., BPPARAM = NULL, parallel = TRUE) {
  if (parallel && requireNamespace("BiocParallel", quietly = TRUE)) {
    if (is.null(BPPARAM)) {            
      BPPARAM <- BiocParallel::SerialParam()
    }
    BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
  } else {
    lapply(X, FUN, ...)
  }
}

utils::globalVariables(c(
  "ES", "grp", "x", "y", "xend", "yend", "group", "value", "variable",
  "gene.set.query", "index"
))


