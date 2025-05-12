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
    df <- cbind(t(cnts[gene.set, , drop = FALSE]), meta)
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
    # use generic accessor to avoid tight coupling to Seurat internals
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      cnts <- SeuratObject::GetAssayData(obj, assay = assay, slot = type)
    } else {
      cnts <- obj@assays[[assay]][type]
    }
  } else if (.is_sce(obj)) {
    pos <- if (assay == "RNA") "counts" else assay
    cnts <- if (assay == "RNA") SummarizedExperiment::assay(obj, pos)
    else SummarizedExperiment::assay(SingleCellExperiment::altExp(obj), pos)
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
      } else  {
        SeuratObject::CreateAssayObject
      }
      suppressWarnings(sc[[name]] <- fn(data = as.matrix(t(enrichment))))
    }
  } else if (.is_sce(sc)) {
    altExp(sc, name) <- SummarizedExperiment::SummarizedExperiment(assays = list(data = t(enrichment)))
  }
  sc
}

.pull.Enrich <- function(sc, name) {
  if (.is_seurat(sc)) {
    Matrix::t(sc[[name]]["data"])
  } else if (.is_sce(sc)) {
    t(SummarizedExperiment::assay(SingleCellExperiment::altExp(sc)[[name]]))
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
    out <- data.frame(sc[[]], ident = SeuratObject::Idents(sc))
  } else if (.is_sce(sc)) {
    out <- data.frame(SummarizedExperiment::colData(sc))
    rownames(out) <- SummarizedExperiment::colData(sc)@rownames
    if ("ident" %!in% colnames(out))
      out$ident <- NA
  } else {
    stop("Unsupported object type")
  }
  out
}

.grabDimRed <- function(sc, dimRed) {
  if (.is_seurat(sc)) {
    list(PCA = sc[[dimRed]]@cell.embeddings, sc[[dimRed]]@misc)
  } else if (.is_sce(sc)) {
    list(PCA = SingleCellExperiment::reducedDim(sc, dimRed),
         sc@metadata[c("eigen_values", "contribution", "rotation")])
  }
}

# -----------------------------------------------------------------------------
#  Underlying Enrichment Calculations
# -----------------------------------------------------------------------------

#─ Ensures a package is present and attaches quietly; 
.load_backend <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " not installed – install or choose a different `method`.",
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
           t(UCell::ScoreSignatures_UCell(matrix  = expr,
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


