# -----------------------------------------------------------------------------
#  CLASS HELPERS 
# -----------------------------------------------------------------------------
.is_seurat      <- function(x) inherits(x, "Seurat")
.is_sce         <- function(x) inherits(x, "SummarizedExperiment")
.is_seurat_or_sce <- function(x) .is_seurat(x) || .is_sce(x)

.checkSingleObject <- function(obj) {
  if (!.is_seurat_or_sce(obj))
    stop("Expecting a Seurat or SummarizedExperiment object")
}

# `base::%||%` only exists in R >= 4.4.0 and escape declares R (>= 4.1).
`%||%` <- function(x, y) if (is.null(x)) y else x

# Single error message for every "you asked for an assay that isn't there" case,
# so a missing name never reaches assay() as NULL.
.stop_missing_assay <- function(requested, available,
                                what = "assay", where = "the object") {
  stop("Could not find ", what, " '", requested, "' in ", where, ". ",
       if (length(available))
         paste0("Available: ", paste(available, collapse = ", "), ".")
       else "No assays are present.",
       call. = FALSE)
}

# -----------------------------------------------------------------------------
#  ORDERING UTILITY 
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

.alphanumericalSort <- function(x, ignore.case = TRUE) {
  if (!is.character(x)) {
    message("Input 'x' is not a character vector. Attempting to convert.")
    x <- as.character(x)
  }
  unique_elements <- unique(x)
  alpha_part <- gsub("[0-9]+", "", unique_elements, perl = TRUE)
  if (ignore.case) {
    alpha_part <- tolower(alpha_part)
  }
  numeric_part <- suppressWarnings(as.numeric(gsub("[^0-9]", "", unique_elements)))
  sorted_levels <- unique_elements[
    order(alpha_part, numeric_part)
  ]
  
  return(sorted_levels)
}

# -----------------------------------------------------------------------------
#  DATA.frame BUILDERS 
# -----------------------------------------------------------------------------
.makeDFfromSCO <- function(input.data, assay = "escape", gene.set = NULL,
                           group.by = NULL, split.by = NULL, facet.by = NULL, color.by = NULL) {
  if (is.null(assay))
    stop("Please provide assay name")
  
  # Pull count matrix (features) and metadata
  cnts <- .cntEval(input.data, assay = assay, type = "data")
  features <- rownames(cnts)
  meta <- .grabMeta(input.data)
  meta.cols <- colnames(meta)
  
  # All potential column-like arguments
  cols <- unique(c(group.by, split.by, facet.by, color.by))
  
  # Check that each is either metadata or a feature
  bad.cols <- cols[!(cols %in% meta.cols | cols %in% features)]
  if (length(bad.cols) > 0) {
    stop("The following variables are not found in either metadata or features: ", paste(bad.cols, collapse = ", "))
  }
  
  # Determine if color.by is a feature or meta
  is_feature_color <- !is.null(color.by) && color.by %in% features
  is_meta_color <- !is.null(color.by) && color.by %in% meta.cols
  
  # Prepare metadata subset
  meta <- meta[, intersect(cols, meta.cols), drop = FALSE]
  
  # Convert gene.set if "all"
  if (length(gene.set) == 1 && gene.set == "all") {
    gene.set <- features
  }
  
  # Build data frame with expression values
  if (length(gene.set) == 1) {
    df <- cbind(value = cnts[gene.set, ], meta)
    colnames(df)[1] <- gene.set
  } else {
    df <- cbind(Matrix::t(cnts[gene.set, , drop = FALSE]), meta)
  }
  
  # Add color.by feature expression if it's a gene but not in gene.set
  if (is_feature_color && !(color.by %in% gene.set)) {
    df[[color.by]] <- cnts[color.by, ]
  }
  
  return(df)
}


.prepData <- function(input.data, assay, gene.set, group.by, split.by, facet.by, color.by) {
  if (.is_seurat_or_sce(input.data)) {
    df <- .makeDFfromSCO(input.data, assay, gene.set, group.by, split.by, facet.by, color.by)
    
    if (identical(gene.set, "all")) {
      meta_cols <- c(group.by, split.by, facet.by)
      # Do not remove color.by if it's also a feature
      non_gene_color <- if (!is.null(color.by) && color.by %in% colnames(df) && !(color.by %in% gene.set)) color.by else NULL
      gene.set <- setdiff(colnames(df), c(meta_cols, non_gene_color))
    }
    
  } else {
    all.cols <- unique(c(gene.set, group.by, split.by, facet.by, color.by))
    missing.cols <- setdiff(all.cols, colnames(input.data))
    if (length(missing.cols) > 0) {
      stop("The following columns are missing in the input data: ", paste(missing.cols, collapse = ", "))
    }
    
    if (identical(gene.set, "all")) {
      gene.set <- setdiff(colnames(input.data), c(group.by, split.by, facet.by, color.by))
    }
    
    df <- input.data[, unique(c(gene.set, group.by, split.by, facet.by, color.by)), drop = FALSE]
  }
  
  return(df)
}

# -----------------------------------------------------------------------------
#  COLOUR SCALES 
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
#  MATRIX / VECTOR SPLITTERS
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
#  INPUT ASSAY RESOLUTION
# -----------------------------------------------------------------------------
# Translate the user-facing `input.assay` into the layer (Seurat) or main assay
# (SummarizedExperiment) name that .cntEval() should read.
#
#   "auto"       counts for the native backend, logcounts for plaid
#   "counts"     Seurat layer "counts"  / SCE assay "counts"
#   "logcounts"  Seurat layer "data"    / SCE assay "logcounts"
#   <other>      taken literally
.resolve_input_assay <- function(obj, input.assay = "auto",
                                 backend = "native") {
  if (!is.character(input.assay) || length(input.assay) != 1L)
    stop("`input.assay` must be a single string.", call. = FALSE)

  ia <- if (identical(input.assay, "auto")) {
    if (identical(backend, "plaid")) "logcounts" else "counts"
  } else {
    input.assay
  }

  type <- if (.is_seurat(obj)) {
    switch(ia, counts = "counts", logcounts = "data", data = "data", ia)
  } else {
    ia
  }

  .assert_input_assay(obj, type, ia)
  type
}

.assert_input_assay <- function(obj, type, requested) {
  if (.is_seurat(obj)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) return(invisible(TRUE))
    lay <- tryCatch(SeuratObject::Layers(obj, assay = "RNA"),
                    error = function(e) character())
    # Seurat v5 split objects expose "data.1", "data.2", ... rather than "data"
    ok <- length(lay) == 0L ||
      any(lay == type | startsWith(lay, paste0(type, ".")))
    if (!ok)
      stop("input.assay = \"", requested, "\" needs layer '", type,
           "' in the 'RNA' assay, which has: ", paste(lay, collapse = ", "),
           ".\n  Run Seurat::NormalizeData() (and Seurat::JoinLayers() if the ",
           "object is split), or set input.assay = \"counts\".", call. = FALSE)

  } else if (.is_sce(obj)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      return(invisible(TRUE))
    an <- SummarizedExperiment::assayNames(obj)
    if (!type %in% an)
      stop("input.assay = \"", requested, "\" needs assay '", type,
           "', which is not present. Available: ", paste(an, collapse = ", "),
           ".\n  Run scuttle::logNormCounts(), or set input.assay = \"counts\".",
           call. = FALSE)

  } else if (identical(type, "logcounts")) {
    message("`input.data` is a bare matrix; escape assumes it is already ",
            "log-normalized (input.assay = \"logcounts\").")
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
#  EXPRESSION MATRIX EXTRACTOR
# -----------------------------------------------------------------------------
#' @importFrom MatrixGenerics rowSums2
.cntEval <- function(obj, assay = "RNA", type = "counts") {
  if (.is_seurat(obj)) {
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      so_version <- utils::packageVersion("SeuratObject")
      if (so_version >= "5.0.0") {
        # SeuratObject >= 5: always use 'layer' (the 'slot' argument is defunct).
        # For v3 Assay objects inside a v5 environment, the layer API still
        # works — it maps legacy slot names transparently.
        cnts <- SeuratObject::GetAssayData(obj, assay = assay, layer = type)
      } else {
        cnts <- SeuratObject::GetAssayData(obj, assay = assay, slot = type)
      }
    } else {
      cnts <- obj@assays[[assay]][[type]]
    }

  } else if (.is_sce(obj)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      cnts <- if (assay == "RNA") {
        ## `assay = "RNA"` means "the main expression assay"; `type` names it.
        ## Every in-package caller passes type = "counts", so this is inert
        ## unless a caller deliberately asks for logcounts.
        avail <- SummarizedExperiment::assayNames(obj)
        if (!type %in% avail)
          .stop_missing_assay(type, avail, "assay",
                              "the SummarizedExperiment")
        SummarizedExperiment::assay(obj, type)
      } else {
        avail <- SingleCellExperiment::altExpNames(obj)
        if (!assay %in% avail)
          .stop_missing_assay(assay, avail, "altExp",
                              "the SingleCellExperiment")
        SummarizedExperiment::assay(SingleCellExperiment::altExp(obj, assay))
      }
    } else {
      stop("SummarizedExperiment and SingleCellExperiment packages are required but not installed.")
    }
  } else {
    cnts <- obj
  }
  # Conditionally require DelayedMatrixStats if cnts is dgCMatrix
  if (inherits(cnts, "dgCMatrix")) {
    if (!requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
      stop("Package 'DelayedMatrixStats' is required to handle sparse matrices. Please install it with BiocManager::install('DelayedMatrixStats').")
    }
    loadNamespace("DelayedMatrixStats")
  }
  cnts[MatrixGenerics::rowSums2(cnts, na.rm = TRUE) != 0, , drop = FALSE]
}

# -----------------------------------------------------------------------------
#  ATTACH / PULL ENRICHMENT MATRICES 
# -----------------------------------------------------------------------------
.adding.Enrich <- function(sc, enrichment, name) {
  if (.is_seurat(sc)) {
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      so_version <- utils::packageVersion("SeuratObject")
      fn <- if (so_version >= "5.0.0") {
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
  if (is.null(name) || !is.character(name) || length(name) != 1L)
    stop("`assay` must be a single enrichment assay name.", call. = FALSE)

  if (.is_seurat(sc)) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      avail <- SeuratObject::Assays(sc)
      if (!name %in% avail)
        .stop_missing_assay(name, avail, "enrichment assay",
                            "the Seurat object")
      so_version <- utils::packageVersion("SeuratObject")
      if (so_version >= "5.0.0") {
        Matrix::t(SeuratObject::GetAssayData(sc, assay = name, layer = "data"))
      } else {
        Matrix::t(SeuratObject::GetAssayData(sc, assay = name, slot = "data"))
      }
    } else {
      stop("Matrix package is required to transpose Seurat assay data.")
    }

  } else if (.is_sce(sc)) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      avail <- SingleCellExperiment::altExpNames(sc)
      if (!name %in% avail)
        .stop_missing_assay(name, avail, "enrichment altExp",
                            "the SingleCellExperiment")
      # altExp(sc, name), not altExp(sc)[[name]] - the latter indexes colData of
      # the *first* altExp and silently returns NULL.
      Matrix::t(SummarizedExperiment::assay(
        SingleCellExperiment::altExp(sc, name)))
    } else {
      stop("SummarizedExperiment and SingleCellExperiment packages are required to pull enrichment from SCE object.")
    }

  } else {
    stop("Unsupported object type for pulling enrichment.")
  }
}

# -----------------------------------------------------------------------------
#  GENE-SET / META HELPERS 
# -----------------------------------------------------------------------------
.GS.check <- function(gene.sets) {
  if (is.null(gene.sets))
    stop("Please supply 'gene.sets'")
  if (inherits(gene.sets, "GeneSetCollection"))
    return(GSEABase::geneIds(gene.sets))
  gene.sets
}

# Align a gene-set list to the columns of an enrichment matrix.
#
# Seurat coerces underscores in feature names to hyphens when an assay is built,
# so scores pulled back off a Seurat object have mangled column names while
# scores held as a plain matrix (or in an SCE altExp) keep the originals. Match
# the literal names first and fall back to the mangled form only for the columns
# that are still unmatched - mangling unconditionally silently drops every
# HALLMARK_/GO_/REACTOME_ set for non-Seurat input.
#
# Returns the sets re-keyed and re-ordered to `cols`, one per column.
.match_sets_to_cols <- function(gene.sets, cols) {
  nm <- names(gene.sets)
  if (is.null(nm))
    stop("`gene.sets` must be a named list.", call. = FALSE)
  if (is.null(cols))
    stop("Enrichment matrix has no column names; cannot match gene sets.",
         call. = FALSE)

  idx  <- match(cols, nm)
  miss <- is.na(idx)
  if (any(miss))
    idx[miss] <- match(cols[miss], gsub("_", "-", nm, fixed = TRUE))

  if (all(is.na(idx)))
    stop("None of the supplied gene sets match enrichment columns.",
         call. = FALSE)
  if (anyNA(idx))
    stop("No gene set supplied for enrichment column(s): ",
         paste(cols[is.na(idx)], collapse = ", "),
         ". Supply the same `gene.sets` used to compute the scores.",
         call. = FALSE)

  out <- gene.sets[idx]
  names(out) <- cols
  out
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
    if (!"ident" %in% colnames(out)) {
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

# -----------------------------------------------------------------------------
#  THEME HELPER
# -----------------------------------------------------------------------------
#' @importFrom ggplot2 %+replace%
.themeEscape <- function(base_size = 12,
                         base_family = "sans",
                         grid_lines = "Y",
                         axis_lines = FALSE,
                         legend_position = "right") {

  t <- ggplot2::theme_bw(base_size = base_size, base_family = base_family)

  t <- t %+replace%
    ggplot2::theme(
      # Plot titles and caption
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1.2), hjust = 0, face = "bold",
        margin = ggplot2::margin(b = base_size / 2)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, face = "italic",
        margin = ggplot2::margin(b = base_size)
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.8), hjust = 1
      ),

      # Axis titles and text
      axis.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold"),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.85)),

      # Facet strips
      strip.text = ggplot2::element_text(
        size = ggplot2::rel(0.9),
        face = "bold",
        margin = ggplot2::margin(base_size / 2.5, base_size / 2.5, base_size / 2.5, base_size / 2.5)
      ),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),

      # Panel border and background
      panel.border = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),

      # Legend styling
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = legend_position,
      legend.key = ggplot2::element_rect(fill = "white"),

      # Plot spacing/margins
      plot.margin = ggplot2::margin(base_size, base_size, base_size, base_size)
    )

  # Handle grid lines
  grid_lines <- toupper(grid_lines)
  if (grid_lines == "NONE") {
    t <- t + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else if (grid_lines == "X") {
    t <- t + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey85", linewidth = 0.25)
    )
  } else if (grid_lines == "Y") {
    t <- t + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.25)
    )
  } else {
    t <- t + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank()
    )
  }

  # Handle axis lines
  if (isTRUE(axis_lines)) {
    t <- t + ggplot2::theme(
      axis.line = ggplot2::element_line(color = "grey30", linewidth = 0.5)
    )
  }

  return(t)
}

# helper to match summary function
.match_summary_fun <- function(fun) {
  if (is.function(fun)) return(fun)
  if (!is.character(fun) || length(fun) != 1)
    stop("'summary.stat' must be a single character keyword or a function")
  kw <- tolower(fun)
  fn <- switch(kw,
               mean      = base::mean,
               median    = stats::median,
               sum       = base::sum,
               sd        = stats::sd,
               max       = base::max,
               min       = base::min,
               geometric = function(x) exp(mean(log(x + 1e-6))),
               stop("Unsupported summary keyword: ", fun))
  
  # Attach keyword as attribute
  attr(fn, "keyword") <- kw
  fn
}

