#' Classical GSEA-style Running-Enrichment Plot
#'
#' Produces the familiar two-panel GSEA graphic—running enrichment score
#' (RES) plus a “hit” rug—for a **single gene-set** evaluated across
#' multiple biological groups (clusters, conditions, samples, ...).  
#'
#' **Algorithm (Subramanian _et al._, PNAS 2005)**  
#' 1. Within every group, library-size-normalise counts to CPM.  
#' 2. Collapse gene expression with `summary.fun` (mean/median/…).  
#' 3. Rank genes (descending) to obtain one ordered list per group.  
#' 4. Compute the weighted Kolmogorov–Smirnov running score  
#'    (weight = \|stat\|^*p*).  
#' 5. ES = maximum signed deviation of the curve.  
#'
#' @param input.data  A \link[SeuratObject]{Seurat} object or a
#' \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param gene.set.use Character(1). Name of the gene set to display.
#' @param gene.sets A named list of character vectors, the result of
#' [getGeneSets()], or the built-in data object [escape.gene.sets].
#' @param group.by Metadata column. Defaults to the Seurat/SCE `ident` 
#' slot when `NULL`.
#' @param summary.fun Method used to collapse expression within each
#' group **before** ranking: one of `"mean"` (default), `"median"`, `"max"`,
#'`"sum"`, or `"geometric"`.
#' @param p Weighting exponent in the KS statistic (classical GSEA uses `p = 1`).
#' @param nperm Integer ≥ 0. Gene-label permutations per group (default 1000). 
#' `0` value will skip NES/*p* calculation.
#' @param rug.height Vertical spacing of the hit rug as a fraction of the
#' y-axis (default `0.02`).
#' @param digits Number of decimal places displayed for ES in the
#' legend (default `2`).
#' @param BPPARAM A \pkg{BiocParallel} parameter object describing the
#' parallel backend. 
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#'  
#' @examples
#' pbmc_small <- SeuratObject::pbmc_small
#'
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' gseaEnrichment(pbmc_small,
#'                gene.set.use = "Bcells",
#'                gene.sets    = gs,
#'                group.by     = "groups",
#'                summary.fun  = "mean",
#'                digits       = 3)
#'
#' @seealso \code{\link{escape.matrix}}, \code{\link{densityEnrichment}}
#' @importFrom stats na.omit
#' @return A single `patchwork`/`ggplot2` object
#' @export
gseaEnrichment <- function(input.data,
                           gene.set.use,
                           gene.sets,
                           group.by    = NULL,
                           summary.fun = "mean",
                           p           = 1,
                           nperm       = 1000,
                           rug.height  = 0.02,
                           digits      = 2,
                           BPPARAM     = NULL,
                           palette     = "inferno") {
  
  ## ---- 0.  Checks ----------------------------------------------------------
  gene.sets <- .GS.check(gene.sets)
  if (length(gene.set.use) != 1L)
    stop("'gene.set.use' must be length 1")
  if (!gene.set.use %in% names(gene.sets))
    stop("Unknown gene-set")
  
  if (is.null(group.by)) group.by <- "ident"
  meta <- .grabMeta(input.data)
  if (!group.by %in% colnames(meta))
    stop("'", group.by, "' not found in metadata")
  
  groups <- stats::na.omit(unique(meta[[group.by]]))
  if (length(groups) < 2)
    stop("Need 2 or more groups")
  
  summary.fun <- .match_summary_fun(summary.fun)
  
  ## ---- 1.  Expression matrix & rankings ------------------------------------
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts") |>
    .filterFeatures()
  
  gs.genes <- intersect(gene.sets[[gene.set.use]], rownames(cnts))
  if (!length(gs.genes))
    stop("Gene-set has no overlap with the matrix")
  
  getStats <- function(mat) {
    keyword <- attr(summary.fun, "keyword")
    switch(keyword,
           mean      = MatrixGenerics::rowMeans2(mat),
           median    = MatrixGenerics::rowMedians(mat),
           max       = MatrixGenerics::rowMaxs(mat),
           sum       = MatrixGenerics::rowSums2(mat),
           geometric = exp(MatrixGenerics::rowMeans2(log1p(mat))))  # log1p is sparse-safe
  }
  
  ranking.list <- lapply(groups, function(g) {
    idx  <- which(meta[[group.by]] == g)
    lib  <- Matrix::colSums(cnts[, idx, drop = FALSE]) / 1e6  # CPM scale
    sub  <- cnts[, idx, drop = FALSE]
    
    # Sparse-safe column normalization using Diagonal
    norm <- sub %*% Matrix::Diagonal(x = 1 / lib)
    
    stat <- getStats(norm)
    sort(stat, decreasing = TRUE)
  })
  names(ranking.list) <- groups
  n.genes <- length(ranking.list[[1L]])
  
  ## ---- 2.  ES, NES, p-value per group --------------------------------------
  es      <- nes <- pval <- numeric(length(groups))
  curves  <- vector("list", length(groups))
  
  for (i in seq_along(groups)) {
    rvec        <- ranking.list[[i]]
    weight      <- abs(rvec[gs.genes])^p
    curves[[i]] <- .computeRunningES(names(rvec), gs.genes, weight)
    es[i]       <- ifelse(max(abs(curves[[i]])) == abs(max(curves[[i]])),
                          max(curves[[i]]), min(curves[[i]]))
    
    ## ---- permutation null --------------------------------------------------
    if (nperm > 0) {
      nullES <- .plapply(
        seq_len(nperm),
        function(xx) {
          hits   <- sample.int(n.genes, length(gs.genes))
          weight <- abs(rvec[hits])^p
          cur    <- .computeRunningES(names(rvec), names(rvec)[hits], weight)
          ifelse(max(abs(cur)) == abs(max(cur)), max(cur), min(cur))
        },
        BPPARAM  = BPPARAM,   # will be ignored in serial mode
        parallel = TRUE       # set FALSE to force serial execution
      )
      nullES <- unlist(nullES, use.names = FALSE)
      
      nes[i]  <- es[i] / mean(abs(nullES))
      pval[i] <- (sum(abs(nullES) >= abs(es[i])) + 1) / (nperm + 1)
    } else {
      nes[i]  <- NA_real_
      pval[i] <- NA_real_
    }
  }
  
  ## ---- 3.  Legend labels ----------------------------------------------------
  labES  <- formatC(es,  digits = digits, format = "f")
  labNES <- formatC(nes, digits = digits, format = "f")
  labP   <- ifelse(is.na(pval), "NA",
                   formatC(pval, digits = 2, format = "e"))
  pretty.grp <- paste0(groups,
                       " (NES = ", labNES,
                       ", p = ", labP, ")")
  
  ## ---- 4.  Data frames for ggplot ------------------------------------------
  running.df <- data.frame(
    rank = rep(seq_len(n.genes), times = length(groups)),
    ES   = unlist(curves, use.names = FALSE),
    grp  = factor(rep(pretty.grp, each = n.genes), levels = pretty.grp)
  )
  
  rug.df <- do.call(rbind, lapply(seq_along(groups), function(i) {
    data.frame(
      x    = which(names(ranking.list[[i]]) %in% gs.genes),
      y    = -(i-1)*rug.height,
      xend = which(names(ranking.list[[i]]) %in% gs.genes),
      yend = -(i)*rug.height,
      grp  = pretty.grp[i])
  }))
  
  ## ---- 5.  Plot -------------------------------------------------------------
  cols <- .colorizer(palette, length(groups))
  
  p_top <- ggplot2::ggplot(running.df, ggplot2::aes(rank, ES, colour = grp)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(y = "Running Enrichment Score") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  
  p_mid <- ggplot2::ggplot(rug.df) +
    ggplot2::geom_segment(ggplot2::aes(x, y, xend = xend, yend = yend,
                                       colour = grp)) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    theme_classic() +
    ggplot2::ylim(-length(groups)*rug.height, 0) + 
    theme(axis.title = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5))
   
  p_top / p_mid + patchwork::plot_layout(heights = c(3, 0.4))
}

