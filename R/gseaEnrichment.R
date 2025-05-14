#' Classical GSEA-style Running-Enrichment Plot
#'
#' Produces the familiar two-panel GSEA graphic—running enrichment score
#' (RES) plus a “hit” rug—for a **single gene-set** evaluated across
#' multiple biological groups (clusters, conditions, samples, …).  
#' The maximal signed deviation of each running-score curve is taken as
#' the enrichment score (**ES**) and printed directly inside the legend
#' label, e.g. `Cluster-A (ES = 1.42)`.  
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
#' @param gene.set.use Character(1).  Name of the gene set to display.
#' @param gene.sets A named list of character vectors, the result of
#' [getGeneSets()], or the built-in data object [escape.gene.sets].
#' @param group.by Metadata column. Defaults to the Seurat/SCE `ident` 
#' slot when `NULL`.
#' @param summary.fun  Method used to collapse expression within each
#* group **before** ranking: one of `"mean"` (default), `"median"`, `"max"`,
#*`"sum"`, or `"geometric"`
#* @param p Weighting exponent in the KS statistic (classical GSEA uses `p = 1`).
#' @param rug.height   Vertical spacing of the hit rug as a fraction of the
#' y-axis (default `0.02`).
#' @param digits       Number of decimal places displayed for ES in the
#' legend (default `2`).
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.

#'
#' @return A single `patchwork`/`ggplot2` object that can be further
#'   modified with `+` (e.g. `+ ggtitle()`).
#'
#' @examples
#' data(pbmc_small)
#'
#' GS <- list(Immune = c("CD3D","CD3E","CD3G","MS4A1","CD79A","CD79B"))

#' gseaEnrichment(pbmc_small,
#'                gene.set.use = "Immune",
#'                gene.sets    = GS,
#'                group.by     = "groups",
#'                summary.fun  = "median",
#'                digits       = 3)
#'
#' @seealso \code{\link{escape.matrix}}, \code{\link{densityEnrichment}}
gseaEnrichment <- function(input.data,
                           gene.set.use,
                           gene.sets,
                           group.by    = NULL,
                           summary.fun = "mean",
                           p           = 1,
                           rug.height  = 0.02,
                           digits      = 2,
                           palette     = "inferno") {
  
  ## ---------- 0  Checks (unchanged) ----------------------------------------
  gene.sets <- .GS.check(gene.sets)
  if (length(gene.set.use) != 1L)
    stop("'gene.set.use' must be length 1")
  if (!gene.set.use %in% names(gene.sets))
    stop("Unknown gene-set")
  
  if (is.null(group.by)) group.by <- "ident"
  meta <- .grabMeta(input.data)
  if (!group.by %in% colnames(meta))
    stop("'", group.by, "' not found in metadata")
  
  groups <- na.omit(unique(meta[[group.by]]))
  if (length(groups) < 2)
    stop("Need 2 groups or more")
  
  summary.fun <- .match_summary_fun(summary.fun)
  
  ## ---------- 1  Expression & ranking vectors ------------------------------
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  cnts <- .filterFeatures(cnts)
  
  gene.order <- rownames(cnts)
  gs.genes   <- intersect(gene.sets[[gene.set.use]], gene.order)
  if (!length(gs.genes))
    stop("Gene-set has no overlap with the matrix")
  
  getStats <- function(mat) {
    switch(attr(summary.fun, "keyword"),
           mean      = MatrixGenerics::rowMeans2(mat),
           median    = matrixGenerics::rowMedians(mat),
           max       = matrixGenerics::rowMaxs(mat),
           sum       = matrixGenerics::rowSums2(mat),
           geometric = exp(matrixGenerics::rowMeans2(log(mat + 1e-6))),
           summary.fun(mat))
  }
  
  ranking.list <- lapply(groups, function(g) {
    idx  <- which(meta[[group.by]] == g)
    lib  <- Matrix::colSums(cnts[, idx, drop = FALSE])
    stat <- getStats(t(t(cnts[, idx, drop = FALSE]) / lib) * 1e6)
    sort(stat, decreasing = TRUE)
  })
  names(ranking.list) <- groups
  
  ## ---------- 2  Running ES & add ES to legend ------------------------------
  es.vec  <- numeric(length(groups))
  curves  <- vector("list", length(groups))
  
  for (i in seq_along(groups)) {
    rvec        <- ranking.list[[i]]
    weight      <- abs(rvec[gs.genes])^p
    curves[[i]] <- .computeRunningES(names(rvec), gs.genes, weight)
    es.vec[i]   <- ifelse(max(abs(curves[[i]])) == abs(max(curves[[i]])),
                          max(curves[[i]]), min(curves[[i]]))
  }
  
  # Build pretty legend labels: Group (ES = 1.23)
  pretty.grp <- paste0(groups,
                       " (ES = ", formatC(es.vec, digits = digits, format = "f"),
                       ")")
  
  ## ---------- 3  Data frames for ggplot -------------------------------------
  running.df <- data.frame(
    rank = rep(seq_along(ranking.list[[1]]), times = length(groups)),
    ES   = unlist(curves, use.names = FALSE),
    grp  = factor(rep(pretty.grp, each = length(curves[[1]])),
                  levels = pretty.grp)
  )
  
  rug.df <- do.call(rbind, lapply(seq_along(groups), function(i) {
    data.frame(x    = which(names(ranking.list[[i]]) %in% gs.genes),
               y    = -(i-1)*rug.height,
               xend = which(names(ranking.list[[i]]) %in% gs.genes),
               yend = -(i)*rug.height,
               grp  = pretty.grp[i])
  }))
  
  ## ---------- 4  Plot -------------------------------------------------------
  cols <- .colorizer(palette, length(groups))
  
  p_top <- ggplot(running.df, aes(rank, ES, colour = grp)) +
    geom_step(linewidth = 0.8) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(y = "Running Enrichment Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  p_mid <- ggplot(rug.df) +
    geom_segment(aes(x, y, xend = xend, yend = yend, colour = grp)) +
    scale_colour_manual(values = cols, guide = "none") +
    theme_void() +
    ylim(-length(groups)*rug.height, 0)
  
  p_top / p_mid + patchwork::plot_layout(heights = c(3, 0.4))
}

#---------------- Helper: wrap summary.fun keyword ---------------------------
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

#------------ Helper: running ES (unchanged) ---------------------------------
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
