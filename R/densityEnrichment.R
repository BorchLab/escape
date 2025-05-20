#' Visualize Mean Density Ranking of Genes Across Gene Sets
#' 
#' This function allows to the user to examine the mean ranking
#' within the groups across the gene set. The visualization uses
#' the density function to display the relative position and distribution
#' of rank.
#'
#' @param input.data A \link[SeuratObject]{Seurat} object or a
#' \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param gene.set.use Character(1).  Name of the gene set to display.
#' @param gene.sets A named list of character vectors, the result of
#' [getGeneSets()], or the built-in data object [escape.gene.sets].
#' @param group.by Metadata column. Defaults to the Seurat/SCE `ident` 
#' slot when `NULL`.
#' @param rug.height Vertical spacing of the hit rug as a fraction of the
#' y-axis (default `0.02`).
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc_small <- SeuratObject::pbmc_small
#'                         
#' densityEnrichment(pbmc_small, 
#'                   gene.set.use = "Tcells",
#'                   gene.sets = gs)
#'
#' @return A `patchwork`/`ggplot2` object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats na.omit
#' @importFrom MatrixGenerics rowMeans2
densityEnrichment <- function(input.data,
                              gene.set.use,
                              gene.sets,
                              group.by = NULL,
                              rug.height  = 0.02,
                              palette  = "inferno") {
  ## -------- 0  Input checks --------------------------------------------------
  .checkSingleObject(input.data)
  if (is.null(group.by)) group.by <- "ident"
  
  gene.sets <- .GS.check(gene.sets)
  if (!gene.set.use %in% names(gene.sets))
    stop("'gene.set.use' not found in 'gene.sets'")
  
  ## -------- 1  Counts & grouping --------------------------------------------
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts") |>
    .filterFeatures()
  
  meta   <- .grabMeta(input.data)
  groups <- na.omit(unique(meta[[group.by]]))
  
  ## -------- 2  Fast rank computation per group ------------------------------
  n.genes <- nrow(cnts)
  weights <- abs(seq(n.genes, 1) - n.genes/2)       # fixed triangular weight
  rank.mat <- matrix(NA_integer_, n.genes, length(groups),
                     dimnames = list(rownames(cnts),
                                     paste0(group.by, ".", groups)))
  
  compute.cdf <- utils::getFromNamespace("compute.gene.cdf", "GSVA")
  
  for (i in seq_along(groups)) {
    cols <- which(meta[[group.by]] == groups[i])
    tmp  <- cnts[, cols, drop = FALSE]
    
    dens <- suppressWarnings(
      compute.cdf(tmp, seq_len(ncol(tmp)), TRUE, FALSE)
    )
    ord  <- apply(dens, 2, order, decreasing = TRUE)          # genes × cells
    scores <- vapply(seq_len(ncol(ord)),
                     function(j) weights[ord[, j]],
                     numeric(n.genes))
    
    mean.score     <- rowMeans2(scores)
    rank.mat[, i]  <- round(rank(mean.score, ties.method = "average") / 2)
  }
  
  ## -------- 3  Long data.frame w/o extra deps -------------------------------
  in.set <- rownames(rank.mat) %in% gene.sets[[gene.set.use]]
  long.df <- data.frame(
    value          = as.vector(rank.mat),
    variable       = rep(colnames(rank.mat), each = n.genes),
    gene.set.query = rep(ifelse(in.set, "yes", NA_character_), times = length(groups)),
    stringsAsFactors = FALSE
  )
  
  ## -------- 4  Plots ---------------------------------------------------------
  cols <- .colorizer(palette, length(groups))
  plot.df <- subset(long.df, gene.set.query == "yes" & is.finite(value))
  
  p1 <- ggplot(plot.df,
               aes(x = value, fill = variable)) +
    geom_density(alpha = 0.4, colour = "black") +
    scale_fill_manual(values = cols, name = "Group") +
    labs(y = "Rank density") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
  
  ## simple segment plot for mean-rank positions
  offset <- rug.height
  seg.df <- within(plot.df, {
    ord   <- match(variable, unique(variable))
    y     <- -(ord * offset - offset)
    yend  <- y - offset
  })
  
  p2 <- ggplot(seg.df, aes(x = value, xend = value,
                           y = y, yend = yend,
                           colour = variable)) +
    geom_segment(linewidth = 1) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(x = "Mean rank order") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3,1))
}