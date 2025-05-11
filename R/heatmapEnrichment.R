#' Generate a heatmap to visualize enrichment values
#' 
#' This function allows to the user to examine the heatmap with the mean
#' enrichment values by group. The heatmap will have the gene sets as rows
#' and columns will be the grouping variable.
#'
#' @param input.data  Output of \code{\link{escape.matrix}} or a single‑cell
#'   object previously processed by \code{\link{runEscape}}.
#' @param assay       Name of the assay containing enrichment data when
#'   `input.data` is a single‑cell object.
#' @param group.by    Metadata column used to define columns in the heatmap.
#'   Defaults to the Seurat/SCE `ident` slot.
#' @param gene.set.use Vector of gene‑set names to plot, or \code{"all"}
#'   (default) to show every available gene set.
#' @param cluster.rows,cluster.columns Logical; if \code{TRUE}, rows/columns
#'   are ordered by Ward‑linkage hierarchical clustering (Euclidean distance).
#' @param facet.by    Optional metadata column to facet the heatmap.
#' @param scale       If \code{TRUE}, Z‑transforms each gene‑set column _after_
#'   summarisation.
#' @param summary.stat Character keyword (\code{"mean"}, \code{"median"},
#'   \code{"sum"}, \code{"sd"}, \code{"max"}, \code{"min"},
#'   \code{"geometric"}) **or** a custom function to collapse scores within
#'   each group.  Defaults to \code{"mean"}.
#' @param palette     Any palette from \link[grDevices]{hcl.pals}; default
#'   \code{"inferno"}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' gs <- list(B  = c("MS4A1","CD79B","CD79A"),
#'            T  = c("CD3D","CD3E","CD3G"))
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#' heatmapEnrichment(pbmc, assay = "escape", palette = "viridis")
heatmapEnrichment <- function(input.data,
                              assay          = NULL,
                              group.by       = NULL,
                              gene.set.use   = "all",
                              cluster.rows   = FALSE,
                              cluster.columns= FALSE,
                              facet.by       = NULL,
                              scale          = FALSE,
                              summary.stat   = "mean",
                              palette        = "inferno") {
  #---------------------- helper: match/validate summary function -------------
  .match_summary_fun <- function(fun) {
    if (is.function(fun)) return(fun)
    if (!is.character(fun) || length(fun) != 1L)
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
    attr(fn, "keyword") <- kw
    fn
  }
  summary_fun <- .match_summary_fun(summary.stat)
  
  #---------------------- defaults & data extraction --------------------------
  if (is.null(group.by)) group.by <- "ident"
  df <- .prepData(input.data, assay, gene.set.use, group.by, NULL, facet.by)
  
  # determine gene‑set columns ------------------------------------------------
  if (identical(gene.set.use, "all")) {
    gene.set <- setdiff(colnames(df), c(group.by, facet.by))
  } else {
    gene.set <- gene.set.use
  }
  if (!length(gene.set))
    stop("No gene‑set columns found to plot.")
  
  #---------------------- summarise ------------------------------------------
  if (is.null(facet.by)) {
    grp <- df[[group.by]]
    agg <- aggregate(df[gene.set], by = list(!!group.by := grp), FUN = summary_fun)
  } else {
    grp  <- df[[group.by]]
    fac  <- df[[facet.by]]
    agg  <- aggregate(df[gene.set],
                      by = list(!!group.by := grp, !!facet.by := fac),
                      FUN = summary_fun)
  }
  
  # option: Z‑transform per gene‑set column ----------------------------------
  if (scale) {
    agg[gene.set] <- lapply(agg[gene.set], function(col) as.numeric(scale(col)))
  }
  
  #---------------------- reshape for ggplot (base R) -------------------------
  long <- data.frame(
    variable = rep(gene.set, each = nrow(agg)),
    value    = as.vector(t(agg[gene.set])),
    group    = rep(agg[[group.by]], times = length(gene.set)),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet.by)) long[[facet.by]] <- rep(agg[[facet.by]], times = length(gene.set))
  
  #---------------------- optional clustering --------------------------------
  if (cluster.rows) {
    ord <- hclust(dist(t(agg[gene.set])), method = "ward.D2")$order
    long$variable <- factor(long$variable, levels = gene.set[ord])
  }
  if (cluster.columns) {
    ord <- hclust(dist(agg[gene.set]), method = "ward.D2")$order
    col_levels <- agg[[group.by]][ord]
    long$group <- factor(long$group, levels = col_levels)
  }
  
  #---------------------- build ggplot ---------------------------------------
  p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = variable, fill = value)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(colours = .colorizer(palette, 11),
                                  name   = "Enrichment") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title      = ggplot2::element_blank(),
                   axis.ticks      = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.direction= "horizontal")
  if (!is.null(facet.by)) {
    p <- p + ggplot2::facet_grid(stats::as.formula(paste(". ~", facet.by)))
  }
  p
}
