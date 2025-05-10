#Developing split violin plot
#Code from: https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), 
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = 
                                                                  if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], 
                                              newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- 
                               round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), 
                                         all(draw_quantiles <= 1))
                               quantiles <- 
                                 ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), 
                                                  setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", 
                                                grid::grobTree(GeomPolygon$draw_panel(newdata, ...), 
                                                               quantile_grob))
                             } else {
                               ggplot2:::ggname("geom_split_violin", 
                                                GeomPolygon$draw_panel(newdata, ...))}
                           })

#Defining new geometry
#Code from: https://stackoverflow.com/a/45614547
geom_split_violin <- 
  function(mapping = NULL, data = NULL, 
           stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, 
           trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, 
           inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, 
          inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, 
                                                   draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }

#' Plot Enrichment Distributions Using Split or Dodged Violin Plots
#'
#' Visualize the distribution of gene set enrichment scores across groups using
#' violin plots. When `split.by` contains exactly two levels, the function draws
#' split violins for easy group comparison within each `group.by` category. If
#' `split.by` has more than two levels, standard dodged violins are drawn instead.
#'
#' @param input.data A matrix or single-cell object (e.g., Seurat or 
#'   SingleCellExperiment) containing enrichment scores from 
#'   \code{\link{escape.matrix}} or \code{\link{runEscape}}.
#' @param assay Name of the assay containing enrichment scores if `input.data` 
#'   is a single-cell object.
#' @param split.by A metadata column used to split or color violins. Must contain 
#'   at least two levels. If it contains more than two, dodged violins are used.
#' @param group.by Metadata column used for the x-axis grouping. If not specified, 
#'   defaults to \code{"ident"}.
#' @param gene.set Name of the gene set to visualize on the y-axis.
#' @param order.by Method to order the x-axis: either \code{"mean"} to order by 
#'   mean enrichment, \code{"group"} for alphanumerical order, or \code{NULL} 
#'   to retain the original order.
#' @param facet.by Optional metadata column used to facet the plot into multiple panels.
#' @param scale Logical; if \code{TRUE}, enrichment values are Z-transformed 
#'   prior to plotting.
#' @param palette Color palette to use for fill aesthetics. Must be a valid 
#'   palette from \code{\link[grDevices]{hcl.pals}}.
#'
#' @return A \code{ggplot2} object displaying enrichment score distributions by group.
#'
#' @import ggplot2
#' @importFrom grDevices hcl.pals
#'
#' @examples
#' gene.sets <- list(
#'   Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'   Tcells = c("CD3E", "CD3D", "CD3G", "CD7", "CD8A")
#' )
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, gene.sets = gene.sets)
#'
#' splitEnrichment(
#'   input.data = pbmc_small,
#'   assay = "escape",
#'   split.by = "groups",
#'   gene.set = "Tcells"
#' )
#'
#' @export
splitEnrichment <- function(input.data,
                            assay = NULL,
                            split.by = NULL,
                            group.by = NULL,
                            gene.set = NULL,
                            order.by = NULL,
                            facet.by = NULL,
                            scale = TRUE,
                            palette = "inferno") {
  
  if (is.null(split.by)) stop("Please specify a variable for 'split.by'.")
  if (is.null(group.by)) group.by <- "ident"
  
  enriched <- .prepData(input.data, assay, gene.set, group.by, split.by, facet.by)
  
  split.levels <- unique(enriched[[split.by]])
  n.levels <- length(split.levels)
  
  if (n.levels < 2) stop("split.by must have at least two levels.")
  
  if (scale) {
    enriched[[gene.set]] <- scale(enriched[[gene.set]])
  }
  
  if (!is.null(order.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  }
  
  plot <- ggplot(enriched, aes(x = .data[[group.by]],
                               y = .data[[gene.set]],
                               fill = .data[[split.by]])) +
    xlab(group.by) +
    ylab(paste0(gene.set, "\n Enrichment Score")) +
    labs(fill = split.by) +
    scale_fill_manual(values = .colorizer(palette, n.levels)) +
    theme_classic()
  
  # Use split violin for binary factors; dodge otherwise
  if (n.levels == 2) {
    plot <- plot +
      geom_split_violin(alpha = 0.8, lwd = 0.25)
  } else {
    plot <- plot +
      geom_violin(position = position_dodge(width = 0.8), alpha = 0.8, lwd = 0.25)
  }
  
  # Add a central boxplot
  plot <- plot +
    geom_boxplot(width = 0.1,
                 fill = "grey",
                 alpha = 0.5,
                 outlier.shape = NA,
                 position = if (n.levels == 2) position_identity() else position_dodge(width = 0.8),
                 notch = TRUE)
  
  # Add faceting if specified
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(as.formula(paste(". ~", facet.by)))
  }
  
  return(plot)
}
