#' Plot Enrichment Distributions Using Split or Dodged Violin Plots
#'
#' Visualize the distribution of gene set enrichment scores across groups using
#' violin plots. When `split.by` contains exactly two levels, the function draws
#' split violins for easy group comparison within each `group.by` category. If
#' `split.by` has more than two levels, standard dodged violins are drawn instead.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single‑cell
#' object previously processed by \code{\link{runEscape}}.
#' @param assay Name of the assay holding enrichment scores when
#' `input.data` is a single‑cell object. Ignored otherwise.
#' @param split.by A metadata column used to split or color violins. Must contain 
#'   at least two levels. If it contains more than two, dodged violins are used.
#' @param group.by Metadata column plotted on the *x*‑axis.  Defaults to the
#' Seurat/SCE `ident` slot when `NULL`.
#' @param gene.set.use Character(1).  Name of the gene set to display.
#' @param order.by How to arrange the x‑axis:
#'   *`"mean"`* – groups ordered by decreasing group mean;
#'   *`"group"`* – natural sort of group labels;
#'   *`NULL`* – keep original ordering.
#' @param facet.by Optional metadata column used to facet the plot.
#' @param scale Logical; if `TRUE` scores are centred/scaled (Z‑score) prior
#' to plotting.
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#'
#' @return A [ggplot2] object.
#'
#' @import ggplot2
#' @importFrom grDevices hcl.pals
#' @importFrom stats as.formula
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'
#' splitEnrichment(input.data = pbmc,
#'                 assay = "escape",
#'                 split.by = "groups",
#'                 gene.set.use = "Tcells")
#'
#' @export
splitEnrichment <- function(input.data,
                            assay = NULL,
                            split.by = NULL,
                            group.by = NULL,
                            gene.set.use = NULL,
                            order.by = NULL,
                            facet.by = NULL,
                            scale = TRUE,
                            palette = "inferno") {
  
  if (is.null(split.by)) stop("Please specify a variable for 'split.by'.")
  if (is.null(group.by)) group.by <- "ident"
  
  # Prepare tidy data with relevant metadata columns
  enriched <- .prepData(input.data, assay, gene.set.use, group.by, split.by, facet.by)
  
  # Determine the number of levels in the splitting variable
  split.levels <- unique(enriched[[split.by]])
  n.levels     <- length(split.levels)
  if (n.levels < 2) stop("split.by must have at least two levels.")
  
  # Optional Z-score scaling of enrichment values
  if (scale) {
    enriched[[gene.set.use]] <- scale(enriched[[gene.set.use]])
  }
  
  # Optional reordering of x-axis categories
  if (!is.null(order.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  }
  
  # Create a composite group for proper boxplot dodging
  enriched$group_split <- interaction(enriched[[group.by]], enriched[[split.by]], sep = "_")
  dodge <- position_dodge(width = 0.8)
  
  # Base plot
  plot <- ggplot(enriched, aes(x = .data[[group.by]],
                               y = .data[[gene.set.use]],
                               fill = .data[[split.by]])) +
    xlab(group.by) +
    ylab(paste0(gene.set.use, "\nEnrichment Score")) +
    labs(fill = split.by) +
    scale_fill_manual(values = .colorizer(palette, n.levels)) +
    theme_classic()
  
  # Split violin if binary, otherwise dodge standard violins
  if (n.levels == 2) {
    plot <- plot +
      geom_split_violin(alpha = 0.8, lwd = 0.25)  +
      geom_boxplot(width = 0.1,
                   fill = "grey",
                   alpha = 0.6,
                   outlier.shape = NA,
                   position = position_identity(),
                   notch = FALSE)
  } else {
    plot <- plot +
      geom_violin(position = dodge, alpha = 0.8, lwd = 0.25) +
      geom_boxplot(width = 0.1,
                   fill = "grey",
                   alpha = 0.6,
                   outlier.shape = NA,
                   position = dodge,
                   notch = FALSE,
                   aes(group = .data$group_split))
  }
  

  
  # Optional faceting
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(as.formula(paste(". ~", facet.by)))
  }
  
  return(plot)
}

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
