#' Visualize Enrichment Distributions Using Geyser Plots
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a geyser plot.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single‑cell
#' object previously processed by \code{\link{runEscape}}.
#' @param assay Name of the assay holding enrichment scores when
#' `input.data` is a single‑cell object. Ignored otherwise.
#' @param group.by Metadata column plotted on the *x*‑axis.  Defaults to the
#' Seurat/SCE `ident` slot when `NULL`.
#' @param gene.set Character(1). Gene‑set to plot (must exist in the
#' enrichment matrix).
#' @param color.by Aesthetic mapped to point color. Use either
#' *"group"* (default = `group.by`) for categorical coloring or the
#' *name of a gene‑set* (e.g. same as `gene.set`) to obtain a numeric
#  gradient. Any other metadata or column present in the data is also
#' accepted.
#' @param order.by How to arrange the x‑axis:
#'   *`"mean"`* – groups ordered by decreasing group mean;
#'   *`"group"`* – natural sort of group labels;
#'   *`NULL`* – keep original ordering.
#' @param facet.by Optional metadata column used to facet the plot.
#' @param scale Logical; if `TRUE` scores are centred/scaled (Z‑score) prior
#' to plotting.
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs,
#'             min.size = NULL)
#'
#' geyserEnrichment(pbmc, 
#'                  assay = "escape", 
#'                  gene.set = "Tcells")
#'
#' @import ggplot2
#' @importFrom ggdist stat_pointinterval
geyserEnrichment <- function(input.data,
                             assay     = NULL,
                             group.by  = NULL,
                             gene.set,
                             color.by  = "group",
                             order.by  = NULL,
                             scale     = FALSE,
                             facet.by  = NULL,
                             palette   = "inferno") {
  ## ---- 0) Sanity checks -----------------------------------------------------
  if (missing(gene.set) || length(gene.set) != 1L)
    stop("Please supply exactly one 'gene.set' to plot.")
  
  if (is.null(group.by))
    group.by <- "ident"
  
  if (identical(color.by, "group"))
    color.by <- group.by
  
  ## ---- 1) Build tidy data.frame -------------------------------------------
  enriched <- .prepData(input.data, assay, gene.set, group.by,
                        split.by = NULL, facet.by = facet.by)
  
  ## Optionally Z‑transform ----------------------------------------------------
  if (scale)
    enriched[[gene.set]] <- as.numeric(scale(enriched[[gene.set]]))
  
  ## Optionally reorder groups -------------------------------------------------
  if (!is.null(order.by))
    enriched <- .orderFunction(enriched, order.by, group.by)
  
  ## ---- 2) Plot --------------------------------------------------------------
  plt <- ggplot(enriched, aes(x = .data[[group.by]],
                              y = .data[[gene.set]],
                              colour = .data[[color.by]])) +
    # Raw points --------------------------------------------------------------
  geom_jitter(width = 0.25, size = 1.5, alpha = 0.6, na.rm = TRUE) +
    
    # White base interval + median point -------------------------------------
  stat_pointinterval(interval_size_range = c(2, 3), fatten_point = 1.4,
                     interval_colour = "white", point_colour = "white",
                     position = position_dodge(width = 0.6), show.legend = FALSE) +
    
    # Black outline for clarity ----------------------------------------------
  stat_pointinterval(interval_size_range = c(1, 2), fatten_point = 1.4,
                     interval_colour = "black", point_colour = "black",
                     position = position_dodge(width = 0.6), show.legend = FALSE) +
    
    labs(x        = group.by,
         y        = paste0(gene.set, "\nEnrichment Score"),
         colour   = color.by) +
    theme_classic() +
    theme(legend.direction = "horizontal",
          legend.position  = "bottom")
  
  ## ---- 3) Colour scale ------------------------------------------------------
  plt <- .colorby(enriched, plt, color.by, palette, type = "color")
  
  ## ---- 4) Facetting ---------------------------------------------------------
  if (!is.null(facet.by))
    plt <- plt + facet_grid(as.formula(paste(".~", facet.by)))
  
  plt
}
