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
#' @param summarise.by Optional metadata column used to summarise data.
#' @param summary.stat Optional method used to summarize expression within each
#'   group defined by \code{summarise.by}. One of: \code{"mean"} (default),
#'   \code{"median"}, \code{"max"}, \code{"sum"}, or \code{"geometric"}.
#' @param scale Logical; if `TRUE` scores are centered/scaled (Z‑score) prior
#' to plotting.
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
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
#' @importFrom stats as.formula
#' @return A \pkg{ggplot2} object.
#' @export
geyserEnrichment <- function(input.data,
                             assay     = NULL,
                             group.by  = NULL,
                             gene.set,
                             color.by  = "group",
                             order.by  = NULL,
                             scale     = FALSE,
                             facet.by  = NULL,
                             summarise.by = NULL,
                             summary.stat   = "mean",
                             palette   = "inferno") {
  ## ---- 0) Sanity checks -----------------------------------------------------
  if (missing(gene.set) || length(gene.set) != 1L)
    stop("Please supply exactly one 'gene.set' to plot.")
  
  if (is.null(group.by))
    group.by <- "ident"
  
  if (identical(color.by, "group"))
    color.by <- group.by
  
  if (!is.null(summarise.by) && (identical(summarise.by, group.by) || 
      identical(summarise.by, facet.by)))
    stop("'summarise.by' cannot be the same as 'group.by' or 'facet.by'. 
         Please choose a different metadata column.")
  
  # ---- 1) helper to match summary function -------------------------
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
    fn
  }
  summary_fun <- .match_summary_fun(summary.stat)
  
  ## ---- 2) Build tidy data.frame -------------------------------------------
  enriched <- .prepData(input.data, assay, gene.set, group.by,
                        split.by = summarise.by, facet.by = facet.by, color.by = color.by)
  
  ## Optionally summarise data with **base aggregate()** ----------------------
  if (!is.null(summarise.by)) {
    grp_cols   <- c(summarise.by, group.by, facet.by, color.by)
    enriched <- aggregate(enriched[gene.set],
                          by   = enriched[grp_cols],
                          FUN  = summary_fun,
                          SIMPLIFY = FALSE)
  }

  ## Optionally Z‑transform ----------------------------------------------------
  if (scale)
    enriched[[gene.set]] <- as.numeric(scale(enriched[[gene.set]]))
  
  ## Optionally reorder groups -------------------------------------------------
  if (!is.null(order.by))
    enriched <- .orderFunction(enriched, order.by, group.by)
  
  ## ---- 3) Plot --------------------------------------------------------------
  if (!is.null(color.by))
    plt <- ggplot(enriched, aes(x = .data[[group.by]],
                                y = .data[[gene.set]],
                                group = .data[[group.by]],
                                colour = .data[[color.by]]))
  else
    plt <- ggplot(enriched, aes(x = .data[[group.by]],
                                y = .data[[gene.set]]),
                                group = .data[[group.by]])

    # Raw points --------------------------------------------------------------
  plt <- plt + geom_jitter(width = 0.25, size = 1.5, alpha = 0.6, na.rm = TRUE) +
    
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
  
  ## ---- 4) Colour scale ------------------------------------------------------
  if (!is.null(color.by)) 
    plt <- .colorby(enriched, plt, color.by, palette, type = "color")
  
  ## ---- 5) Facetting ---------------------------------------------------------
  if (!is.null(facet.by))
    plt <- plt + facet_grid(as.formula(paste(".~", facet.by)))
  
  plt
}
