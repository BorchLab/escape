#' Plot 2D Enrichment Distributions With Density or Hexplots
#'
#' Visualize the relationship between *two* enrichment scores at single-cell
#' resolution. By default points are shaded by local 2-D density
#' (`color.by = "density"`), but users can instead color by a metadata column
#' (discrete) or by the raw gene-set scores themselves (continuous).
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single‑cell
#' object previously processed by \code{\link{runEscape}}.
#' @param assay Name of the assay holding enrichment scores when
#' `input.data` is a single‑cell object. Ignored otherwise.
#' @param x.axis,y.axis Gene-set names to plot on the *x* and *y* axes.
#' @param facet.by Optional metadata column used to facet the plot.
#' @param group.by Metadata column plotted.  Defaults to the
#' Seurat/SCE `ident` slot when `NULL`.
#' @param color.by Aesthetic mapped to point color. Use 
#' `"density"` (default), `"group"`, `"x"`, or `"y"`.  The latter two apply a 
#' continuous gradient to the corresponding axis.
#' @param style `"point"` (density-aware points) or `"hex"` (hex-bin).
#' @param scale Logical; if `TRUE` scores are centered/scaled (Z‑score) prior
#' to plotting.
#' @param bins Number of hex bins along each axis when `style = "hex"`.
#' @param point.size,alpha  Aesthetic tweaks for `style = "point"`.
#' @param add.corr Logical. Add Pearson and Spearman correlation
#' coefficients (top-left corner of the first facet).
#' @param palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#' 
#' @examples
#' gs <- list(
#'   Bcells = c("MS4A1","CD79B","CD79A","IGH1","IGH2"),
#'   Tcells = c("CD3E","CD3D","CD3G","CD7","CD8A")
#' )
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'
#' scatterEnrichment(
#'   pbmc,
#'   assay     = "escape",
#'   x.axis    = "Tcells",
#'   y.axis    = "Bcells",
#'   color.by  = "group",        
#'   group.by  = "groups",
#'   add.corr  = TRUE,
#'   point.size = 1
#' )
#'
#' @return A \pkg{ggplot2} object.
#' @importFrom stats as.formula
#' @export
scatterEnrichment <- function(input.data,
                              assay      = NULL,
                              x.axis,
                              y.axis,
                              facet.by   = NULL,
                              group.by   = NULL,
                              color.by   = c("density", "group", "x", "y"),
                              style      = c("point", "hex"),
                              scale      = FALSE,
                              bins       = 40,
                              point.size = 1.2,
                              alpha      = 0.8,
                              palette    = "inferno",
                              add.corr   = FALSE) {
  
  ## ---- 0  Argument sanity checks -------------------------------------------
  style <- match.arg(style, choices = c("point", "hex"))
  color.by <- match.arg(color.by, choices = c("density", "group", "x", "y"))
  if (is.null(group.by)) group.by <- "ident"
  gene.set <- c(x.axis, y.axis)
  
  ## ---- 1  Assemble long data-frame -----------------------------------------
  enriched <- .prepData(input.data, assay, gene.set, group.by, NULL, facet.by,
                        color.by = NULL)
  
  if (scale) {
    enriched[, gene.set] <- apply(enriched[, gene.set, drop = FALSE], 2, scale)
  }
  
  ## ---- 2  Base ggplot2 object ----------------------------------------------
  aes_base <- ggplot2::aes(x = .data[[x.axis]], y = .data[[y.axis]])
  
  ## ---- 3  Choose colouring strategy ----------------------------------------
  
  if (color.by == "density") {
    aes_combined <- aes_base  # no color aesthetic
  } else if (color.by == "group") {
    aes_combined <- ggplot2::aes(
      x = .data[[x.axis]], 
      y = .data[[y.axis]], 
      color = .data[[group.by]]
    )
  } else {  # "x" or "y"
    sel <- if (color.by == "x") x.axis else y.axis
    aes_combined <- ggplot2::aes(
      x = .data[[x.axis]], 
      y = .data[[y.axis]], 
      color = .data[[sel]]
    )
  }
  
  # Now build the plot
  plt <- ggplot2::ggplot(enriched, aes_combined)
  
  ## ---- 4  Geometry ---------------------------------------------------------
  if (style == "point") {
    if (color.by == "density") {
      plt <- plt +
        ggpointdensity::geom_pointdensity(size = point.size, alpha = alpha) +
        ggplot2::scale_color_gradientn(
          colors = .colorizer(palette, 11),
          name   = "Local density")
    } else {
      geom <- ggplot2::geom_point(size = point.size, alpha = alpha)
      plt  <- plt + geom
    }
  } else {                                # hex-bin
    plt <- plt +
      ggplot2::stat_binhex(bins = bins, alpha = alpha) +
      ggplot2::scale_fill_gradientn(
        colors = .colorizer(palette, 11),
        name   = "Cells / bin")
  }
  
  ## ---- 5  Colour scaling for non-density modes -----------------------------
  if (color.by != "density") {
    sel <- switch(color.by,
                  group = group.by,
                  x     = x.axis,
                  y     = y.axis)
    
    plt <- .colorby(enriched, plt,
                    color.by = sel,
                    palette  = palette,
                    type     = "color")
  }
  
  ## ---- 6  Axes, theme, faceting -------------------------------------------
  plt <- plt +
    ggplot2::labs(x = paste0(x.axis, "\nEnrichment score"),
                  y = paste0(y.axis, "\nEnrichment score")) +
    ggplot2::theme_classic()
  
  if (!is.null(facet.by)) {
    plt <- plt + ggplot2::facet_grid(as.formula(paste(". ~", facet.by)))
  }
  
  ## ---- 7  Optional correlation overlay -------------------------------------
  if (add.corr) {
    cor_pears   <- stats::cor(enriched[[x.axis]], enriched[[y.axis]],
                              method = "pearson", use = "pairwise.complete.obs")
    cor_spear   <- stats::cor(enriched[[x.axis]], enriched[[y.axis]],
                              method = "spearman", use = "pairwise.complete.obs")
    lbl <- sprintf("Pearson rho = %.2f\nSpearman rho = %.2f", cor_pears, cor_spear)
    plt <- plt +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = lbl,
                        hjust = 0, vjust = 1, size = 3.5,
                        fontface = "italic")
  }
  
  plt
}
