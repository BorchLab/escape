#' Visualize enrichment results with a ridge plot
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param input.data  Enrichment output from [escape.matrix()] or
#'   a single-cell object produced by [runEscape()].
#' @param gene.set    Gene-set (column) to plot **(length 1)**.
#' @param assay       Assay name if `input.data` is a single-cell object.
#' @param group.by    Metadata column for the y-axis groups
#'   (default `"ident"` in Seurat / SCE).
#' @param color.by    Either `"group"` (use `group.by` colours) or the
#'   name of a numeric column to map to a fill gradient.
#' @param order.by    `"mean"` | `"group"` | `NULL`.  Re-orders `group.by`
#'   factor by mean score or alphanumerically.
#' @param scale       Logical.  Z-transform the selected `gene.set`.
#' @param facet.by    Optional column to facet (`. ~ facet.by`).
#' @param add.rug     Draw per-cell tick marks underneath each ridge.
#' @param palette     Palette passed to [grDevices::hcl.colors()].
#'
#' @return A [ggplot2] object.
#' @export
#'
#' @examples
#' gs <- list(
#'   B = c("MS4A1","CD79A","CD79B"),
#'   T = c("CD3D","CD3E","CD3G","CD7")
#' )
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'
#' ridgeEnrichment(pbmc, assay = "escape",
#'                 gene.set = "T",
#'                 group.by = "groups")
ridgeEnrichment <- function(input.data,
                            gene.set,
                            assay      = NULL,
                            group.by   = NULL,
                            color.by   = "group",
                            order.by   = NULL,
                            scale      = FALSE,
                            facet.by   = NULL,
                            add.rug    = FALSE,
                            palette    = "inferno")
{
  ## ---- 0  sanity -------------------------------------------------------
  if (!requireNamespace("ggridges", quietly = TRUE))
    stop("Package 'ggridges' is required for ridge plots; please install it.")
  if (length(gene.set) != 1L)
    stop("'gene.set' must be length 1.")
  if (is.null(group.by)) group.by <- "ident"
  if (identical(color.by, "group")) color.by <- group.by
  
  ## ---- 1  build long data.frame ---------------------------------------
  df <- .prepData(input.data, assay, gene.set, group.by,
                  split.by = NULL, facet.by = facet.by)
  
  ## optional scaling (Z-transform per gene-set) -------------------------
  if (scale)
    df[[gene.set]] <- as.numeric(scale(df[[gene.set]], center = TRUE))
  
  ## optional re-ordering of the y-axis factor ---------------------------
  if (!is.null(order.by))
    df <- .orderFunction(df, order.by, group.by)
  
  ## detect “gradient” mode (numeric colour mapped to x) -----------------
  gradient.mode <-
    is.numeric(df[[color.by]]) && identical(color.by, gene.set)
  
  ## ---- 2  base ggplot --------------------------------------------------
  aes_base <- ggplot2::aes(
    x    = .data[[gene.set]],
    y    = .data[[group.by]],
    fill = if (gradient.mode) ggplot2::after_stat(x) else .data[[color.by]]
  )
  
  p <- ggplot2::ggplot(df, aes_base)
  
  ## choose ridge geometry + rug -----------------------------------------
  ridge_fun <- if (gradient.mode)
    ggridges::geom_density_ridges_gradient else ggridges::geom_density_ridges
  p <- p + do.call(ridge_fun, c(
    list(
      jittered_points = add.rug,
      point_shape     = '|',
      point_size      = 2.5,
      point_alpha     = 1,
      alpha           = 0.8,
      quantile_lines  = TRUE,
      quantile_fun    = median,
      vline_width     = 0.9
    ),
    if (add.rug) list(
      position = ggridges::position_points_jitter(width = 0.05, height = 0)
    )
  ))
  
  ## ---- 3  scales & labels ---------------------------------------------
  p <- p +
    ylab(group.by) +
    xlab(paste0(gene.set, "\nEnrichment Score")) +
    ggplot2::theme_classic(base_size = 11)
  
  p <- .colorby(df, p, color.by, palette, type = "fill")
  
  ## facetting ------------------------------------------------------------
  if (!is.null(facet.by))
    p <- p + ggplot2::facet_grid(stats::as.formula(paste(". ~", facet.by)))
  
  p
}