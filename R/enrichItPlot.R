#' Adaptive visualisation of enrichIt results
#'
#' @param res `data.frame` returned by [enrichIt()].
#' @param plot.type `"bar"`, `"dot"`, or `"cnet"`.
#' @param top Integer. Keep the top *n* terms **per database**
#' (ranked by adjusted *p*). Set to `Inf` to keep all.
#' @param x.measure A column in `res` mapped to the *x*-axis
#' (ignored for `"cnet"`). Default `"-log10(padj)"`.
#' @param color.measure Column mapped to color (dot plot only).
#' Default same as `x.measure`.
#' @param show.counts Logical. Annotate bar plot with the `Count` (number of genes).
#' @param palette palette Character. Any palette from \code{\link[grDevices]{hcl.pals}}.
#' @param ... Further arguments passed to **ggplot2** geoms (e.g. 
#' `alpha`, `linewidth`).
#'
#' @return A **patchwork** object (bar / dot) or **ggraph** object (cnet).
#' @export
#'
#' @examples
#' \dontrun{
#' ranks <- setNames(markers$avg_log2FC, rownames(markers))
#' gs    <- getGeneSets("Homo sapiens", library = c("H", "C2"))
#' res   <- enrichIt(ranks, gs)
#'
#' enrichItPlot(res)               
#' enrichItPlot(res, "dot", top=10) 
#' enrichItPlot(res, "cnet", top=5) 
#' }
enrichItPlot <- function(res,
                         plot.type      = c("bar", "dot", "cnet"),
                         top            = 20,
                         x.measure      = "-log10(padj)",
                         color.measure = x.measure,
                         show.counts    = TRUE,
                         palette        = "inferno",
                         ...) {
  
  stopifnot(is.data.frame(res))
  plot.type <- match.arg(plot.type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install 'ggplot2'.")
  
  ## 0  housekeeping ----------------------------------------------------
  res <- res[order(res$padj, res$pval), , drop = FALSE]
  ## use Count if present, otherwise fall back on leadingEdge length
  if (!"Count" %in% names(res))
    res$Count <- vapply(strsplit(res$leadingEdge, ";"), length, integer(1))
  
  # Convert Database to factor
  if ("Database" %in% names(res)) {
    res$Database[is.na(res$Database)] <- "Unknown"
  } else {
    res$Database <- "Unknown"
  }
  res$Database <- factor(res$Database)
  res$Term     <- with(res, reorder(pathway, -padj))
  
  res$`-log10(padj)` <- -log10(res$padj + 1e-300)
  
  ## top-n per library -------------------------------------------------------
  if (is.finite(top)) {
    res <- do.call(rbind, lapply(split(res, res$Database), head, n = top))
    res$Term <- factor(res$Term, levels = unique(res$Term))
  }
  
  ## Bar Plot
  if (plot.type == "bar") {                 
    p <- ggplot2::ggplot(res,
                         ggplot2::aes(x = .data[[x.measure]], y = .data$Term)) +
      ggplot2::geom_col(fill = .colorizer(palette, n = 1)) +
      ggplot2::labs(x = x.measure, y = NULL) +
      ggplot2::theme_classic()
    
    if (isTRUE(show.counts)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$Count,
                     x = .data[[x.measure]] + max(.data[[x.measure]])*0.02),
        hjust = 0, size = 3)
    }
    p <- p + ggplot2::coord_cartesian(clip = "off")
  ## Dot Plot  
  } else if (plot.type == "dot") {         
    
    p <- ggplot2::ggplot(res,
                         ggplot2::aes(x = .data$geneRatio, y = .data$Term,
                                      color = .data[[color.measure]],
                                      size   = .data$size*.data$geneRatio)) +
      ggplot2::geom_point(...) +
      ggplot2::scale_size_continuous(name = "Core Count") +
      ggplot2::labs(x = "geneRatio", y = NULL,
                    color = color.measure) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.box = "vertical")
    
    if (!is.null(palette))
      p <- p + ggplot2::scale_color_gradientn(colors = .colorizer(palette, 11))
  # Network Plot
  } else {                                  
    if (!requireNamespace("ggraph", quietly = TRUE))
      stop("Install 'ggraph' for the cnet option.")
    if (!requireNamespace("igraph", quietly = TRUE))
      stop("Install 'igraph' for the cnet option.")
    
    # keep leading-edge genes only -> explode rows
    le_df <- res[seq_len(top), c("Database", "pathway", "leadingEdge")]
    le_df <- within(le_df, {
      leadingEdge <- strsplit(leadingEdge, ";")
    })
    edges <- do.call(rbind, lapply(1:nrow(le_df), function(i) {
      data.frame(pathway = le_df$pathway[i],
                 gene    = le_df$leadingEdge[[i]],
                 Database = le_df$Database[i],
                 stringsAsFactors = FALSE)
    }))
    
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% res$pathway, "pathway", "gene")
    igraph::V(g)$size <- ifelse(igraph::V(g)$type == "pathway", 8, 3)
    
    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(aes(alpha = after_stat(index)), show.legend = FALSE) +
      ggraph::geom_node_point(aes(size = .data$size,
                                  color = .data$type)) +
      ggraph::geom_node_text(aes(label = .data$name),
                             repel = TRUE, size = 3,
                             vjust = 1.5, check_overlap = TRUE) +
      ggplot2::scale_color_manual(values = .colorizer(palette, n = 2)) +
      ggplot2::theme_void()
    
    
  }
  return(p)
}
