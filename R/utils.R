"%!in%" <- Negate("%in%")

is_seurat_object <- function(obj) inherits(obj, "Seurat")
is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
is_seurat_or_se_object <- function(obj) {
  is_seurat_object(obj) || is_se_object(obj)
}

.checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#' @importFrom dplyr group_by summarise_at
#' @importFrom stringr str_sort
.orderFunction <- function(dat, order.by, group.by){
  if(order.by %!in% c("mean", "group.by")) {
    stop(paste0("Please select either 'mean' or 'group.by' for ordering."))
  }
  if(order.by == "mean") {
    summary <- dat %>%
                  group_by(dat[,group.by]) %>%
                  summarise_at(.vars = colnames(.)[1], mean) %>%
                  as.data.frame()
    summary <- summary[order(summary[,2], decreasing = TRUE),]
    dat[,group.by] <- factor(dat[,group.by], levels = summary[,1])
  }
  else if (order.by == "group.by") {
    dat[,group.by] <- factor(dat[,group.by], str_sort(unique(dat[,group.by]), numeric = TRUE))
  }
  return(dat)
}

.makeDFfromSCO <- function(input.data, 
                           assay = "escape", 
                           gene.set = NULL,
                           group.by = NULL, 
                           split.by = NULL, 
                           facet.by = NULL) {
  if(is.null(assay)){
    stop("Please add the assay name in which to plot from")
  }
  columns <- unique(c(group.by, split.by, facet.by))
  cnts <- .cntEval(input.data, 
                   assay = assay, 
                   type = "data")
  if(length(gene.set) == 1 && gene.set == "all") {
    gene.set <- rownames(cnts)
  }
  meta <- .grabMeta(input.data)
  if(length(gene.set) == 1) {
    enriched <- data.frame(cnts[gene.set,], meta[,columns])
  } else {
    enriched <- data.frame(t(cnts[gene.set,]), meta[,columns])
  }
  colnames(enriched) <- c(gene.set, columns)
  return(enriched)
}

#Prepare Data
.prepData <- function(input.data, assay, gene.set, group.by, split.by, facet.by) {
  
  if (inherits(x=input.data, what ="Seurat") || 
      inherits(x=input.data, what ="SummarizedExperiment")) {
    enriched <- .makeDFfromSCO(input.data, assay, gene.set, group.by, split.by, facet.by)
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(enriched)[colnames(enriched) %!in% c(group.by, split.by, facet.by)]
      gene.set <- gene.set[!grepl("meta", gene.set)]
    }
  } else if (!is_seurat_or_se_object(input.data)) {
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(input.data)
      gene.set <- gene.set[gene.set %!in% c(group.by, split.by, facet.by)]
    } 
      enriched <- data.frame(input.data[,c(gene.set,group.by, split.by, facet.by)])
    }
    
  colnames(enriched) <- c(gene.set, group.by, split.by, facet.by)
  return(enriched)
}

#' @importFrom stringr str_sort
.colorby <- function(enriched,
                     plot, 
                     color.by,
                     palette, 
                     type = "fill") { 
  if (inherits(enriched[,color.by], c("factor", "character"))) {
    grouping <- str_sort(unique(enriched[,color.by]), numeric = TRUE)
  }

  if(type == "fill") {
    if(inherits(enriched[,color.by], "numeric")) {
      plot <- plot +
              scale_fill_gradientn(colors = .colorizer(palette, 11)) + 
              labs(fill = color.by) 
    } else {
      col <- length(unique(enriched[,color.by]))
      col.pal <- .colorizer(palette, col)
      names(col.pal) <- grouping
      plot <- plot + 
              scale_fill_manual(values = col.pal) +
              labs(fill = color.by) 
    }
  } else if (type == "color") {
    if(inherits(enriched[,color.by], "numeric")) {
      plot <- plot +
              scale_color_gradientn(colors = .colorizer(palette, 11)) + 
              labs(color = color.by) 
    } else {
      col <- length(unique(enriched[,color.by]))
      col.pal <- .colorizer(palette, col)
      names(col.pal) <- grouping
      plot <- plot + 
              scale_color_manual(values = col.pal) +
              labs(color = color.by) 
    }
  }
  return(plot)
}


#Pulling a color palette for visualizations
#' @importFrom grDevices hcl.colors
#' @keywords internal
.colorizer <- function(palette = "inferno", 
                       n= NULL) {
  colors <- hcl.colors(n=n, palette = palette, fixup = TRUE)
  return(colors)
}

#split data matrix into cell chunks
#modified this from https://github.com/carmonalab/UCell
.split_data.matrix <- function(matrix, chunk.size = 1000) {
  ncols <- dim(matrix)[2]
  nchunks <- ceiling(ncols / chunk.size)  # Total number of chunks
  
  split.data <- vector("list", nchunks)  # Preallocate list for efficiency
  for (i in seq_len(nchunks)) {
    min <- (i - 1) * chunk.size + 1
    max <- min(i * chunk.size, ncols)
    split.data[[i]] <- matrix[, min:max, drop = FALSE]  # Ensure consistent structure
  }
  return(split.data)
}

#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom MatrixGenerics rowSums2
.cntEval <- function(obj, 
                     assay = "RNA", 
                     type = "counts") {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

#Add the values to single cell object
#' @importFrom SeuratObject CreateAssayObject CreateAssay5Object
#' @importFrom SummarizedExperiment SummarizedExperiment assays<-
#' @importFrom SingleCellExperiment altExps altExp<- 
.adding.Enrich <- function(sc, enrichment, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    if (as.numeric(substr(sc@version,1,1)) == 5) {
      new.assay <- suppressWarnings(CreateAssay5Object(
                                    data = as.matrix(t(enrichment))))
    } else {
      new.assay <- suppressWarnings(CreateAssayObject(
                                    data = as.matrix(t(enrichment))))
    }
    
    suppressWarnings(sc[[enrichment.name]] <- new.assay)
  } else if (inherits(sc, "SingleCellExperiment")) {
    altExp(sc, enrichment.name) <- SummarizedExperiment(assays = t(enrichment))
    names(assays(altExp(sc, enrichment.name))) <- enrichment.name
  }
  return(sc)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp
#' @importFrom Matrix t
.pull.Enrich <- function(sc, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    values <- Matrix::t(sc[[enrichment.name]]["data"])
  } else if (inherits(sc, "SingleCellExperiment")) {
    if(length(assays(altExp(sc))) == 1) {
      values <- t(assay(altExps(sc)[[enrichment.name]]))
    }
  }
}

#' @importFrom GSEABase geneIds
.GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  }
  egc <- gene.sets
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData 
#' @importFrom SeuratObject Idents
#' @importFrom methods slot
#' @keywords internal
.grabMeta <- function(sc) {
  if (is_seurat_object(sc)) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    colnames(meta)[length(meta)] <- "ident"
    
  } else if (is_se_object(sc)){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    colnames(meta)[clu] <- "ident"
  } else {
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.")
  }
  return(meta)
}

#' @importFrom SingleCellExperiment reducedDim 
.grabDimRed <- function(sc, dimRed) {
  if (is_seurat_object(sc)) {
    values <- c(list(PCA = sc[[dimRed]]@cell.embeddings),
                            sc[[dimRed]]@misc)
    
  } else if (is_se_object(sc)){
    values <- c(list(PCA = reducedDim(sc, dimRed)),
                   sc@metadata[c("eigen_values","contribution","rotation")])
    
  }
  return(values)
}

#function to split matrices by row
#adopted from ucells split_data.matrix
split_rows <- function (matrix, chunk.size = 1000) 
{
  nrows <- dim(matrix)[1]
  if(is.vector(matrix)){
    nrows <- length(matrix)
  }
  nchunks <- (nrows - 1)%/%chunk.size + 1
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks - 1) {
      left <- nrows - (i - 1) * chunk.size
      max <- min + round(left/2) - 1
    }
    else {
      max <- min(i * chunk.size, nrows)
    }
    split.data[[i]] <- matrix[min:max,]
    min <- max + 1
  }
  return(split.data)
}
#function to split vector
#adopted from ucells split_data.matrix
split_vector <- function (vector, chunk.size = 1000) 
{

  nrows <- length(vector)
  nchunks <- (nrows - 1)%/%chunk.size + 1
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks - 1) {
      left <- nrows - (i - 1) * chunk.size
      max <- min + round(left/2) - 1
    }
    else {
      max <- min(i * chunk.size, nrows)
    }
    split.data[[i]] <- vector[min:max]
    min <- max + 1
  }
  return(split.data)
}


