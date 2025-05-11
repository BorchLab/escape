# create a cache environment.
.msigdb_cache <- new.env(parent = emptyenv())

# Helper function: Retrieve (or download and cache) the msigdb object.
getMsigdbCached <- function(org, id, version) {
  cache_key <- paste(org, id, version, sep = "_")
  
  if (exists(cache_key, envir = .msigdb_cache)) {
    message("Loading msigdb object from cache")
    msigdb_obj <- get(cache_key, envir = .msigdb_cache)
  } else {
    message("Downloading msigdb object")
    msigdb_obj <- suppressMessages(getMsigdb(org = org, id = id, version = version))
    msigdb_obj <- suppressMessages(suppressWarnings(appendKEGG(msigdb_obj)))
    assign(cache_key, msigdb_obj, envir = .msigdb_cache)
  }
  return(msigdb_obj)
}

#' Get a collection of gene sets from the msigdb
#'
#' This function retrieves gene sets from msigdb and caches the downloaded object 
#' for future calls. It allows subsetting by main collection (library), 
#' subcollection, or specific gene sets, and only supports human 
#' ("Homo sapiens") and mouse ("Mus musculus").
#'
#' @param species The scientific name of the species of interest; only 
#' "Homo sapiens" or "Mus musculus" are supported.
#' @param library A character vector of main collections (e.g. "H", "C5"). 
#' If provided, only gene sets in these collections are returned.
#' @param subcategory A character vector specifying sub-collection abbreviations 
#' (e.g. "CGP", "CP:REACTOME") to further subset the gene sets.
#' @param gene.sets A character vector of specific gene set names to select. 
#' This filter is applied after other subsetting.
#' @param version The version of MSigDB to use (default "7.4").
#' @param id The gene identifier type to use (default "SYM" for gene symbols).
#'
#' @return A named list of gene identifiers for each gene set.
#'
#' @examples
#' \dontrun{
#' # Get all hallmark gene sets from human.
#' gs <- getGeneSets(species = "Homo sapiens", 
#'                   library = "H")
#'
#' # Get a subset based on main collection and subcollection.
#' gs <- getGeneSets(species = "Homo sapiens", 
#'                   library = c("C2", "C5"), 
#'                   subcategory = "GO:BP")
#' }
#' @importFrom GSEABase GeneSet GeneSetCollection geneIds
#' @importFrom msigdb getMsigdb appendKEGG
#' @importFrom stringr str_replace_all
#' @export
getGeneSets <- function(species = "Homo sapiens", 
                        library = NULL, 
                        subcategory = NULL,
                        gene.sets = NULL,
                        version = "7.4",
                        id = "SYM") {
  # Only support human and mouse.
  if (!(species %in% c("Homo sapiens", "Mus musculus"))) {
    stop("Supported species are only 'Homo sapiens' and 'Mus musculus'.")
  }
  
  # Map species name to the organism code used by msigdb.
  org <- ifelse(species == "Homo sapiens", "hs", "mm")
  
  # Retrieve the msigdb object, from cache if available.
  msigdb_obj <- getMsigdbCached(org = org, id = id, version = version)
  
  # Filter by main collection using the S4 slot:
  if (!is.null(library)) {
    msigdb_obj <- msigdb_obj[sapply(msigdb_obj, function(x) 
      toupper(x@collectionType@category) %in% toupper(library))]
  }
  
  # Filter by subcollection using the S4 slot:
  if (!is.null(subcategory)) {
    msigdb_obj <- msigdb_obj[sapply(msigdb_obj, function(x) 
      x@collectionType@subCategory %in% toupper(subcategory))]
  }
  
  # Optional filtering by specific gene set names.
  if (!is.null(gene.sets)) {
    msigdb_obj <- msigdb_obj[sapply(msigdb_obj, function(x) x@setName %in% gene.sets)]
  }
  
  if (length(msigdb_obj) == 0) {
    warning("No gene sets found for the specified parameters.")
    return(NULL)
  }
  
  # Build the gene set list.
  gs_names <- unique(sapply(msigdb_obj, function(x) x@setName))
  gene_set_list <- vector("list", length(gs_names))
  for (i in seq_along(gs_names)) {
    genes <- unique(unlist(lapply(msigdb_obj, function(x) {
      if (x@setName == gs_names[i]) {
        return(x@geneIds)
      }
    })))
    gene_set_list[[i]] <- GSEABase::GeneSet(genes, setName = gs_names[i])
  }
  
  # Create a GeneSetCollection and return as a named list.
  gsc <- GSEABase::GeneSetCollection(gene_set_list)
  mod.names <- stringr::str_replace_all(names(gsc), "_", "-")
  gene_list <- GSEABase::geneIds(gsc)
  names(gene_list) <- mod.names
  
  return(gene_list)
}
