.msigdb_cache_dir <- tools::R_user_dir("escape", "cache")
dir.create(.msigdb_cache_dir, showWarnings = FALSE, recursive = TRUE)

# Function to cache and retrieve MSigDB gene sets
.msigdb_cached <- function(org, id = "SYM", version = "7.4") {
  key <- paste(org, id, version, sep = "_")
  file_path <- file.path(.msigdb_cache_dir, paste0(key, ".rds"))
  
  if (file.exists(file_path)) {
    gs <- readRDS(file_path)
  } else {
    if (!requireNamespace("msigdb", quietly = TRUE))
      stop("Package 'msigdb' must be installed to download MSigDB resources")
    
    gs <- suppressMessages(
      msigdb::getMsigdb(org = org, id = id, version = version)
    )
    
    # Optionally append KEGG pathways, but fail gracefully
    gs <- tryCatch(
      suppressWarnings(msigdb::appendKEGG(gs)),
      error = function(e) gs
    )
    
    saveRDS(gs, file_path)
  }
  
  gs
}

#' Get a collection of gene sets from the msigdb
#'
#' This function retrieves gene sets from msigdb and caches the downloaded object 
#' for future calls. It allows subsetting by main collection (library), 
#' subcollection, or specific gene sets, and only supports human 
#' ("Homo sapiens") and mouse ("Mus musculus").
#'
#' @param species   `"Homo sapiens"` (default) or `"Mus musculus"`.
#' @param library   Optional vector of main collection codes (e.g. `"H"`, `"C5"`).
#' @param subcategory Optional vector of sub-collection codes (e.g. `"GO:BP"`).
#' @param gene.sets Optional vector of specific gene-set names.
#' @param version   MSigDB version (character, default `"7.4"`).
#' @param id        Identifier type (default `"SYM"` for symbols).
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
#'
#' @return A named `list` of character vectors (gene IDs).  
#'         If **GSEABase** is installed, the function also returns (invisibly)
#'         a `GeneSetCollection` with the same content.
#' @export
getGeneSets <- function(species      = c("Homo sapiens", "Mus musculus"),
                        library      = NULL,
                        subcategory  = NULL,
                        gene.sets    = NULL,
                        version      = "7.4",
                        id           = "SYM")
{
  species <- match.arg(species)
  org     <- if (species == "Homo sapiens") "hs" else "mm"
  
  ## download or fetch from cache ------------------------------------------------
  msig <- .msigdb_cached(org, id, version)
  
  ## helper to interrogate S4 slots without formal import ------------------------
  .get_slot_nested <- function(x, outer_slot, inner_slot) {
    outer <- methods::slot(x, outer_slot)
    methods::slot(outer, inner_slot)
  }
  
  ## apply successive filters in one pass ---------------------------------------
  keep <- rep(TRUE, length(msig))
  
  if (!is.null(library)) {
    keep <- keep & vapply(msig,
                          \(x) toupper(.get_slot_nested(x, "collectionType", "category")),
                          "", USE.NAMES = FALSE) %in% toupper(library)
  }
  
  if (!is.null(subcategory)) {
    keep <- keep & vapply(msig,
                          function(x) {
                            ct <- methods::slot(x, "collectionType")
                            toupper(methods::slot(ct, "subCategory"))
                          },
                          "", USE.NAMES = FALSE) %in% toupper(subcategory)
  }
  
  if (!is.null(gene.sets)) {
    keep <- keep & vapply(msig, \(x) x@setName, "", USE.NAMES = FALSE) %in% gene.sets
  }
  
  msig <- msig[keep]
  if (!length(msig)) {
    warning("No gene sets matched the requested filters.")
    return(NULL)
  }
  
  ## build simple list -----------------------------------------------------------
  g.list <- lapply(msig, function(x) x@geneIds)
  names(g.list) <- vapply(msig, function(x) x@setName, "", USE.NAMES = FALSE)
  names(g.list) <- gsub("_", "-", names(g.list), fixed = TRUE)
  
  ## optionally attach GeneSetCollection invisibly ------------------------------
  if (requireNamespace("GSEABase", quietly = TRUE)) {
    gsc <- GSEABase::GeneSetCollection(
      Map(GSEABase::GeneSet, g.list, setName = names(g.list))
    )
    invisible(gsc)
  }
  
  g.list
}