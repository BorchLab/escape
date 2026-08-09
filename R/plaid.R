# -----------------------------------------------------------------------------
#  PLAID BACKEND
#
#  plaid (Zito et al., Bioinformatics 2025, btaf621) reimplements several
#  single-sample enrichment scores on top of a single sparse crossprod. escape
#  exposes it as an opt-in backend rather than a replacement: the replaid.*
#  functions are fast *approximations* of the methods they are named after, not
#  drop-in numerical equivalents.
# -----------------------------------------------------------------------------

# canonical key -> label used in messages and errors
.METHOD_LABELS <- c(SSGSEA    = "ssGSEA",
                    GSVA      = "GSVA",
                    UCELL     = "UCell",
                    AUCELL    = "AUCell",
                    PLAID     = "PLAID",
                    SINGSCORE = "singscore",
                    SCSE      = "scSE")

# canonical key -> exported plaid function name
.PLAID_FUNS <- c(SSGSEA    = "replaid.ssgsea",
                 GSVA      = "replaid.gsva",
                 UCELL     = "replaid.ucell",
                 AUCELL    = "replaid.aucell",
                 SINGSCORE = "replaid.sing",
                 SCSE      = "replaid.scse",
                 PLAID     = "plaid")

# methods that exist only through plaid
.PLAID_ONLY <- c("PLAID", "SINGSCORE", "SCSE")

# methods escape can compute itself
.NATIVE_METHODS <- c("SSGSEA", "GSVA", "UCELL", "AUCELL")

# one-time-per-session message bookkeeping
.escape_env <- new.env(parent = emptyenv())

# -----------------------------------------------------------------------------
#  METHOD / BACKEND RESOLUTION
# -----------------------------------------------------------------------------
.resolve_backend <- function(method, backend = c("native", "plaid")) {
  backend <- match.arg(backend)

  if (!is.character(method) || length(method) != 1L)
    stop("`method` must be a single string.", call. = FALSE)

  key <- toupper(method)
  if (!key %in% names(.METHOD_LABELS))
    stop("Unknown `method`: '", method, "'. One of: ",
         paste(.METHOD_LABELS, collapse = ", "), ".", call. = FALSE)

  ## PLAID / singscore / scSE have no native implementation, so `backend` is
  ## not a meaningful choice for them - route to plaid and say so.
  if (key %in% .PLAID_ONLY && backend != "plaid") {
    message("method = '", .METHOD_LABELS[[key]], "' is provided by the plaid ",
            "backend; using backend = \"plaid\".")
    backend <- "plaid"
  }

  list(method = key, backend = backend)
}

# -----------------------------------------------------------------------------
#  DEPENDENCY GUARD
# -----------------------------------------------------------------------------
.require_plaid <- function() {
  if (!requireNamespace("plaid", quietly = TRUE))
    stop("plaid not installed. Install it with ",
         "BiocManager::install(\"plaid\") to use backend = \"plaid\".",
         call. = FALSE)
  invisible(TRUE)
}

.plaid_fun <- function(key) {
  .require_plaid()
  utils::getFromNamespace(.PLAID_FUNS[[key]], "plaid")
}

# Warn once per session that plaid scores are approximations.
.plaid_fidelity_note <- function(key) {
  if (isTRUE(.escape_env$plaid_noted)) return(invisible(NULL))
  .escape_env$plaid_noted <- TRUE
  message("Scoring with plaid::", .PLAID_FUNS[[key]], "(). The replaid.* ",
          "functions are fast reimplementations and are not guaranteed to ",
          "reproduce the native escape scores exactly - see ",
          "?escape.matrix for the per-method fidelity notes.")
  invisible(NULL)
}

# -----------------------------------------------------------------------------
#  DISPATCH
# -----------------------------------------------------------------------------
# expr      : genes x cells (sparse is kept sparse)
# gene_sets : named list of character vectors, already min.size-filtered
# fn        : injection point - supply a stub to unit test the dispatch logic
#             without plaid installed
# returns   : cells x gene-sets
.plaid_dispatch <- function(expr, gene_sets, method,
                            min.size     = 5,
                            backend.args = list(),
                            groups       = NULL,
                            fn           = NULL) {
  key <- toupper(method)
  if (!key %in% names(.PLAID_FUNS))
    stop("Unknown `method`: '", method, "'.", call. = FALSE)

  if (is.null(fn)) {
    fn <- .plaid_fun(key)
    .plaid_fidelity_note(key)
  }
  fmls <- names(formals(fn))

  if (!length(gene_sets))
    stop("No gene sets to score.", call. = FALSE)
  if (is.null(names(gene_sets)))
    stop("`gene.sets` must be a named list.", call. = FALSE)

  ## ---- defaults escape owns -------------------------------------------------
  ## plaid defaults to max.genes = 500, which silently drops most HALLMARK, GO
  ## and REACTOME sets. Set a real cap above anything scoreable rather than
  ## relying on a sentinel whose meaning lives inside plaid.
  cap  <- max(c(lengths(gene_sets), nrow(expr), 1L))
  args <- list()
  if ("min.genes" %in% fmls)
    args$min.genes <- max(1L, as.integer(min.size %||% 1L))
  if ("max.genes" %in% fmls)
    args$max.genes <- as.integer(cap)
  ## `chunk` is documented only on plaid(); do not guess it through replaid dots
  if (key == "PLAID" && "chunk" %in% fmls && !is.null(groups))
    args$chunk <- as.integer(groups)
  ## plaid's `assay=` is deliberately never set - escape always hands plaid a
  ## bare matrix, so plaid never has to reach into an object.

  ## ---- user overrides -------------------------------------------------------
  if (length(backend.args)) {
    if (is.null(names(backend.args)) || any(!nzchar(names(backend.args))))
      stop("`backend.args` must be a fully named list.", call. = FALSE)
    bad <- setdiff(names(backend.args), fmls)
    if (length(bad))
      stop("Unknown `backend.args` for method '", .METHOD_LABELS[[key]],
           "' (plaid::", .PLAID_FUNS[[key]], "): ",
           paste(bad, collapse = ", "), ".\n  Accepted: ",
           paste(setdiff(fmls, c("X", "matG", "...")), collapse = ", "),
           call. = FALSE)
    if ("normalize" %in% names(backend.args))
      message("`backend.args$normalize` is plaid's median normalization of ",
              "the scores. escape's own `normalize =` argument is post-hoc ",
              "drop-out scaling (performNormalization()). They are independent.")
    args <- utils::modifyList(args, backend.args, keep.null = TRUE)
  }

  ## ---- call -----------------------------------------------------------------
  out <- do.call(fn, c(list(X = expr, matG = gene_sets), args))
  out <- as.matrix(out)

  ## ---- validate + orient (plaid returns gene sets x cells) ------------------
  ## an empty result loses its (zero-length) rownames, so check it before the
  ## dimnames guard or the user gets a misleading message
  if (!nrow(out))
    stop("plaid scored none of the ", length(gene_sets), " gene set(s). ",
         "Check that identifiers in `gene.sets` match rownames of the ",
         "expression matrix, and inspect backend.args$min.genes / ",
         "backend.args$max.genes.", call. = FALSE)
  if (is.null(rownames(out)) || is.null(colnames(out)))
    stop("plaid returned a matrix without dimnames; cannot align results.",
         call. = FALSE)
  if (!setequal(colnames(out), colnames(expr)))
    stop("plaid returned ", ncol(out), " cell(s) but ", ncol(expr),
         " were supplied.", call. = FALSE)
  out <- out[, colnames(expr), drop = FALSE]

  dropped <- setdiff(names(gene_sets), rownames(out))
  if (length(dropped) == length(gene_sets))
    stop("plaid scored none of the ", length(gene_sets), " gene set(s). ",
         "Check that identifiers in `gene.sets` match rownames of the ",
         "expression matrix.", call. = FALSE)
  if (length(dropped))
    warning("plaid dropped ", length(dropped), " of ", length(gene_sets),
            " gene set(s): ", paste(utils::head(dropped, 5L), collapse = ", "),
            if (length(dropped) > 5L) ", ..." else "",
            ". Inspect backend.args$min.genes / backend.args$max.genes.",
            call. = FALSE)

  out <- out[intersect(names(gene_sets), rownames(out)), , drop = FALSE]
  t(out)                                    # cells x gene-sets
}

# Record which engine produced a score matrix. .adding.Enrich() drops
# attributes, so runEscape() also stashes this on the object itself.
.stamp_backend <- function(mat, method, backend) {
  attr(mat, "escape.backend") <- list(
    backend = backend,
    method  = unname(.METHOD_LABELS[[toupper(method)]]),
    plaid.version = if (identical(backend, "plaid") &&
                        requireNamespace("plaid", quietly = TRUE))
      as.character(utils::packageVersion("plaid")) else NA_character_
  )
  mat
}

# Persist provenance on the object, since assay containers drop attributes.
.record_backend <- function(sc, name, prov) {
  if (is.null(prov)) return(sc)
  if (.is_seurat(sc)) {
    if (requireNamespace("SeuratObject", quietly = TRUE))
      SeuratObject::Misc(sc, slot = paste0(name, "_backend")) <- prov
  } else if (.is_sce(sc)) {
    if (requireNamespace("S4Vectors", quietly = TRUE))
      S4Vectors::metadata(sc)[[paste0(name, "_backend")]] <- prov
  }
  sc
}
