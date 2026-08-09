# 2.9.1

## NEW FEATURES
* **`plaid` backend**: `escape.matrix()` and `runEscape()` gain `backend = c("native", "plaid")`. Setting `backend = "plaid"` routes `ssGSEA`, `GSVA`, `UCell`, and `AUCell` to the `replaid.*` family for large speed and memory gains. `plaid` is a `Suggests`, so the native path is unaffected if it is not installed. **These are approximations, and the gap is wider than the plaid documentation suggests** — measured pooled Pearson agreement with the native scores was roughly 0.85 (ssGSEA), 0.83 (UCell, AUCell) and 0.73 (GSVA) on a 2000-gene simulation, and lower on narrower matrices. See `?escape.matrix` and the vignette before substituting one backend for the other.
* **Three new methods**: `method = "PLAID"`, `"singscore"`, and `"scSE"`. These have no native implementation and select the `plaid` backend automatically.
* **`backend.args`**: a named list carrying method-specific tuning (`alpha`, `tau`, `rowtf`, `aucMaxRank`, `rmax`, `nsmooth`, `stats`, `chunk`, `removeLog2`, `scoreMean`). Names are validated against the target function, so typos error instead of being silently swallowed. Note `backend.args$normalize` is `plaid`'s median normalization of the scores and is unrelated to escape's `normalize`.
* **`input.assay`**: choose which expression matrix to score. `"auto"` (default) keeps raw counts for the native backend and picks log-normalized values for `plaid`. `"counts"` and `"logcounts"` map to the correct layer for both Seurat and `SummarizedExperiment`-derived objects.
* **`SpatialExperiment` support** is now explicit and tested, including `spatialCoords()`, `imgData()`, and `sample_id` preservation through `runEscape()`.
* Enrichment scores carry an `escape.backend` attribute recording the engine, method, and `plaid` version; `runEscape()` stores the same record in `metadata()` (SCE) or `Misc()` (Seurat).

## BUG FIXES
* **`escape.matrix(normalize = TRUE)` failed on every `SummarizedExperiment`-derived input** with `unable to find an inherited method for function 'assay' for signature 'x = "NULL"'`. `.pull.Enrich()` used `altExp(sc)[[name]]`, which indexes the `colData` of the *first* altExp rather than selecting the named one, and returned `NULL`. Now uses `altExp(sc, name)`. This also fixes `performNormalization()` and `performPCA()` on `SingleCellExperiment` and `SpatialExperiment` objects. Reported by Cathal King, who also provided the Xenium test data.
* **Underscored gene-set names were silently dropped during normalization.** `performNormalization()` rewrote every `_` to `-` to match Seurat's feature-name coercion, which meant `HALLMARK_*`, `GO_*`, and `REACTOME_*` sets matched nothing for matrix, `SingleCellExperiment`, or `SpatialExperiment` input. Matching is now literal first with the Seurat-mangled form as a fallback, and per-set scale factors are aligned to the enrichment columns explicitly so a partial match can never divide the wrong column.
* Requesting an assay that does not exist now reports what was asked for and what is available, instead of failing with an S4 dispatch error.
* `performNormalization()` gives a supplied `enrichment.data` precedence over scores stored on the object, and warns when both are present.
* `%||%` was used but never defined; it only exists in base R >= 4.4 while the package declares R >= 4.1.
* `runEscape()` now calls `escape.matrix()` with named rather than positional arguments.

## ENHANCEMENTS
* Removed duplicate internal definitions of `.split_cols()`, `.match_summary_fun()`, and `.filter_genes()`.
* `escape.matrix()` validates `method` up front rather than failing inside the scoring switch.
* Test suite: added a toy `SpatialExperiment` fixture built at run time (no new package data), regression coverage for all three bug fixes, and `plaid` dispatch tests that run whether or not `plaid` is installed.

# 2.7.3

## BUG FIXES
* **Robust Seurat v5/v3 assay handling**: Fixed compatibility with latest SeuratObject where the `slot` argument to `GetAssayData()` is now defunct. All internal accessors (`.cntEval()`, `.pull.Enrich()`) now use the `layer` API for SeuratObject >= 5.0.0 and fall back to `slot` only for older versions.
* Updated `.adding.Enrich()` to use the SeuratObject package version (not `sc@version`) when selecting `CreateAssay5Object` vs `CreateAssayObject`, preventing mismatches on objects created across Seurat versions.
* Fixed `performPCA()` error ("Enrichment matrix must be numeric") when enrichment data is returned as a sparse matrix from the Seurat v5 layer API.

## ENHANCEMENTS
* **CI/CD**: Added `devel` branch to R-CMD-check and test-coverage workflow triggers
* **Automated releases**: Added GitHub Actions workflow to create releases from version tags (push `v*` tags to trigger)

# 2.6.2

## NEW FEATURES
* Added `.themeEscape()` internal theme function for consistent visualization styling across all plotting functions

## ENHANCEMENTS
* **Seurat v5 compatibility**: Updated `.cntEval()` to detect SeuratObject version and use `layer` argument instead of deprecated `slot` argument for SeuratObject >= 5.0.0
* **Consistent theming**: Applied unified theme styling across all visualization functions (`ridgeEnrichment()`, `splitEnrichment()`, `geyserEnrichment()`, `heatmapEnrichment()`, `scatterEnrichment()`, `pcaEnrichment()`, `densityEnrichment()`, `gseaEnrichment()`, `enrichItPlot()`)
* **Improved `densityEnrichment()`**: Added plot title showing gene set name, alphanumeric sorting of group labels, and improved rug segment styling

## DOCUMENTATION
* Reformatted roxygen2 documentation across all exported functions for consistency
* Standardized use of `\code{}`, `\itemize{}`, `\enumerate{}`, `\strong{}`, and `\emph{}` tags
* Replaced Unicode characters with ASCII equivalents for better portability

# 2.6.1 (2025-10-31)

Bioconductor Release 3.22

## BUG FIXES
* Fixed `densityEnrichment()` interaction with GSVA package through `compute.gene.cdf` — corrected boolean argument for internal CDF computation

# 2.5.5 (2025-06-11)

## NEW FEATURES
* `geyserEnrichment()` gains `summarise.by` argument to collapse data (e.g., by patient ID) before plotting
* `color.by` now accepts both metadata columns and features (other gene sets) across `geyserEnrichment()`, `heatmapEnrichment()`, `ridgeEnrichment()`, `scatterEnrichment()`, and `splitEnrichment()`

## BUG FIXES
* Fixed `gseaEnrichment()` missing attribute error and plotting issues
* Fixed `color.by` checks to correctly determine if a variable is a feature or metadata column
* Enabled `dgCMatrix` scaling input in `geyserEnrichment()` — sparse matrices were not handled before

# 2.5.4 (2025-06-05)

## BUG FIXES
* Fixed wide-to-long format conversion for `heatmapEnrichment()` (issue #160)
* Fixed `t()` calls on sparse matrices by properly handling `dgCMatrix` in internal utilities
* `rowSums2()` now ignores NAs, fixing errors when `NA` values were present in the expression matrix
* Added `DelayedMatrixStats` as a conditional dependency for `dgCMatrix` count matrices

# 2.5.3 (2025-05-19)

## Highlights
* **Streamlined code-base** -- major internal refactor for clarity, speed, and a ~20% smaller dependency tree
* **Consistent, flexible visualisation API** across all plotting helpers
* **Robust unit-test suite** (>250 expectations) now ships with the package

## New features
* Added `enrichIt()` -- rank-based GSEA wrapper via `fgsea` for differential gene expression analysis
* Added `enrichItPlot()` -- visualization helper for `enrichIt()` results using network-style gene set plots
* Added `gseaEnrichment()` -- classic GSEA enrichment plot with running enrichment score, rug marks, and NES/p-value display

## ENHANCEMENTS
| Area | Function(s) | What changed |
|------|-------------|--------------|
| **Visualisation** | `ridgeEnrichment()` | True gradient coloring for numeric `color.by`; optional per-cell rugs; quantile median line; fixed grey-fill bug |
| | `densityEnrichment()` | New `rug.height` parameter; ~4x faster ranking via `MatrixGenerics::rowMeans2`; cleaner two-panel layout via **patchwork** |
| | `gseaEnrichment()` | New `rug.height` parameter; legend shows ES/NES/p; vectorised ES calculation |
| | `splitEnrichment()` | Rewritten: split violins when `split.by` has 2 levels, dodged violins otherwise; inline boxplots; auto Z-scaling |
| | `scatterEnrichment()` | Density-aware points via **ggpointdensity**; hex-bin alternative; optional Pearson/Spearman overlay; continuous or discrete color mapping |
| **Dimensionality reduction** | `performPCA()` / `pcaEnrichment()` | Uses `irlba::prcomp_irlba()` for large matrices; stores eigen-values/contribution in `misc`; `add.percent.contribution` now always respected |
| **Scoring backend** | `escape.matrix()` / `.compute_enrichment()` | Lazy loading of heavy back-ends (GSVA, UCell, AUCell); unified `.build_gsva_param()`; drops empty gene-sets up-front |
| **Normalization** | `performNormalization()` | Chunk-wise expressed-gene scaling (memory-friendly); accepts external `scale.factor`; optional signed log-transform; returns object with assay `<assay>_normalized` |
| **Gene-set retrieval** | `getGeneSets()` | Downloads cached under `tools::R_user_dir("escape", "cache")`; graceful KEGG append; clearer error for non-human/mouse requests |

## Performance and dependency changes
* Replaced **plyr**, **stringr**, **rlang** usage with base-R helpers; these packages are now **Suggests** only
* Common color and label utilities (`.colorizer()`, `.colorby()`, `.orderFunction()`) removed redundant tidyverse imports
* Internal matrices split/chunked with new `.split_*` helpers to cap memory during parallel scoring/normalization

## BUG FIXES
* Gradient mode in `ridgeEnrichment()` no longer produces grey fills when the chosen gene-set is mapped to `color.by`
* `pcaEnrichment()` axis labels correctly include variance contribution when `display.factors = FALSE`
* `.grabDimRed()` handles both Seurat v5 and <v5 slot structures; fixes missing eigen-values for SCE objects
* `escape.matrix()` respects `min.size = NULL` (no filtering) and handles zero-overlap gene-sets gracefully
* Global variable declarations consolidated -- eliminates R CMD check NOTES regarding `na.omit`, `value`, etc.

## Documentation
* DESCRIPTION rewritten -- heavy packages moved to **Suggests**; added explicit `Config/reticulate` for BiocParallel
* `escape.gene.sets` data object now fully documented with source, usage, and reference

# 2.4.1 (2025-03-05)

* Version bump to align with Bioconductor release cycle
* `escape.matrix()` now silently removes gene-sets with zero detected features

# 2.2.4 (2025-01-13)

## CHANGES
* Switched MSigDB dependency from **msigdbr** to **msigdb**
* `getGeneSets()` gains local caching; supports only *Homo sapiens* / *Mus musculus*

# 2.2.3 (2024-12-15)

## BUG FIXES
* Fixed `groups` parameter handling and data splitting in `escape.matrix()`
* Fixed conditional statements in `performNormalization()` per-gene-set rescaling
* Updated internal GSVA function calls to match new GSVA function names
* Imported `Matrix::t()` explicitly to resolve sparse matrix transposition issues

# 2.2.2 (2024-11-30)

## BUG FIXES
* Patched `performNormalization()` conditional logic and per-gene-set rescaling

# 2.2.1 (2024-11-18)

* Version bump for Bioconductor

# 2.1.5 (2024-10-23)

## ENHANCEMENTS
* Seurat v5 compatibility
* Added mean/median options for `heatmapEnrichment()`

# 2.1.4 (2024-09-13)

## BUG FIXES
* Updated `densityEnrichment()` GSVA function pull to match new API

# 2.1.3 (2024-09-13)

## CHANGES
* Updated `densityEnrichment()` for new GSVA function name
* Parallelization of `performNormalization()`
* Refactored `getGeneSets()` to prevent issues with `m_df` error

# 2.0.1 (2024-07-26)

## BUG FIXES
* Fixed `performNormalization()` errors when `input.data` was a matrix; now requires single-cell object and enrichment data
* Passing parallel processing properly to `runEscape()` function

# 1.99.1 (2024-02-29)

## CHANGES
* Ordering by mean values no longer changes the color order
* Added explicit `BPPARAM` argument to `runEscape()` and `escape.matrix()`
* Added additional details in `runEscape()` and `escape.matrix()` for `make.positive`
* Removed plotting of `splitEnrichment()` for `group.by = NULL`
* Separated AUC calculation to rankings and AUC for consistent scores

# 1.99.0 (2024-02-27)

## NEW FEATURES
* Added `runEscape()` -- convenience wrapper to compute and attach enrichment as a new assay
* Added `geyserEnrichment()` -- geyser-style enrichment visualization
* Added `scatterEnrichment()` -- scatter plot for pairwise gene set comparison
* Added `heatmapEnrichment()` -- heatmap of enrichment scores across groups
* Renamed `enrichIt()` to `escape.matrix()`
* Renamed `enrichmentPlot()` to `densityEnrichment()`
* `performPCA()` now works with a matrix or single-cell object
* `pcaEnrichment()` combines biplot-like functions

## CHANGES
* Updated interaction with GSVA package
* Added support for GSVA, AUCell, and single-cell object calculations/visualizations
* Modified `getGeneSets()` to output a list of gene set objects with reformatted names following the Seurat "-" convention

## DEPRECATED AND DEFUNCT
* Deprecated `getSignificance()`
* Deprecated `masterPCAPlot()`

# 1.9.0

* Releveling version for new Bioconductor release
* Removed UCell internal functions to just import the Bioconductor UCell package

# 1.4.2

* Fixed `masterPCAPlot` `top_n()` call to `slice_max` by `top.contributions`

# 1.4.1

* Version number and small edits for Bioconductor compliance
* Removed singscore method
* Added UCell functions internally for Bioconductor compatibility
* Fixed `performPCA`, eliminated merge call

# 1.3.4

* Normalization for ssGSEA no longer uses the range of all gene sets, but columns, normalizing to 0 to 1
* Added Kruskal-Wallis test for multi-group comparison support

# 1.3.3

* Added Wilcoxon and LR for `getSignificance()`
* Median calculated and appended to `getSignificance()` output
* ANOVA in `getSignificance()` returns p-values for each comparison using `TukeyHSD()`
* New `gene.sets` parameter in `getSignificance()` to select specific gene sets
* `enrichmentPlot()` beta release
* Added `subcategory` to `getGeneSets()` for library subset selection
* `min.size` filtering now works correctly

# 1.3.2

* Added `min.size` parameter to `enrichIt()` for gene set size filtering
* Added UCell and singScore support
* New `gene.sets` parameter in `masterPCAPlot()` and `performPCA()` for column selection

# 1.3.1

* Aligning versions to current Bioconductor release
* Added `DietSeurat()` call in vignette to prevent issues
* Added internal gene sets (`escape.gene.sets`)
* Removed `lm.fit` using limma from `getSignificance()`

# 1.0.1

* Removed ggrepel, rlang, and factoextra dependencies
* Updated Seurat package switch
* Optimized count processing: eliminate zero-expression rows in sparse matrix before converting to dense

# 0.99.9

* Changed Seurat dependency, updated vignette

# 0.99.8

* Edited `getSignificance()` ANOVA model call

# 0.99.7

* Edited `getSignificance()` fit call to match documentation

# 0.99.6

* Edited `match.args()` in `getSignificance()`

# 0.99.5

* Edited `match.args()` in `getSignificance()`

# 0.99.4

* Added `match.args()` to `getSignificance()`
* Changed `stop()` to `message()`
* Modified `getSignificance()` to allow ANOVA and T-test

# 0.99.3

* Updated link in description of `getGeneSets()`

# 0.99.2

* Fixed parenthesis in `enrichIt()` call

# 0.99.1

* Removed parallel call in `gsva()` and added BiocParallel
* Changed `cores = 4` to `cores = 2` in the vignette

# 0.99.0

* Preparing for Bioconductor submission
