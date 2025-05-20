# 2.5.0  (2025-05-19)

## Highlights
* **Streamlined code-base** – major internal refactor for clarity, speed and a ~20 % smaller dependency tree.
* **Consistent, flexible visualisation API** across all plotting helpers.
* **Robust unit-test suite** (>250 expectations) now ships with the package.

## New & enhanced functionality
| Area | Function(s) | What changed |
|------|-------------|--------------|
| **Visualisation** | `ridgeEnrichment()` | *True gradient* coloring mode for numeric `color.by`; optional per-cell rugs; quantile median line; fixed grey-fill bug |
| | `densityEnrichment()` | accepts new `rug.height`; ~4× faster ranking routine using `MatrixGenerics::rowMeans2`; cleaner two-panel layout via **patchwork** |
| | `gseaEnrichment()` | new `rug.height`; clearer legend showing ES/NES/ *p*; internal vectorised ES calculation |
| | `splitEnrichment()` | rewritten: split violins when `split.by` has 2 levels, dodged violins otherwise; inline boxplots; auto Z-scaling; palette helper |
| | `scatterEnrichment()` | density-aware points (via **ggpointdensity**), hex-bin alternative, optional Pearson/Spearman overlay, continuous or discrete color mapping |
| **Dimensionality reduction** | `performPCA()` / `pcaEnrichment()` | uses `irlba::prcomp_irlba()` automatically for large matrices; stores eigen-values/contribution in `misc`; `add.percent.contribution` now always respected |
| **Scoring backend** | `escape.matrix()` / `.compute_enrichment()` | lazy loading of heavy back-ends (*GSVA*, *UCell*, *AUCell*); unified `.build_gsva_param()`; drops empty gene-sets up-front |
| **Normalization** | `performNormalization()` | chunk-wise expressed-gene scaling (memory-friendly); accepts external `scale.factor`; optional signed log-transform; returns object with assay `<assay>_normalized` |
| **Gene-set retrieval** | `getGeneSets()` | downloads now cached under `tools::R_user_dir("escape", "cache")`; graceful KEGG append; clearer error for non-human/mouse requests |

## Performance & dependency reductions
* Replaced *plyr*, *stringr*, *rlang* usage with base-R helpers; these packages 
are now **Suggests** only.
* Common color and label utilities (`.colorizer()`, `.colorby()`, `.orderFunction()`) 
removed redundant tidyverse imports.
* Internal matrices split/chunked with new `.split_*` helpers to cap memory 
during parallel scoring/normalization.

## Bug fixes
* Gradient mode in `ridgeEnrichment()` no longer produces grey fills when the 
chosen gene-set is mapped to `color.by`.
* `pcaEnrichment()` axis labels correctly include variance contribution 
when `display.factors = FALSE`.
* `.grabDimRed()` handles both Seurat v5 and <v5 slot structures; fixes missing 
eigen-values for SCE objects.
* `escape.matrix()` respects `min.size = NULL` (no filtering) and handles 
zero-overlap gene-sets gracefully.
* Global variable declarations consolidated – eliminates *R CMD check* NOTES 
regarding `na.omit`, `value`, etc.

## Documentation 
* DESCRIPTION rewritten – heavy packages moved to *Suggests*; added explicit 
`Config/reticulate` for BiocParallel.
* `escape.gene.sets` data object now fully documented with source, usage, and reference.

# 2.4.1  (2025-03-05)
* Version bump to align with Bioconductor release cycle.
* **escape.matrix()** now silently removes gene-sets with zero detected features.

# 2.2.4  (2025-01-13)
## Underlying changes
* Switched MSigDB dependency from **msigdbr** ➜ **msigdb**.
* `getGeneSets()` gains local caching; supports only *Homo sapiens* / *Mus musculus*.

# 2.2.3  (2024-12-15)
* Fixed `groups` parameter handling and data splitting in `escape.matrix()`.
* Improved efficiency of internal `.split_data.matrix()`.

# 2.2.2  (2024-11-30)
* Patched `performNormalization()` conditional logic and per-gene-set rescaling.

# 2.2.1  (2024-11-18)
* Version bump for Bioconductor.

# 2.1.5  (2024-10-23)
* Seurat v5 compatibility; mean/median options for `heatmapEnrichment()`.

# 2.1.4 (2024-09-13)
* update ```densityEnrichment()``` GSVA function pull

# 2.1.3 (2024-09-13)

## UNDERLYING CHANGES
* update ```densityEnrichment()``` for new GSVA function name
* Parallelization of ```performNormalization()```
* Refactor of ```getGeneSets()``` to prevent issues with m_df error. 

# 2.0.1 (2024-07-26)

## UNDERLYING CHANGES
* fixed ```performNormalziation()``` errors when input.data was a matrix, now requires single-cell object and enrichment data
* passing parallel processing properly to ```runEscape()``` function.

# 1.99.1 (2024-02-29)

## UNDERLYING CHANGES
* ordering by mean values no longer changes the color order
* add explicit BPPARAM argument to ```runEscape()``` and ```escape.matrix()```
* added additional details in ```runEscape()``` and ```escape.matrix()``` for make.positive.
* removed plotting of ```splitEnrichment()``` for group.by = NULL
* separated AUC calculation to rankings and AUC, this was only method found to get consistent scores.

# 1.99.0 (2024-02-27)

## NEW FEATURES
* Added ```runEscape()```
* Added ```geyserEnrichment()```
* Added ```scatterEnrichment()```
* Added ```heatmapEnrichment()```
* Changed enrichIt to ```escape.matrix()```
* Changed enrichmentPlot to ```densityEnrichment()```
* ```performPCA()``` now works with a matrix or single-cell object
* ```pcaEnrichment()``` combines biplot-like functions

## UNDERLYING CHANGES

* Updated interaction with gsva package
* Added support for GSVA calculation
* Added support for AUCell calculation
* Added support of visualizations and calculations for single-cell objects
* Modified ```getGeneSets()``` to output a list of gene set objects with reformatted names following the Seurat "-" convention

## DEPRECATED AND DEFUNCT
* Deprecate getSignificance()
* Deprecate masterPCAPlot()

# CHANGES IN VERSION 1.9.0
* Releveling version for commit to new Bioconductor release
* Removed UCell internal functions to just import the Bioconductor UCell package

# CHANGES IN VERSION 1.4.2
* Fixed masterPCAPlot top_n() call to slice_max by top.contributions.

# CHANGES IN VERSION 1.4.1
* Version number and small edits for bioconductor compliance
* Removed singscore method
* Added UCell functions internally so they are compatible with Bioconductor
* Fixed performPCA, eliminated merge call. 

# CHANGES IN VERSION 1.3.4
* Normalization for ssGSEA no longer uses the range of all gene sets, but columns, normalizing it to 0 to 1.
* Added Kruskal-Wallis test for additional support of multi-group comparison

# CHANGES IN VERSION 1.3.3
* Added wilcoxon and LR for getSignificance
* median calculated and appended to the getSignificance() function output
* ANOVA in getSignificance() returns p-values for each comparison using TukeyHSD()
* new parameter gene.sets getSignificance() to select only gene sets or subsets.
* enrichmentPlot() beta release
* Added subcategory to getGeneSets to select subsets of libraries.
* Added default coloring to enrichmentPlot()
* enrichmentPlot() now imports calculations partially from GSVA internal functions to facilitate use of C
* Filtering based on min.size now works instead of not working. 

# CHANGES IN VERSION 1.3.2
* Added removal of gene sets with less than x features parameter in enrichIt - min.size
* Added UCell and singScore support
* new parameter gene.sets in MasterPCAPlot() and performPCA() to allow for selecting specific columns and prevent using other numeric vectors in meta data. 

# CHANGES IN VERSION 1.3.1
* Aligning versions to the current bioconductor release
* Added DietSeurat() call in vignette to prevent issues
* Adding internal gene sets - escape.gene.sets
* Removed lm.fit using limma from getSignificance

# CHANGES IN VERSION 1.0.1
* Removed ggrepel, rlang, and factoextra dependencies.
* Updated Seurat package switch
* Switch the way counts are processed by first eliminating rows with 0 expression in the sparse matrix before converting to a full matrix

# CHANGES IN VERSION 0.99.9
* Changing Seurat dependency, updated vignette

# CHANGES IN VERSION 0.99.8
* Edited getSignificance ANOVA model call

# CHANGES IN VERSION 0.99.7
* Edited getSignificance fit call to match documentation

# CHANGES IN VERSION 0.99.6
* Edited match.args() in getSignificance

# CHANGES IN VERSION 0.99.5
* Edited match.args() in getSignificance

# CHANGES IN VERSION 0.99.4
* Added match.args() to getSignificance
* Changed stop() to message()
* Modified getSignficance to allow for ANOVA and T.test

# CHANGES IN VERSION 0.99.3
* Updated link in description of getGeneSets.

# CHANGES IN VERSION 0.99.2
*Fixed a parenthesis, yeah a parenthesis. (In enrichIt() call I edited for 99.1)

# CHANGES IN VERSION 0.99.1
* Removed parallel call in gsva() and added biocparallel
* Changed cores = 4 to cores = 2 in the vignette

# CHANGES IN VERSION 0.99.0
* Preparing for bioconductor submission