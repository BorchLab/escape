# escape VERSION 2.4.1

Version bump to be consistent with Bioconductor release

## UNDERLYING CHANGES

* Automatically remove gene sets with 0 in ```escape.matrix()```

# escape VERSION 2.2.4

## UNDERLYING CHANGES

* moved dependency from msigdbr to msigdb
* ```getGeneSets()``` now locally caches the gene sets to improve speed of repeated use
* ```getGeneSets()``` now only supports Mouse or Human

# escape VERSION 2.2.3

## UNDERLYING CHANGES

* fixed handling of *groups* parameter and data splitting in ```escape.matrix()```
* improved efficiency of internal ```.split_data.matrix()```

# escape VERSION 2.2.2

## UNDERLYING CHANGES

* fix ```performNormalization()``` conditional statements
* fix ```performNormalization()``` rescaling for per gene set calculations

# escape VERSION 2.2.1

#VERSION BUMP FOR BIOCONDUCTOR

# escape VERSION 2.1.5 (2024-10-23)

* update handling of v5 Seurat versus <v5 Seurat Objects
* added mean and median visualization for ```heatmapEnrichment()```

# escape VERSION 2.1.4 (2024-09-13)

* update ```densityEnrichment()``` GSVA function pull

# escape VERSION 2.1.3 (2024-09-13)

#VERSION BUMP FOR BIOCONDUCTOR

## UNDERLYING CHANGES

* update ```densityEnrichment()``` for new GSVA function name
* Parallelization of ```performNormalization()```
* Refactor of ```getGeneSets()``` to prevent issues with m_df error. 

# escape VERSION 2.0.1 (2024-07-26)

## UNDERLYING CHANGES

* fixed ```performNormalziation()``` errors when input.data was a matrix, now requires single-cell object and enrichment data
* passing parallel processing properly to ```runEscape()``` function.

# escape VERSION 1.99.1 (2024-02-29)

## UNDERLYING CHANGES

* ordering by mean values no longer changes the color order
* add explicit BPPARAM argument to ```runEscape()``` and ```escape.matrix()```
* added additional details in ```runEscape()``` and ```escape.matrix()``` for make.positive.
* removed plotting of ```splitEnrichment()``` for group.by = NULL
* separated AUC calculation to rankings and AUC, this was only method found to get consistent scores.


# escape VERSION 1.99.0 (2024-02-27)

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