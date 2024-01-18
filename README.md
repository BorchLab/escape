# escape
#### Easy single cell analysis platform for enrichment

<!-- badges: start -->
 [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/escape.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/escape)
[![R-CMD-check](https://github.com/ncborcherding/escape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/escape/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/escape/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/escape?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ncborcherding.github.io/vignettes/escape_vignette.html)
<!-- badges: end -->

<img align="right" src="https://github.com/ncborcherding/escape/blob/dev/www/escape_hex.png" width="352" height="352">

### Introduction
Single-cell sequencing (SCS) is an emerging technology in the across the diverse array of biological fields. Part of the struggle with the high-resolution approach of SCS, is distilling the data down to meaningful scientific hypotheses. escape was created to bridge SCS results, either from raw counts or from the popular Seurat R package, with gene set enrichment analyses (GSEA), allowing users to simply and easily graph outputs. The package accesses the entire [Molecular Signature Database v7.0](https://www.gsea-msigdb.org/gsea/msigdb/search.jsp) and enables users to select single, multiple gene sets, and even libraries to perform enrichment analysis on. 

### Installation

```devtools::install_github("ncborcherding/escape")```

##### Most up-to-date version

```devtools::install_github("ncborcherding/escape@dev")```

### Learning To Use escape:

Vignette available [here](https://ncborcherding.github.io/vignettes/escape_vignette.html), includes 2,000 malignant and nonmalignant peripheral blood T cells from a patient with cutaenous T cell lymphoma.

### Citation 
If using escape, please cite the [article](https://www.nature.com/articles/s42003-020-01625-6): Borcherding, N., Vishwakarma, A., Voigt, A.P. et al. Mapping the immune environment in clear cell renal carcinoma by single-cell genomics. Commun Biol 4, 122 (2021). https://doi.org/10.1038/s42003-020-01625-6. 

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 
