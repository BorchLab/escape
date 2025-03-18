# escape
#### Easy single cell analysis platform for enrichment

<!-- badges: start -->
 [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/escape.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/escape)
[![R-CMD-check](https://github.com/ncborcherding/escape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/escape/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/escape/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/escape?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ncborcherding.github.io/vignettes/escape_vignette.html)
<!-- badges: end -->

<img align="right" src="https://github.com/BorchLab/escape/blob/master/www/escape_hex.png" width="352" height="352">

### Introduction
Single-cell sequencing (SCS) is a fundamental technology in investigating a diverse array of biological fields. Part of the struggle with the high-resolution approach of SCS, is distilling the data down to meaningful scientific hypotheses. *escape* was created to bridge SCS results, either from raw counts or from popular R-based single-cell pipelines, like [Seurat](https://satijalab.org/seurat/) or [SingleCellExperiment](https://bioconductor.org/books/release/OSCA/book-contents.html#basics), with gene set enrichment analyses (GSEA). The *escape* package allows users to easily incorporate multiple methods of GSEA and offers several visualization and analysis methods. The package accesses the entire [Molecular Signature Database v7.0](https://www.gsea-msigdb.org/gsea/msigdb/search.jsp) and enables users to select single, multiple gene sets, and even libraries to perform enrichment analysis on. 

#### Methods of GSEA Available

* Gene Set Variation Analysis (GSVA) - [citation](https://pubmed.ncbi.nlm.nih.gov/23323831/)
* Single-sample Gene Set Enrichment Analysis (ssGSEA) - [citation](https://pubmed.ncbi.nlm.nih.gov/19847166/)
* AUCell - [citation](https://pubmed.ncbi.nlm.nih.gov/28991892/)
* UCell -[citation](https://pubmed.ncbi.nlm.nih.gov/34285779/) 

More information on each method is available in the *escape* manual for ```escape.matrix()``` and the citation links. If using these methods, users should cite the original works as well.

### Installation

#### Install Via GitHub

```r
devtools::install_github("BorchLab/escape")
```

#### Install via Bioconductor

[escape v2](https://www.bioconductor.org/packages/release/bioc/html/escape.html) is available for Bioconductor users. It can be installed with the following:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("escape")
```

### Learning To Use escape:

#### Basic Usage

#### Running escape.matrix()

The basic function of enrichment analysis is done using the ```escape.matrix()``` function, with the user defining the **method** and **gene.sets** to use. 

```r
#Defining Gene Set To Use:
GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
           Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))

#Using Seurat Built-In Example:
pbmc_small <- SeuratObject::pbmc_small
  
#Running Enrichment
trial.ssGSEA <- escape.matrix(pbmc_small, 
                              method = "ssGSEA", 
                              gene.sets = GS, 
                              min.size = NULL)
```

#### runEscape()

Alternatively, ```runEscape()``` will perform the enrichment calculations as above, but also automatically amend the single-cell object with the values added as an assay, which is named via the **new.assay.name** parameter. This facilitates easy downstream visualization and analysis. 

```r
pbmc_small <- runEscape(pbmc_small,
                        method = "ssGSEA", 
                        new.assay.name = "escape.ssGSEA",
                        gene.sets = GS, 
                        min.size = NULL)
```                

### Vignette

A more comprehensive vignette including the visualizations, principal component analysis and differential testing is available [here](https://www.borch.dev/uploads/screpertoire/articles/Running_Escape.html).


## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/escape/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the the Seurat pbmc_small object would 
be extremely helpful.

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/escape/issues).

#### [Pull Requests](https://github.com/BorchLab/escape/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation 
If using escape, please cite the [article](https://www.nature.com/articles/s42003-020-01625-6): Borcherding, N., Vishwakarma, A., Voigt, A.P. et al. Mapping the immune environment in clear cell renal carcinoma by single-cell genomics. Commun Biol 4, 122 (2021). 
