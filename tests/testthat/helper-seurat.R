# Attach Seurat packages *only when they are available*.
# If they are missing, skip all Seurat-dependent tests gracefully.

skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
skip_if_not_installed("Seurat")        # remove if you do not use Seurat proper

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
})