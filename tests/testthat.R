# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(escape)

skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
skip_if_not_installed("Seurat")        

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
})


test_check("escape")
