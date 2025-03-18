# test script for getGeneSets.R - testcases are NOT comprehensive!

context("Testing getGeneSets and caching behavior")

# Define fake S4 classes to mimic the msigdb gene set objects.
setClass("FakeCollectionType", slots = c(category = "character", subCategory = "character"))
setClass("FakeGeneSet", 
         slots = c(setName = "character", 
                   geneIds = "character", 
                   collectionType = "FakeCollectionType"))

# Create two fake gene set objects.
fake1 <- new("FakeGeneSet",
             setName = "HALLMARK_TEST_ONE",
             geneIds = c("geneA", "geneB"),
             collectionType = new("FakeCollectionType", 
                                  category = "H", 
                                  subCategory = "CGP"))

fake2 <- new("FakeGeneSet",
             setName = "TEST_SET",
             geneIds = c("geneC", "geneD"),
             collectionType = new("FakeCollectionType", 
                                  category = "C5", 
                                  subCategory = "GO:BP"))

# Combine into a list to simulate the msigdb object.
fake_list <- list(fake1, fake2)

# Clear the package-level cache before running tests.
rm(list = ls(envir = .msigdb_cache), envir = .msigdb_cache)

# Insert the fake object into the cache for human (key: "hs_SYM_7.4").
assign("hs_SYM_7.4", fake_list, envir = .msigdb_cache)

test_that("Unsupported species throws an error", {
  expect_error(
    getGeneSets(species = "Pan troglodytes"),
    "Supported species are only 'Homo sapiens' and 'Mus musculus'."
  )
})

test_that("Filtering by library (main collection) works", {
  gs <- getGeneSets(species = "Homo sapiens", library = "H")
  # Only fake1 has library "H".
  expect_equal(names(gs), "HALLMARK-TEST-ONE")
  expect_equal(gs[["HALLMARK-TEST-ONE"]], c("geneA", "geneB"))
})

test_that("Filtering by subcategory works", {
  gs <- getGeneSets(species = "Homo sapiens", subcategory = "GO:BP")
  # Only fake2 has subcategory "GO:BP".
  expect_equal(names(gs), "TEST-SET")
  expect_equal(gs[["TEST-SET"]], c("geneC", "geneD"))
})

test_that("Filtering by specific gene.sets works", {
  gs <- getGeneSets(species = "Homo sapiens", gene.sets = "HALLMARK_TEST_ONE")
  expect_equal(names(gs), "HALLMARK-TEST-ONE")
  expect_equal(gs[["HALLMARK-TEST-ONE"]], c("geneA", "geneB"))
})

test_that("Combined filtering by library and subcategory works", {
  gs <- getGeneSets(species = "Homo sapiens", library = "C5", subcategory = "GO:BP")
  expect_equal(names(gs), "TEST-SET")
  expect_equal(gs[["TEST-SET"]], c("geneC", "geneD"))
})

test_that("No gene sets found triggers a warning and returns NULL", {
  expect_warning(
    result <- getGeneSets(species = "Homo sapiens", library = "NONEXISTENT"),
    "No gene sets found for the specified parameters."
  )
  expect_null(result)
})

test_that("Caching behavior works for a new species (Mus musculus)", {
  # Remove any existing mouse object from the cache.
  if (exists("mm_SYM_7.4", envir = .msigdb_cache)) {
    rm("mm_SYM_7.4", envir = .msigdb_cache)
  }
  
  # Capture messages on the first call (should simulate a download).
  msgs_download <- character()
  withCallingHandlers({
    getGeneSets(species = "Mus musculus", library = "H")
  }, message = function(m) {
    msgs_download <<- c(msgs_download, m$message)
    invokeRestart("muffleMessage")
  })
  expect_true(any(grepl("Downloading msigdb object", msgs_download)))
  
  # Now the mouse object should be cached.
  msgs_cache <- character()
  withCallingHandlers({
    getGeneSets(species = "Mus musculus", library = "H")
  }, message = function(m) {
    msgs_cache <<- c(msgs_cache, m$message)
    invokeRestart("muffleMessage")
  })
  expect_true(any(grepl("Loading msigdb object from cache", msgs_cache)))
})
