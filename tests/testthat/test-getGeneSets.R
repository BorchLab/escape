# test script for getGeneSets.R - testcases are NOT comprehensive!

test_that("species argument is validated", {
  expect_error(getGeneSets("Pan troglodytes"), 
               regexp =  "Homo sapiens")
})

test_that("filtering by library works", {
  gs <- getGeneSets("Homo sapiens", library = "H")
  expect_named(gs, "HALLMARK-TEST-ONE")
  expect_identical(gs[[1]], c("geneA", "geneB"))
})

test_that("filtering by sub-category works", {
  gs <- getGeneSets("Homo sapiens", subcategory = "GO:BP")
  expect_named(gs, "TEST-SET")
  expect_identical(gs[[1]], c("geneC", "geneD"))
})

test_that("filtering by explicit gene.sets works", {
  gs <- getGeneSets("Homo sapiens", gene.sets = "HALLMARK_TEST_ONE")
  expect_named(gs, "HALLMARK-TEST-ONE")
  expect_identical(gs[[1]], c("geneA", "geneB"))
})

test_that("combined filters (library + subcategory) work", {
  gs <- getGeneSets("Homo sapiens", library = "C5", subcategory = "GO:BP")
  expect_named(gs, "TEST-SET")
})

test_that("requesting an empty subset warns and returns NULL", {
  expect_warning(
    out <- getGeneSets("Homo sapiens", library = "NONEXISTENT"),
    regexp = "matched the requested filters."
  )
  expect_null(out)
})
