# ------------------------------------------------------------------ #
# 1. build a tiny fake MSigDB object --------------------------------
# ------------------------------------------------------------------ #
setClass("FakeCollectionType",
         slots = c(category = "character", subCategory = "character"))
setClass("FakeGeneSet",
         slots = c(setName = "character",
                   geneIds = "character",
                   collectionType = "FakeCollectionType"))

.fake_msigdb <- list(
  new("FakeGeneSet",
      setName        = "HALLMARK_TEST_ONE",
      geneIds        = c("geneA", "geneB"),
      collectionType = new("FakeCollectionType",
                           category    = "H",
                           subCategory = "CGP")),
  new("FakeGeneSet",
      setName        = "TEST_SET",
      geneIds        = c("geneC", "geneD"),
      collectionType = new("FakeCollectionType",
                           category    = "C5",
                           subCategory = "GO:BP"))
)

# ------------------------------------------------------------------ #
# 2. overwrite .msigdb_cached() inside the escape namespace ---------
# ------------------------------------------------------------------ #
ns <- asNamespace("escape")
unlockBinding(".msigdb_cached", ns)
assign(".msigdb_cached",
       function(org, id = "SYM", version = "7.4") .fake_msigdb,
       envir = ns)
lockBinding(".msigdb_cached", ns)