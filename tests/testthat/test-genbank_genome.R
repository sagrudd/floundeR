context("Genbank File Parsing")



test_that("TB_NC_000962", {
  TB_reference = flnDr("NC_000962")
  tb <- GenbankGenome$new(TB_reference)

  expect_equal(tb$accession, "NC_000962")
  expect_equal(tb$version, "NC_000962.2")
  expect_equal(tb$definition, "Mycobacterium tuberculosis H37Rv chromosome, complete genome.")
}
)
