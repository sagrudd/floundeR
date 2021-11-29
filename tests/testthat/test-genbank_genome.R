context("Genbank File Parsing")



test_that("TB_NC_000962", {
  TB_reference = flnDr("TB_H37Rv.gb.gz")
  tb <- GenbankGenome$new(TB_reference)

  expect_equal(tb$accession, "AL123456 BX842572-BX842584")
  expect_equal(tb$version, "AL123456.3")
  expect_equal(tb$definition, "Mycobacterium tuberculosis H37Rv complete genome.")
}
)
