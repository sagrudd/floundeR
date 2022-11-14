context("Genbank File Parsing")


TB_reference = flnDr("TB_H37Rv.gb.gz")
tb <- GenbankGenome$new(TB_reference)
recA <- tb$get_cds("recA")

test_that("TB_H37Rv", {
  expect_equal(tb$accession, "AL123456 BX842572-BX842584")
  expect_equal(tb$version, "AL123456.3")
  expect_equal(tb$definition, "Mycobacterium tuberculosis H37Rv complete genome.")
})

test_that("listcds", {
  cds <- tb$list_cds()
  testthat::expect_s4_class(cds, "GenomicRanges")
  testthat::expect_s4_class(cds, "GRanges")

  tibble <- tb$as_tibble()
  testthat::expect_s3_class(tibble, "tbl")
  testthat::expect_s3_class(tibble, "data.frame")
})

test_that("single cds", {
  # there should be a single gene
  expect_equal(length(recA), 1)
  testthat::expect_s4_class(recA, "GenomicRanges")
})

test_that("multi cds", {
  multi <- tb$get_cds(c("rpoB", "rrs", "inh", "recA"))
  # there should be a single gene
  expect_equal(length(multi), 2) # a couple of the genes above are rRNAs ...
  testthat::expect_s4_class(recA, "GenomicRanges")
})
