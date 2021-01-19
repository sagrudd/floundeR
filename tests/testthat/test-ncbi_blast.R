context("FAST5 file testing")

test_that("single.fast5", {
  blast_results <- flnDr("drosophila_uniref100.blastx.gz")
  blast <- Blast$new(blast_file=blast_results)

  expect_type(blast, "environment")
  expect_true(all(c("Blast", "FloundeR", "R6") %in% class(blast)))

  expect_equal(blast$count(), 248)
  }
)
