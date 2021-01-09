context("FASTA file testing")

test_that("single.fast5", {
  
  canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
  fasta <- Fasta$new(canonical_fasta)
  expect_type(fasta, "environment")
  expect_true(all(c("Fasta", "FloundeR", "R6") %in% class(fasta)))
  expect_equal(fasta$count(), 18656)
  expect_equal(fasta$sequence_chunks(chunk_size=250), 75)
}
)