
context("Demultiplex object creation")
test_that("SequencingSet creation", {
  
  
  library(floundeR)
  seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))
  expect_equal(class(seqsum$sequencingset)[1], "SequencingSet")
  
  expect_error(
    seqsum$demultiplex, 
    "MultiplexSet requires a `barcode_arrangement` column")
  

  seqsum2 <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"),
                                   flnDr("barcoding_summary.txt.bz2"))
  expect_error(seqsum2$demultiplex, NA)
})
  
