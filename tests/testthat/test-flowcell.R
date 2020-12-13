
context("passive flowcell tests")


test_that("passive flowcell tests", {
  
  expect_type(Flowcell$new(), "environment")
  fc <- Flowcell$new()
  expect_true(is.R6(fc))
  expect_equal(class(fc)[1], "Flowcell")
  
}) 

context("seqsum flowcells")
test_that("sequencing_summary based flowcells", {
  
  library(floundeR)
  seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))
  
  expect_equal(class(seqsum$flowcell)[1], "Flowcell")
  
  expect_error(seqsum$flowcell$density_data, NA)
  expect_equal(class(seqsum$flowcell$density_data)[1], "Angenieux")
  
  # SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))$flowcell$density_data$data
})
