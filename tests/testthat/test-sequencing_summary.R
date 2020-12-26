
context("Sequencing_summary testing")
test_that("sequencing_summary", {

    sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
    # fail with a non-string starter of length == 1
    expect_error(SequencingSummary$new())
    expect_error(SequencingSummary$new(123))
    # check construction with non_existent file
    expect_error(SequencingSummary$new("thisfile.txt"))
    # and check that the specified extant path is actually a file (not dir)
    expect_error(SequencingSummary$new(".."))
    expect_error(seqsum <- SequencingSummary$new(sequencing_summary), NA)
})


context("Sequencing_summary + barcoding_summary testing")
test_that("sequencing_summary + barcode", {

    sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
    barcodes_summary <- flnDr("barcoding_summary.txt.bz2")

    # fail if barcode is provided but not seq_sum
    expect_error(SequencingSummary$new(barcoding_summary_file=barcodes_summary))
    expect_error(
        SequencingSummary$new(
            sequencing_summary_file=sequencing_summary,
            barcoding_summary_file=123))
    expect_error(
        SequencingSummary$new(
            sequencing_summary_file=sequencing_summary,
            barcoding_summary_file="missing.csv"))
    expect_error(
        SequencingSummary$new(
            sequencing_summary_file=sequencing_summary,
            barcoding_summary_file=".."))
    expect_error(
        seqsumBC <- SequencingSummary$new(
            sequencing_summary_file=sequencing_summary,
            barcoding_summary_file=barcodes_summary), NA)
})



context("Flowcell object creation")
test_that("Flowcell object creation", {
    
    sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
    barcodes_summary <- flnDr("barcoding_summary.txt.bz2")

    seqsum <- SequencingSummary$new(
        sequencing_summary_file=sequencing_summary,
        barcoding_summary_file=barcodes_summary)
    
    expect_equal(class(seqsum$flowcell)[1], "Flowcell")
    
    flowcell <- seqsum %>% to_flowcell()
    expect_equal(class(flowcell)[1], "Flowcell")
    

    # please note that Flowcell testing is continued in the Flowcell unit tests
})
