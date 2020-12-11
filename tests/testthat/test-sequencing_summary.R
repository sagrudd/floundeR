
context("Sequencing_summary testing")
test_that("sequencing_summary", {

    sequencing_summary <- system.file("extdata",
                                      "sequencing_summary.txt.bz2",
                                      package="floundeR")
    # fail with a non-string starter of length == 1
    expect_error(SequencingSummary$new())
    expect_error(SequencingSummary$new(sequencing_summary_file=123))
    # check construction with non_existent file
    expect_error(SequencingSummary$new(sequencing_summary_file="thisfile.txt"))
    # and check that the specified extant path is actually a file (not dir)
    expect_error(SequencingSummary$new(sequencing_summary_file=".."))
    expect_error(SequencingSummary$new(sequencing_summary_file=sequencing_summary), NA)
})


context("Sequencing_summary + barcoding_summary testing")
test_that("sequencing_summary + barcode", {

    sequencing_summary <- system.file("extdata",
                                      "sequencing_summary.txt.bz2",
                                      package="floundeR")
    barcodes_summary <- system.file("extdata",
                                    "barcoding_summary.txt.bz2",
                                    package="floundeR")

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
        SequencingSummary$new(
            sequencing_summary_file=sequencing_summary,
            barcoding_summary_file=barcodes_summary), NA)
})
