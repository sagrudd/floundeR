
context("Demultiplex object creation")
test_that("SequencingSet creation", {


    library(floundeR)
    seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))

    expect_error(
        seqsum$demultiplex,
        "MultiplexSet requires a `barcode_arrangement` column")

    seqsum2 <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"),
                                     flnDr("barcoding_summary.txt.bz2"))
    expect_error(seqsum2$demultiplex, NA)

    expect_equal(class(seqsum2$demultiplex)[1], "MultiplexSet")
    expect_true(tibble::is_tibble(seqsum2$demultiplex$as_tibble()))

    expect_true(tibble::is_tibble(seqsum2$demultiplex$enumerate$data))
})

