
test_that("SequencingSet creation", {


    expect_type(SequencingSet$new("Species", tibble::as_tibble(iris)), "environment")
    sqset <- SequencingSet$new("Species", tibble::as_tibble(iris))
    expect_true(is.R6(sqset))
    expect_equal(class(sqset)[1], "SequencingSet")

})


test_that("SequencingSet from SequencingSummary", {

    library(floundeR)
    seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))

    expect_equal(class(seqsum$sequencingset)[1], "SequencingSet")
    expect_error(seqsum$sequencingset$enumerate, NA)
})
