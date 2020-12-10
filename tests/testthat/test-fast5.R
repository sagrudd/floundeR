context("FAST5 file testing")

singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
notFast5 <- system.file("extdata", "not_fast5.fast5", package="floundeR")
multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")

test_that("single.fast5", {
    expect_equal(file.isFast5(dirname(singleFast5)), FALSE)
    expect_equal(file.isFast5(singleFast5), TRUE)
    expect_equal(file.isSingleFast5(singleFast5), TRUE)
    expect_equal(file.isMultiFast5(singleFast5), FALSE)
})

test_that("not_fast5.fast5", {
    expect_equal(file.isFast5(notFast5), FALSE)
})

test_that("multi.fast5", {
    expect_equal(file.isFast5(multiFast5), TRUE)
    expect_equal(file.isSingleFast5(multiFast5), FALSE)
    expect_equal(file.isMultiFast5(multiFast5), TRUE)
})

