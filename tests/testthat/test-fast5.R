context("FAST5 file testing")



test_that("single.fast5", {
    singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
    singleF5 <- Fast5$new(fast5_file=singleFast5)
    expect_type(singleF5, "environment")
    expect_true(is.R6(singleF5))
    expect_true(singleF5$is_single_fast5())
    expect_false(singleF5$is_multi_fast5())
    expect_equal(singleF5$get_platform(), "gridion")

    expect_equal(singleF5$get_flowcell_id(), "FAK42335")
    expect_equal(singleF5$get_exp_start_time(),
                 lubridate::ymd_hms("2019-01-15 15:54:24 UTC"))
    expect_equal(singleF5$get_read_count(), 1)
    expect_equal(singleF5$get_sequencing_kit(), "sqk-lsk109")
    expect_equal(singleF5$get_flowcell_type(), "flo-min106")
})

test_that("not_fast5.fast5", {
    notFast5 <- system.file("extdata", "not_fast5.fast5", package="floundeR")
    expect_error(Fast5$new(fast5_file=notFast5))
})

test_that("multi.fast5", {
    multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
    multiF5 <- Fast5$new(fast5_file=multiFast5)
    expect_type(multiF5, "environment")
    expect_true(is.R6(multiF5))
    expect_false(multiF5$is_single_fast5())
    expect_true(multiF5$is_multi_fast5())
    expect_equal(multiF5$get_platform(), "gridion")
    expect_equal(multiF5$get_flowcell_id(), "FAK42335")
    expect_equal(multiF5$get_exp_start_time(),
                 lubridate::ymd_hms("2019-01-15 15:54:24 UTC"))
    expect_equal(multiF5$get_read_count(), 10)
    expect_equal(multiF5$get_sequencing_kit(), "sqk-lsk109")
    expect_equal(multiF5$get_flowcell_type(), "flo-min106")
})

test_that("prom.fast5", {
    promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
    promF5 <- Fast5$new(fast5_file=promFast5)
    expect_type(promF5, "environment")
    expect_true(is.R6(promF5))
    expect_false(promF5$is_single_fast5())
    expect_true(promF5$is_multi_fast5())
    expect_equal(promF5$get_platform(), "promethion")
    expect_equal(is.data.frame(promF5$get_info()), TRUE)
    expect_equal(nrow(promF5$get_info(atomic=TRUE)), 1)
    expect_equal(nrow(promF5$get_info(atomic=FALSE)), 25)
    expect_equal(promF5$get_flowcell_id(), "")
    expect_equal(promF5$get_exp_start_time(),
                 lubridate::ymd_hms("2018-05-08 15:53:32 UTC"))
    expect_equal(promF5$get_read_count(), 25)
    expect_equal(promF5$get_sequencing_kit(), "sqk-lsk109")
    expect_equal(promF5$get_flowcell_type(), "flo-pro001")
})

