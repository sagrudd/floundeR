context("FAST5 file testing")

singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
notFast5 <- system.file("extdata", "not_fast5.fast5", package="floundeR")
multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")

test_that("single.fast5", {
    expect_equal(file.isFast5(dirname(singleFast5)), FALSE)
    expect_equal(file.isFast5(singleFast5), TRUE)
    expect_equal(file.isSingleFast5(singleFast5), TRUE)
    expect_equal(file.isMultiFast5(singleFast5), FALSE)
    expect_equal(fast5.get_platform(singleFast5), "gridion")
    expect_equal(fast5.get_flowcell_id(singleFast5), "FAK42335")
    expect_equal(fast5.get_exp_start_time(singleFast5),
                 lubridate::ymd_hms("2019-01-15 15:54:24 UTC"))
    expect_equal(fast5.get_read_count(singleFast5), 1)
    expect_equal(fast5.get_sequencing_kit(singleFast5), "sqk-lsk109")
    expect_equal(fast5.get_flowcell_type(singleFast5), "flo-min106")
})

test_that("not_fast5.fast5", {
    expect_equal(file.isFast5(notFast5), FALSE)
    expect_equal(file.isSingleFast5(notFast5), FALSE)
    expect_equal(file.isSingleFast5(notFast5), FALSE)
    expect_equal(fast5.get_platform(notFast5), NULL)
    expect_equal(fast5.get_flowcell_id(notFast5), NULL)
    # expect_equal(fast5.get_exp_start_time(notFast5), "2019-01-15 15:54:24 UTC")
    expect_equal(fast5.get_read_count(notFast5), 0)
    expect_equal(fast5.get_sequencing_kit(notFast5), NULL)
    expect_equal(fast5.get_flowcell_type(notFast5), NULL)
})

test_that("multi.fast5", {
    expect_equal(file.isFast5(multiFast5), TRUE)
    expect_equal(file.isSingleFast5(multiFast5), FALSE)
    expect_equal(file.isMultiFast5(multiFast5), TRUE)
    expect_equal(fast5.get_platform(multiFast5), "gridion")
    expect_equal(fast5.get_flowcell_id(multiFast5), "FAK42335")
    expect_equal(fast5.get_exp_start_time(multiFast5),
                 lubridate::ymd_hms("2019-01-15 15:54:24 UTC"))
    expect_equal(fast5.get_read_count(multiFast5), 10)
    expect_equal(fast5.get_sequencing_kit(multiFast5), "sqk-lsk109")
    expect_equal(fast5.get_flowcell_type(multiFast5), "flo-min106")
})

test_that("prom.fast5", {
    expect_equal(file.isFast5(promFast5), TRUE)
    expect_equal(file.isSingleFast5(promFast5), FALSE)
    expect_equal(file.isMultiFast5(promFast5), TRUE)
    expect_equal(fast5.get_platform(promFast5), "promethion")
    expect_equal(is.data.frame(fast5.get_info(promFast5)), TRUE)
    expect_equal(nrow(fast5.get_info(promFast5, atomic=TRUE)), 1)
    expect_equal(nrow(fast5.get_info(promFast5, atomic=FALSE)), 25)
    expect_equal(fast5.get_flowcell_id(promFast5), "")
    expect_equal(fast5.get_exp_start_time(promFast5),
                 lubridate::ymd_hms("2018-05-08 15:53:32 UTC"))
    expect_equal(fast5.get_read_count(promFast5), 25)
    expect_equal(fast5.get_sequencing_kit(promFast5), "sqk-lsk109")
    expect_equal(fast5.get_flowcell_type(promFast5), "flo-pro001")
})

