test_that("ONT Zymo POD5 prefix accessor returns the canonical prefix", {
  expect_identical(
    floundeR::flounder_ont_zymo_pod5_prefix(),
    "zymo_fecal_2025.05/raw/PAU85136/pod5/"
  )
})

test_that("ONT open-data listing normalises S3 object metadata", {
  objects <- list(
    list(
      Key = "zymo_fecal_2025.05/raw/PAU85136/pod5/example-pass.pod5",
      Size = "47077200",
      LastModified = "2025-05-19T23:24:03.000Z",
      ETag = "\"abc123\"",
      StorageClass = "STANDARD"
    ),
    list(
      Key = "zymo_fecal_2025.05/raw/PAU85136/pod5/example-fail.pod5",
      Size = 163007608,
      LastModified = "2025-05-19T14:31:37.000Z",
      ETag = "\"def456\""
    )
  )

  listing <- floundeR:::.flounder_ont_open_data_listing(
    objects,
    bucket = "ont-open-data"
  )

  expect_s3_class(listing, "data.frame")
  expect_named(listing, c(
    "bucket",
    "key",
    "size",
    "last_modified_utc",
    "etag",
    "storage_class"
  ))
  expect_equal(nrow(listing), 2L)
  expect_identical(listing$bucket, c("ont-open-data", "ont-open-data"))
  expect_equal(listing$size, c(47077200, 163007608))
  expect_true(is.na(listing$storage_class[[2L]]))
})

test_that("ONT open-data public wrapper uses anonymous S3 listing parameters", {
  old <- Sys.getenv("AWS_NO_SIGN_REQUEST", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("AWS_NO_SIGN_REQUEST")
    } else {
      Sys.setenv(AWS_NO_SIGN_REQUEST = old)
    }
  })
  Sys.unsetenv("AWS_NO_SIGN_REQUEST")

  observed <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .flounder_ont_open_data_get_bucket = function(bucket, prefix, max, region) {
      observed$bucket <- bucket
      observed$prefix <- prefix
      observed$max <- max
      observed$region <- region
      observed$anonymous <- Sys.getenv("AWS_NO_SIGN_REQUEST", unset = "")
      list(list(
        Key = paste0(prefix, "PAU85136_pass_279c9095_68316534_8289.pod5"),
        Size = "47077200",
        LastModified = "2025-05-19T23:24:03.000Z",
        ETag = "\"pass-etag\"",
        StorageClass = "INTELLIGENT_TIERING"
      ))
    },
    .package = "floundeR"
  )

  listing <- floundeR::ont_open_data_list(max = 1)

  expect_identical(observed$bucket, "ont-open-data")
  expect_identical(observed$prefix, floundeR::flounder_ont_zymo_pod5_prefix())
  expect_identical(observed$max, 1L)
  expect_identical(observed$region, "eu-west-1")
  expect_identical(observed$anonymous, "true")
  expect_identical(Sys.getenv("AWS_NO_SIGN_REQUEST", unset = ""), "")
  expect_equal(listing$size, 47077200)
})

test_that("ONT open-data fetch downloads one explicit object into cache", {
  key <- paste0(
    floundeR::flounder_ont_zymo_pod5_prefix(),
    "PAU85136_pass_279c9095_68316534_8289.pod5"
  )
  cache_dir <- tempfile("flounder-ont-cache-")
  observed <- new.env(parent = emptyenv())
  observed$saves <- 0L

  testthat::local_mocked_bindings(
    .flounder_ont_open_data_get_bucket = function(bucket, prefix, max, region) {
      observed$bucket <- bucket
      observed$prefix <- prefix
      observed$max <- max
      observed$region <- region
      observed$anonymous <- Sys.getenv("AWS_NO_SIGN_REQUEST", unset = "")
      list(list(
        Key = key,
        Size = "47077200",
        LastModified = "2025-05-19T23:24:03.000Z",
        ETag = "\"pass-etag\"",
        StorageClass = "INTELLIGENT_TIERING"
      ))
    },
    .flounder_ont_open_data_save_object = function(key, bucket, file, region, overwrite) {
      observed$saves <- observed$saves + 1L
      observed$saved_key <- key
      observed$saved_bucket <- bucket
      observed$saved_region <- region
      observed$saved_overwrite <- overwrite
      dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
      writeBin(charToRaw("pod5 placeholder"), file)
      TRUE
    },
    .package = "floundeR"
  )

  fetched <- floundeR::ont_open_data_fetch(key = key, cache_dir = cache_dir)

  expect_identical(observed$bucket, "ont-open-data")
  expect_identical(observed$prefix, key)
  expect_identical(observed$max, 1L)
  expect_identical(observed$region, "eu-west-1")
  expect_identical(observed$anonymous, "true")
  expect_identical(observed$saves, 1L)
  expect_identical(observed$saved_key, key)
  expect_true(file.exists(fetched$cache_path))
  expect_identical(fetched$downloaded, TRUE)
  expect_equal(fetched$size, 47077200)
})

test_that("ONT open-data fetch reuses cached object when overwrite is false", {
  key <- paste0(
    floundeR::flounder_ont_zymo_pod5_prefix(),
    "PAU85136_pass_279c9095_68316534_8289.pod5"
  )
  cache_dir <- tempfile("flounder-ont-cache-")
  dir.create(cache_dir)
  cache_path <- file.path(cache_dir, basename(key))
  writeBin(charToRaw("already cached"), cache_path)
  observed <- new.env(parent = emptyenv())
  observed$saves <- 0L

  testthat::local_mocked_bindings(
    .flounder_ont_open_data_get_bucket = function(bucket, prefix, max, region) {
      list(list(
        Key = key,
        Size = "47077200",
        LastModified = "2025-05-19T23:24:03.000Z",
        ETag = "\"pass-etag\"",
        StorageClass = "INTELLIGENT_TIERING"
      ))
    },
    .flounder_ont_open_data_save_object = function(...) {
      observed$saves <- observed$saves + 1L
      TRUE
    },
    .package = "floundeR"
  )

  fetched <- floundeR::ont_open_data_fetch(key = key, cache_dir = cache_dir)

  expect_identical(observed$saves, 0L)
  expect_identical(fetched$cache_path, normalizePath(cache_path, mustWork = FALSE))
  expect_identical(fetched$downloaded, FALSE)
})

test_that("ONT open-data fetch rejects unsafe or missing selections", {
  expect_error(
    floundeR::ont_open_data_fetch(key = ""),
    class = "floundeR_open_data_error"
  )
  expect_error(
    floundeR::ont_open_data_fetch(
      key = paste0(floundeR::flounder_ont_zymo_pod5_prefix(), "example.pod5"),
      file_name = "../example.pod5"
    ),
    class = "floundeR_open_data_error"
  )
})

test_that("ONT open-data listing rejects invalid inputs before network access", {
  expect_error(
    floundeR::ont_open_data_list(prefix = NA_character_),
    class = "floundeR_open_data_error"
  )
  expect_error(
    floundeR::ont_open_data_list(max = 0),
    class = "floundeR_open_data_error"
  )
  expect_error(
    floundeR::ont_open_data_list(anonymous = NA),
    class = "floundeR_open_data_error"
  )
})
