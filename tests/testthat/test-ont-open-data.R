test_that("ONT Zymo POD5 prefix accessor returns the canonical prefix", {
  expect_identical(
    floundeR::flounder_ont_zymo_pod5_prefix(),
    "zymo_fecal_2025.05/raw/PAU85136/pod5/"
  )
})

test_that("ONT Zymo POD5 dataset helper returns canonical dataset metadata", {
  dataset <- floundeR::ont_zymo_pod5_dataset()

  expect_s3_class(dataset, "data.frame")
  expect_named(dataset, c(
    "dataset_name",
    "bucket",
    "region",
    "prefix",
    "s3_uri",
    "total_pod5_objects",
    "pass_pod5_objects",
    "fail_pod5_objects",
    "total_bytes",
    "verified_utc"
  ))
  expect_equal(nrow(dataset), 1L)
  expect_identical(dataset$bucket, "ont-open-data")
  expect_identical(dataset$region, "eu-west-1")
  expect_identical(dataset$prefix, floundeR::flounder_ont_zymo_pod5_prefix())
  expect_identical(
    dataset$s3_uri,
    "s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/"
  )
  expect_identical(dataset$total_pod5_objects, 9107L)
  expect_identical(dataset$pass_pod5_objects, 8290L)
  expect_identical(dataset$fail_pod5_objects, 817L)
  expect_equal(dataset$total_bytes, 2676973535744)
})

test_that("ONT Zymo POD5 example helper returns selected pass and fail objects", {
  objects <- floundeR::ont_zymo_pod5_example_objects()

  expect_s3_class(objects, "data.frame")
  expect_named(objects, c(
    "role",
    "state",
    "file_name",
    "bucket",
    "region",
    "key",
    "s3_uri",
    "size",
    "last_modified_utc",
    "intended_use"
  ))
  expect_equal(nrow(objects), 2L)
  expect_identical(objects$role, c("pass", "fail"))
  expect_identical(objects$bucket, c("ont-open-data", "ont-open-data"))
  expect_identical(objects$region, c("eu-west-1", "eu-west-1"))
  expect_identical(
    objects$file_name,
    c(
      "PAU85136_pass_279c9095_68316534_8289.pod5",
      "PAU85136_fail_279c9095_68316534_0.pod5"
    )
  )
  expect_identical(
    objects$key,
    paste0(floundeR::flounder_ont_zymo_pod5_prefix(), objects$file_name)
  )
  expect_equal(objects$size, c(47077200, 163007608))
  expect_identical(
    objects$last_modified_utc,
    c("2025-05-19T23:24:03.000Z", "2025-05-19T14:31:37.000Z")
  )
})

test_that("ONT Zymo POD5 example helper filters by role", {
  pass <- floundeR::ont_zymo_pod5_example_objects(role = "pass")
  fail <- floundeR::ont_zymo_pod5_example_objects(role = "fail")

  expect_equal(nrow(pass), 1L)
  expect_equal(nrow(fail), 1L)
  expect_identical(pass$role, "pass")
  expect_identical(fail$role, "fail")
  expect_match(pass$key, "PAU85136_pass_279c9095_68316534_8289[.]pod5$")
  expect_match(fail$key, "PAU85136_fail_279c9095_68316534_0[.]pod5$")
  expect_error(
    floundeR::ont_zymo_pod5_example_objects(role = "other"),
    "should be one of"
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

  expect_named(fetched, c(
    "bucket",
    "key",
    "size",
    "last_modified_utc",
    "etag",
    "storage_class",
    "cache_path",
    "downloaded"
  ))
  expect_identical(observed$bucket, "ont-open-data")
  expect_identical(observed$prefix, key)
  expect_identical(observed$max, 1L)
  expect_identical(observed$region, "eu-west-1")
  expect_identical(observed$anonymous, "true")
  expect_identical(observed$saves, 1L)
  expect_identical(observed$saved_key, key)
  expect_true(file.exists(fetched$cache_path))
  expect_identical(fetched$downloaded, TRUE)
  expect_identical(fetched$key, key)
  expect_equal(fetched$size, 47077200)
  expect_identical(fetched$last_modified_utc, "2025-05-19T23:24:03.000Z")
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
