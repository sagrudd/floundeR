test_that("pod5_find returns an empty data frame for folders without POD5 files", {
  skip_if_no_flounder_rust("POD5 discovery")

  root <- tempfile("pod5-find-empty-")
  dir.create(root)
  writeLines("not pod5", file.path(root, "reads.fastq"))

  discovered <- floundeR::pod5_find(root)

  expect_s3_class(discovered, "data.frame")
  expect_named(discovered, c(
    "path",
    "pod5_file_count",
    "total_bytes",
    "oldest_modified_utc",
    "newest_modified_utc"
  ))
  expect_equal(nrow(discovered), 0L)
})

test_that("pod5_find reports direct POD5 files by containing directory", {
  skip_if_no_flounder_rust("POD5 discovery")

  root <- tempfile("pod5-find-direct-")
  dir.create(root)
  writeBin(charToRaw("abc"), file.path(root, "reads-a.pod5"))
  writeBin(charToRaw("abcd"), file.path(root, "reads-b.POD5"))
  writeBin(charToRaw("ignored"), file.path(root, "reads-bam.bam"))

  discovered <- floundeR::pod5_find(root)

  expect_equal(nrow(discovered), 1L)
  expect_identical(discovered$path, root)
  expect_identical(discovered$pod5_file_count, 2L)
  expect_equal(discovered$total_bytes, 7)
  expect_type(discovered$oldest_modified_utc, "character")
  expect_type(discovered$newest_modified_utc, "character")
})

test_that("pod5_find reports nested POD5-containing folders separately", {
  skip_if_no_flounder_rust("POD5 discovery")

  root <- tempfile("pod5-find-nested-")
  dir.create(root)
  nested <- file.path(root, "sample-a")
  dir.create(nested)
  writeBin(charToRaw("abc"), file.path(root, "root.pod5"))
  writeBin(charToRaw("abcd"), file.path(nested, "nested.pod5"))

  discovered <- floundeR::pod5_find(root)

  expect_equal(nrow(discovered), 2L)
  expect_setequal(discovered$path, c(root, nested))
  expect_equal(sum(discovered$pod5_file_count), 2L)
  expect_equal(sum(discovered$total_bytes), 7)
})

test_that("pod5_find maps Rust discovery failures to POD5 conditions", {
  skip_if_no_flounder_rust("POD5 discovery")

  root <- tempfile("pod5-find-missing-")
  dir.create(root)
  path <- file.path(root, "missing")

  expect_error(
    floundeR::pod5_find(path),
    class = "floundeR_pod5_error"
  )
})
