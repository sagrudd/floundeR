pod5_signature_fixture <- function(path, middle_bytes = as.raw(rep(0, 16))) {
  signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
  writeBin(c(signature, middle_bytes, signature), path)
  path
}

test_that("pod5_subdivide_plan returns deterministic file-count chunks", {
  skip_if_no_flounder_rust("POD5 subdivision planning")

  root <- tempfile("pod5-subdivide-file-count-")
  dir.create(root, recursive = TRUE)
  pod5_signature_fixture(file.path(root, "a.pod5"))
  pod5_signature_fixture(file.path(root, "b.pod5"))
  pod5_signature_fixture(file.path(root, "c.pod5"))

  plan <- floundeR::pod5_subdivide_plan(
    root,
    strategy = "file-count",
    files_per_chunk = 2L
  )

  expect_s3_class(plan, "data.frame")
  expect_flounder_schema(plan, "pod5_subdivide_plan")
  expect_equal(unique(plan$schema_version), 1L)
  expect_equal(unique(plan$source), root)
  expect_equal(unique(plan$strategy), "file-count")
  expect_equal(unique(plan$target), "2 file(s) per chunk")
  expect_equal(plan$chunk_index, c(1, 2))
  expect_equal(plan$chunk_label, c("chunk-0001", "chunk-0002"))
  expect_equal(plan$file_count, c(2L, 1L))
  expect_equal(plan$total_bytes, c(64, 32))
  expect_true(all(is.na(plan$read_count)))
  expect_equal(plan$relative_paths, c("a.pod5,b.pod5", "c.pod5"))
  expect_true(all(is.na(plan$warnings)))
})

test_that("pod5_subdivide_plan groups by sample label", {
  skip_if_no_flounder_rust("POD5 subdivision planning")

  root <- tempfile("pod5-subdivide-sample-label-")
  dir.create(file.path(root, "sample-a"), recursive = TRUE)
  dir.create(file.path(root, "sample-b"), recursive = TRUE)
  pod5_signature_fixture(file.path(root, "sample-a", "a.pod5"))
  pod5_signature_fixture(file.path(root, "sample-a", "b.pod5"))
  pod5_signature_fixture(file.path(root, "sample-b", "c.pod5"))

  plan <- floundeR::pod5_subdivide_plan(root, strategy = "sample-label")

  expect_flounder_schema(plan, "pod5_subdivide_plan")
  expect_equal(plan$chunk_label, c("sample-a", "sample-b"))
  expect_equal(plan$file_count, c(2L, 1L))
  expect_equal(plan$relative_paths, c(
    paste(
      file.path("sample-a", "a.pod5"),
      file.path("sample-a", "b.pod5"),
      sep = ","
    ),
    file.path("sample-b", "c.pod5")
  ))
})

test_that("pod5_subdivide_plan reports metadata-dependent placeholder strategies", {
  skip_if_no_flounder_rust("POD5 subdivision planning")

  root <- tempfile("pod5-subdivide-placeholder-")
  dir.create(root, recursive = TRUE)
  pod5_signature_fixture(file.path(root, "a.pod5"))
  pod5_signature_fixture(file.path(root, "b.pod5"))

  elapsed <- floundeR::pod5_subdivide_plan(
    root,
    strategy = "elapsed-time",
    seconds_per_chunk = 600L
  )
  reads <- floundeR::pod5_subdivide_plan(
    root,
    strategy = "read-count",
    reads_per_chunk = 100L
  )

  expect_equal(elapsed$strategy, "elapsed-time")
  expect_equal(elapsed$target, "600 second(s) per chunk")
  expect_equal(elapsed$chunk_label, "elapsed-time-unavailable")
  expect_match(elapsed$warnings, "requires acquisition timestamps")
  expect_equal(reads$strategy, "read-count")
  expect_equal(reads$target, "100 read(s) per chunk")
  expect_equal(reads$chunk_label, "read-count-unavailable")
  expect_match(reads$warnings, "requires read counts")
})

test_that("pod5_subdivide_plan validates arguments before calling Rust", {
  skip_if_no_flounder_rust("POD5 subdivision planning")

  root <- tempfile("pod5-subdivide-validation-")
  dir.create(root, recursive = TRUE)
  pod5_signature_fixture(file.path(root, "a.pod5"))

  expect_error(
    floundeR::pod5_subdivide_plan(root, files_per_chunk = 0L),
    class = "floundeR_pod5_format_error"
  )
  expect_error(
    floundeR::pod5_subdivide_plan(root, seconds_per_chunk = 1.5),
    class = "floundeR_pod5_format_error"
  )
  expect_error(
    floundeR::pod5_subdivide_plan(file.path(root, "missing")),
    class = "floundeR_pod5_path_error"
  )
})
