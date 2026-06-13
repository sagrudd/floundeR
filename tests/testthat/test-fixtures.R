test_that("offline fixtures are available and intentionally small", {
  fixture_files <- c(
    "sequencing_summary_dorado.tsv",
    "barcoding_summary.tsv",
    "reads.fastq",
    "reads.fasta",
    "pod5_manifest.tsv"
  )

  paths <- fixture_path(fixture_files)
  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size < 10 * 1024))
})

test_that("POD5 fixture is metadata only", {
  pod5_manifest <- read.delim(
    fixture_path("pod5_manifest.tsv"),
    stringsAsFactors = FALSE
  )

  expect_named(pod5_manifest, c(
    "file_name",
    "source_bucket",
    "source_key",
    "bytes",
    "last_modified_utc",
    "read_count",
    "state",
    "sha256"
  ))
  expect_true(all(grepl("\\.pod5$", pod5_manifest$file_name)))
  expect_false(any(file.exists(fixture_path(pod5_manifest$file_name))))
})
