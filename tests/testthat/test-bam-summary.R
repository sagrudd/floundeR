bam_fixture <- function(path) {
  skip_if_not_installed("Rsamtools")

  i32 <- function(value) {
    writeBin(as.integer(value), raw(), size = 4, endian = "little")
  }
  u32 <- i32
  bam_string <- function(value) {
    c(charToRaw(value), as.raw(0))
  }
  bam_record <- function(ref_id, pos, read_name, flag, mapq = 60L) {
    name <- bam_string(read_name)
    l_read_name <- length(name)
    bin_mq_nl <- bitwOr(bitwShiftL(as.integer(mapq), 8L), l_read_name)
    flag_nc <- bitwShiftL(as.integer(flag), 16L)

    c(
      i32(32L + l_read_name),
      i32(ref_id),
      i32(pos),
      u32(bin_mq_nl),
      u32(flag_nc),
      i32(0L),
      i32(-1L),
      i32(-1L),
      i32(0L),
      name
    )
  }

  header <- charToRaw("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n")
  reference <- bam_string("chr1")
  payload <- c(
    as.raw(c(charToRaw("BAM"), 0x01)),
    i32(length(header)),
    header,
    i32(1L),
    i32(length(reference)),
    reference,
    i32(100L),
    bam_record(0L, 10L, "mapped", 0L, mapq = 60L),
    bam_record(-1L, -1L, "unmapped", 4L, mapq = 0L)
  )

  raw_path <- tempfile("flounder-bam-fixture-")
  writeBin(payload, raw_path)
  Rsamtools::bgzip(raw_path, dest = path, overwrite = TRUE)
  path
}

test_that("bam_summary returns stable R-native summary tables", {
  skip_if_no_flounder_rust("BAM summary")

  path <- tempfile("flounder-summary-", fileext = ".bam")
  bam_fixture(path)

  summary <- floundeR::bam_summary(path, include_mapq_hist = TRUE)

  expect_s3_class(summary, "floundeR_bam_summary")
  expect_named(summary, c(
    "status",
    "evidence",
    "header",
    "counts",
    "fractions",
    "fractions_observed",
    "mapq",
    "mapping",
    "anomalies",
    "flag_categories",
    "references",
    "index_derived",
    "mapq_histogram"
  ))
  expect_flounder_schema(summary$status, "bam_summary_status")
  expect_flounder_schema(summary$evidence, "bam_summary_evidence")
  expect_flounder_schema(summary$counts, "bam_summary_counts")
  expect_flounder_schema(summary$fractions, "bam_summary_fractions")
  expect_flounder_schema(summary$mapq, "bam_summary_mapq")

  expect_equal(summary$status$command, "summary")
  expect_true(summary$status$ok)
  expect_equal(summary$status$format, "BAM")
  expect_equal(summary$status$mode, "full_scan")
  expect_equal(summary$status$confidence, "high")
  expect_equal(summary$evidence$records_scanned, 2)
  expect_true(summary$evidence$full_file_scanned)
  expect_equal(summary$header$references_defined, 1L)
  expect_equal(summary$header$sort_order, "coordinate")

  expect_equal(summary$counts$records_examined, 2)
  expect_equal(summary$counts$mapped_records, 1)
  expect_equal(summary$counts$unmapped_records, 1)
  expect_equal(summary$counts$primary_records, 2)
  expect_equal(summary$counts$secondary_records, 0)
  expect_equal(summary$fractions$fraction_mapped, 0.5)
  expect_equal(summary$mapq$min, 0)
  expect_equal(summary$mapq$max, 60)
  expect_equal(summary$mapq$zero_count, 1)
  expect_equal(summary$mapping$status, "mapped")
  expect_equal(summary$references$name, "chr1")
  expect_equal(summary$references$length, 100)
  expect_equal(summary$mapq_histogram$mapq, c(0L, 60L))
  expect_equal(summary$mapq_histogram$read_count, c(1, 1))
})

test_that("bam_summary validates arguments and preserves Bamana path errors", {
  skip_if_no_flounder_rust("BAM summary")

  expect_error(
    floundeR::bam_summary(NA_character_),
    class = "floundeR_bam_path_error"
  )
  expect_error(
    floundeR::bam_summary("missing.bam", sample_records = -1L),
    class = "floundeR_bam_argument_error"
  )
  expect_error(
    floundeR::bam_summary("missing.bam"),
    class = "floundeR_bam_path_error"
  )
})
