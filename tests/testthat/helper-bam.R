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

bam_fixture_without_eof <- function(path) {
  bam_fixture(path)
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  writeBin(bytes[seq_len(length(bytes) - 28L)], path)
  path
}
