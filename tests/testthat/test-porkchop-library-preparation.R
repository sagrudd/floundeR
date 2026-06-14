test_that("Porkchop kit candidates carry heuristic score terminology", {
  skip_if_no_flounder_porkchop("Porkchop library-preparation evidence")

  reads <- data.frame(
    read_id = c("read_la", "read_nb01"),
    sequence = c(
      "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
      "AAGGTTAA CACAAAGACACCGACAACTTTCTT CAGCACCT"
    ),
    stringsAsFactors = FALSE
  )
  reads$sequence <- gsub(" ", "", reads$sequence, fixed = TRUE)

  candidates <- floundeR::library_kit_candidates(reads)

  expect_s3_class(candidates, "data.frame")
  expect_true(all(c(
    "schema_version",
    "kit_id",
    "score",
    "normalized_score",
    "score_kind",
    "support_level",
    "validation_status",
    "known_limitations"
  ) %in% names(candidates)))
  expect_true(nrow(candidates) > 0)
  expect_true(any(candidates$total_hits > 0L))
  expect_true(all(candidates$score_kind == "heuristic_evidence_score"))
  expect_false(any(grepl("prob", names(candidates), ignore.case = TRUE)))
})

test_that("Porkchop adapter and barcode evidence wrappers return motif hits", {
  skip_if_no_flounder_porkchop("Porkchop library-preparation evidence")

  adapter <- floundeR::library_adapter_primer_evidence(
    reads = data.frame(
      read_id = "read_la",
      sequence = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
      stringsAsFactors = FALSE
    ),
    kit_id = "LSK114"
  )

  expect_true(nrow(adapter) > 0)
  expect_true(all(adapter$motif_family %in% c("adapter", "primer")))
  expect_true("LA_top" %in% adapter$motif_name)
  expect_true(all(adapter$match_semantics == "ambiguity_aware_exact"))

  barcode <- floundeR::library_barcode_evidence(
    reads = data.frame(
      read_id = "read_nb01",
      sequence = "AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT",
      stringsAsFactors = FALSE
    ),
    kit_id = "NBD114.24"
  )

  expect_true(nrow(barcode) > 0)
  expect_true(all(barcode$motif_family %in% c("barcode", "flank")))
  expect_true(any(barcode$motif_name %in% c("NB01", "NB_flank_fwd")))
})

test_that("Porkchop cDNA primer evidence preserves class vocabulary", {
  skip_if_no_flounder_porkchop("Porkchop cDNA primer evidence")

  cdna <- floundeR::library_cdna_primer_evidence(
    reads = data.frame(
      read_id = "read_cdna",
      sequence = paste0(
        "TTTCTGTTGGTGCTGATATTGCTTTAAAATTAAAATTAAAATTAAAATTTGGG",
        "AAAA",
        "CTTGCCTGTCGCTCTATCTTCAGAGGAG"
      ),
      stringsAsFactors = FALSE
    ),
    kit_id = "PCS114"
  )

  expect_equal(cdna$class, "full_length")
  expect_true(cdna$classified)
  expect_true(cdna$full_length)
  expect_equal(cdna$five_prime_name, "SSPII")
  expect_equal(cdna$three_prime_name, "RTP")
  expect_gte(cdna$primer_hit_count, 2L)
  expect_match(cdna$workflow_support_note, "primer-pair")
})

test_that("Porkchop wrappers reject unsupported cDNA kit rules", {
  skip_if_no_flounder_porkchop("Porkchop cDNA primer evidence")

  expect_error(
    floundeR::library_cdna_primer_evidence("ACGT", kit_id = "LSK114"),
    class = "floundeR_library_preparation_error"
  )
})
