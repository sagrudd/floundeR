# Porkchop library-preparation QC evidence

These functions call the curated `../porkchop` Rust library in-process
and return R-native evidence tables for nanopore library-preparation QC,
review, reporting, and synoptikon handoff. They do not call the Porkchop
CLI, trim reads, demultiplex reads, or write preprocessing outputs.
Public default builds keep Porkchop unlinked until that dependency is
publicly available; in that mode these wrappers raise typed R
conditions.

## Usage

``` r
library_kit_candidates(reads, read_ids = NULL)

library_adapter_primer_evidence(reads, kit_id, read_ids = NULL)

library_barcode_evidence(reads, kit_id, read_ids = NULL)

library_cdna_primer_evidence(reads, kit_id, read_ids = NULL)
```

## Arguments

- reads:

  Character vector of read sequences, or a data frame with `read_id` and
  `sequence` columns.

- read_ids:

  Optional character vector of read identifiers. Ignored when `reads` is
  a data frame with a `read_id` column.

- kit_id:

  Character scalar. Porkchop kit identifier, such as `LSK114`,
  `NBD114.24`, or `PCS114`.

## Value

`library_kit_candidates()` returns a data frame with kit metadata,
provenance, support-level, validation-status, score, `normalized_score`,
and `score_kind` columns.

`library_adapter_primer_evidence()` and `library_barcode_evidence()`
return one row per motif hit.

`library_cdna_primer_evidence()` returns one row per read with cDNA
class and primer-pair evidence columns.

## Details

Porkchop scores are heuristic evidence scores, not calibrated
probabilities. The `score_kind` column is therefore always
`"heuristic_evidence_score"`. Kit registry support, lifecycle status,
support level, provenance, validation status, and known limitations are
carried into the returned tables where applicable.

## Examples

``` r
if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
  library_kit_candidates(c("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"))
}

if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
  library_adapter_primer_evidence(
    reads = c("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"),
    kit_id = "LSK114"
  )
}

if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
  library_barcode_evidence(
    reads = c("AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT"),
    kit_id = "NBD114.24"
  )
}

if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
  library_cdna_primer_evidence(
    reads = c(paste0(
      "TTTCTGTTGGTGCTGATATTGCTTTAAAATTAAAATTAAAATTAAAATTTGGG",
      "AAAA",
      "CTTGCCTGTCGCTCTATCTTCAGAGGAG"
    )),
    kit_id = "PCS114"
  )
}
```
