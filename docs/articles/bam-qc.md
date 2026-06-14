# BAM QC

This vignette shows the curated BAM QC surface in floundeR. BAM parsing
and inspection belong to Bamana-backed Rust library code called in
process from R; floundeR keeps the R API focused on QC, review,
provenance, reporting, and Synoptikon handoff.

No BAM file is bundled with this vignette. The concrete BAM inspection
calls are therefore shown as a real workflow to run against local
alignments, while the offline report-card schema example remains
evaluable without external data.

``` r

library(floundeR)
#> floundeR v0.21.19
```

## Read-Only BAM Evidence

For an aligned nanopore run, start with read-only evidence calls. These
functions return R-native lists and data frames while preserving
Bamana’s distinctions between summary, verification, validation, EOF
evidence, index state, mapping state, sort evidence, and aux-tag
evidence.

``` r

bam_path <- "/path/to/alignments.bam"

bam_summary_evidence <- bam_summary(bam_path, include_mapq_hist = TRUE)
bam_verify_evidence <- bam_verify(bam_path)
bam_validation_evidence <- bam_validate(bam_path)
bam_eof_evidence <- bam_check_eof(bam_path)
bam_index_evidence <- bam_check_index(bam_path)
bam_mapping_evidence <- bam_check_map(bam_path)
bam_sort_evidence <- bam_check_sort(bam_path)
bam_rg_tag_evidence <- bam_check_tag(bam_path, "RG", full_scan = TRUE)
```

Transformation operations such as filtering, subsampling, sorting,
merging, unmapping, reheadering, or annotation are outside this
read-only vignette. When they are exposed through floundeR, they should
require explicit output paths and provenance reporting.

## BAM Report Card

[`bam_qc_report_card()`](https://sagrudd.github.io/floundeR/reference/bam_qc_report_card.md)
consumes Bamana-derived evidence and evaluates it against stable
pass/warn/fail checks. It does not call Bamana directly, which keeps the
report-card contract testable and reusable by Synoptikon and Grammateus
reporting code.

``` r

bam_card <- bam_qc_report_card(
  summary = bam_summary_evidence,
  index = bam_index_evidence,
  mapping = bam_mapping_evidence,
  sorting = bam_sort_evidence,
  tags = bam_rg_tag_evidence,
  eof = bam_eof_evidence,
  validation = bam_validation_evidence,
  provenance_anomalies = 0,
  expected_tags = c("RG")
)

bam_card
```

Without a local BAM file, the same function still shows the complete
report-card schema and marks absent evidence explicitly.

``` r

bam_card_schema <- bam_qc_report_card(expected_tags = c("RG"))
bam_card_schema[, c("check_id", "status", "observed_value", "details")]
#> # A tibble: 11 × 4
#>    check_id                      status observed_value details                  
#>    <chr>                         <chr>           <dbl> <chr>                    
#>  1 bam_mapping_fraction          warn               NA Observed value is missin…
#>  2 bam_duplicate_fraction        warn               NA Observed value is missin…
#>  3 bam_qc_fail_fraction          warn               NA Observed value is missin…
#>  4 bam_mapq_zero_fraction        warn               NA Observed value is missin…
#>  5 bam_missing_or_unusable_index warn               NA Observed value is missin…
#>  6 bam_stale_index               warn               NA Observed value is missin…
#>  7 bam_sorting_mismatch          warn               NA Observed value is missin…
#>  8 bam_missing_expected_tags     warn               NA Observed value is missin…
#>  9 bam_eof_absence               warn               NA Observed value is missin…
#> 10 bam_validation_findings       warn               NA Observed value is missin…
#> 11 bam_provenance_anomalies      warn               NA Observed value is missin…
```

## QC Handoff Shape

BAM evidence should be handed forward as sectioned QC tables rather than
as opaque log text. The report card can be attached to a Synoptikon
payload or turned into semantic Grammateus report elements when the
optional private runtime is available.

``` r

bam_payload <- as_synoptikon_qc(
  bam = list(report_card = bam_card_schema),
  report_cards = list(bam = bam_card_schema),
  input_provenance = list(
    input_id = "example-bam-qc-schema",
    kind = "bam",
    uri = "not-bundled",
    size_bytes = NA_real_,
    sha256 = NULL,
    metadata = list(package_vignette = TRUE)
  ),
  run_identity = list(
    run_id = "example-run",
    sample_id = "example-sample",
    project_id = "flounder-revival"
  ),
  payload_id = "example-bam-qc-payload",
  generated_at_utc = "2026-06-14T14:42:10Z"
)

names(bam_payload$qc_sections)
#> [1] "sequencing_summary"  "flowcell"            "barcode"            
#> [4] "pod5"                "bam"                 "library_preparation"
bam_payload$handoff
#> $consumer
#> [1] "synoptikon"
#> 
#> $contract_version
#> [1] "v1"
#> 
#> $intended_use
#> [1] "qc_ingestion" "review"       "reporting"   
#> 
#> $compatible_with_grammateus_trusted_reports
#> [1] TRUE
```

For an end-to-end handoff that combines run summary, POD5 provenance,
BAM evidence, and library-preparation evidence, see the Synoptikon
handoff vignette.
