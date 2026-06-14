# Build and write Synoptikon QC payloads

`as_synoptikon_qc()` assembles existing floundeR QC evidence into the
versioned Synoptikon payload contract installed at
`inst/schema/synoptikon-qc-payload-v1.schema.json`.
`write_synoptikon_qc()` writes that payload as JSON.

## Usage

``` r
as_synoptikon_qc(
  run_summary = NULL,
  report_cards = NULL,
  sequencing_summary = NULL,
  flowcell = NULL,
  barcode = NULL,
  pod5 = NULL,
  bam = NULL,
  library_preparation = NULL,
  input_provenance = list(),
  run_identity = list(),
  governance_context = list(),
  payload_id = NULL,
  generated_at_utc = Sys.time(),
  limitations = list(),
  intended_use = c("qc_ingestion", "review", "reporting"),
  compatible_with_grammateus_trusted_reports = TRUE
)

write_synoptikon_qc(path, ..., pretty = TRUE)
```

## Arguments

- run_summary:

  Optional
  [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md)
  one-row tibble.

- report_cards:

  Optional report-card tibble, list of report-card tibbles, or already
  wrapped report-card list.

- sequencing_summary, flowcell, barcode, pod5, bam, library_preparation:

  Optional section evidence. Each may be a data frame/tibble or a named
  list of data frames/tibbles.

- input_provenance:

  List of input provenance records.

- run_identity:

  Named list of run, sample, project, flow-cell, kit, and instrument
  identifiers.

- governance_context:

  Named list of Synoptikon/Mneion governance identifiers.

- payload_id:

  Stable payload identifier. When omitted, a timestamped identifier is
  generated.

- generated_at_utc:

  Production timestamp as POSIXct or ISO-8601 string.

- limitations:

  List of limitation records.

- intended_use:

  Character vector of Synoptikon intended-use labels.

- compatible_with_grammateus_trusted_reports:

  Logical compatibility flag.

- path:

  Output JSON file path for `write_synoptikon_qc()`.

- ...:

  Passed from `write_synoptikon_qc()` to `as_synoptikon_qc()`.

- pretty:

  Whether to pretty-print JSON.

## Value

`as_synoptikon_qc()` returns a named list compatible with the v1
Synoptikon QC payload schema. `write_synoptikon_qc()` returns `path`
invisibly after writing JSON.

## Examples

``` r
summary_file <- flnDr("sequencing_summary.txt.bz2")
payload <- as_synoptikon_qc(
  run_summary = qc_run_summary(summary_file),
  report_cards = qc_report_card(summary_file),
  payload_id = "example-qc-payload"
)
```
