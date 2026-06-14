# Synoptikon QC Handoff

This vignette shows the package-level handoff contract from floundeR QC
data to Synoptikon/Mneion. It is an R package vignette, not the
production QC report renderer. Governed HTML/PDF report rendering
remains a Grammateus responsibility and is intentionally optional for
the open-source floundeR package.

The example is network-free. It uses the package sequencing-summary
fixture and metadata for the selected ONT open-data POD5 object without
downloading POD5 bytes.

For focused input-level tutorials, see the POD5 QC and BAM QC vignettes.
This vignette combines those evidence shapes into the Synoptikon JSON
handoff contract.

``` r

library(floundeR)
#> floundeR v0.21.19
```

## Build Run-Level QC Evidence

Start from a sequencing-summary file and derive the stable QC tables
that Synoptikon can ingest directly.

``` r

summary_file <- flnDr("sequencing_summary.txt.bz2")

run_summary <- qc_run_summary(
  summary_file,
  source_id = "example-sequencing-summary"
)

flowcell_density <- qc_channel_density(summary_file)
barcode_composition <- qc_barcode_composition(summary_file)
run_report_card <- qc_report_card(summary_file)

run_summary
#> # A tibble: 1 × 20
#>   schema_version        source_id read_count passed_read_count failed_read_count
#>   <chr>                 <chr>          <int>             <int>             <int>
#> 1 flounder.qc_run_summ… example-…      10000              8355              1645
#> # ℹ 15 more variables: pass_fraction <dbl>, total_bases <dbl>,
#> #   passed_bases <dbl>, failed_bases <dbl>, mean_read_length <dbl>,
#> #   median_read_length <dbl>, n50_read_length <dbl>, mean_qscore <dbl>,
#> #   median_qscore <dbl>, channel_count <int>, first_read_start_time <dbl>,
#> #   last_read_start_time <dbl>, run_duration_seconds <dbl>,
#> #   barcode_count <int>, unclassified_read_count <int>
run_report_card
#> # A tibble: 7 × 9
#>   schema_version       check_id check_label status observed_value warn_threshold
#>   <chr>                <chr>    <chr>       <chr>           <dbl>          <dbl>
#> 1 flounder.qc_report_… pass_fr… Pass read … warn            0.836           0.85
#> 2 flounder.qc_report_… mean_qs… Mean read … warn            9.48           10   
#> 3 flounder.qc_report_… n50_rea… Read lengt… pass        26855            1000   
#> 4 flounder.qc_report_… total_b… Total base… pass    147095512         1000000   
#> 5 flounder.qc_report_… channel… Active cha… pass          502              64   
#> 6 flounder.qc_report_… unclass… Unclassifi… pass            0               0.1 
#> 7 flounder.qc_report_… barcode… Largest ba… warn           NA               0.6 
#> # ℹ 3 more variables: fail_threshold <dbl>, comparator <chr>, details <chr>
```

## Add POD5 Source Metadata

The canonical contemporary POD5 source for floundeR demonstrations is
the ONT Zymo fecal run at
`s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/`. The routine
pass example is the smallest pass POD5 object recorded for that prefix.
This vignette records the selected object as provenance without
downloading it.

``` r

pod5_dataset <- ont_zymo_pod5_dataset()
pod5_example <- ont_zymo_pod5_example_objects(role = "pass")

pod5_handoff <- data.frame(
  schema_version = "flounder.pod5_open_data_example.v1",
  status = "not_checked",
  bucket = pod5_example$bucket,
  region = pod5_example$region,
  key = pod5_example$key,
  s3_uri = pod5_example$s3_uri,
  size_bytes = pod5_example$size,
  last_modified_utc = pod5_example$last_modified_utc,
  intended_use = pod5_example$intended_use,
  stringsAsFactors = FALSE
)

pod5_dataset
#>                           dataset_name        bucket    region
#> 1 ont_zymo_fecal_2025_05_pau85136_pod5 ont-open-data eu-west-1
#>                                  prefix
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/
#>                                                     s3_uri total_pod5_objects
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/               9107
#>   pass_pod5_objects fail_pod5_objects  total_bytes         verified_utc
#> 1              8290               817 2.676974e+12 2026-06-13T00:00:00Z
pod5_handoff
#>                       schema_version      status        bucket    region
#> 1 flounder.pod5_open_data_example.v1 not_checked ont-open-data eu-west-1
#>                                                                              key
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#>                                                                                              s3_uri
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#>   size_bytes        last_modified_utc
#> 1   47077200 2025-05-19T23:24:03.000Z
#>                                  intended_use
#> 1 Primary routine opt-in example POD5 object.
```

To fetch and inspect the selected POD5 object, opt in explicitly and
cache the file outside the repository:

``` r

if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  fetched <- ont_open_data_fetch(
    key = pod5_example$key,
    cache_dir = tools::R_user_dir("floundeR", which = "cache")
  )

  pod5_info <- pod5_file_info(fetched$cache_path)
  pod5_integrity <- pod5_verify(fetched$cache_path)
}
```

## Represent Reserved Evidence Surfaces

BAM/Bamana and Porkchop/library-preparation evidence are first-class
handoff sections. In this minimal vignette they are intentionally marked
as not checked, because no BAM file or library-preparation screen is
bundled with the package fixture.

``` r

bam_evidence <- data.frame(
  schema_version = "flounder.bam_qc_report_card.v1",
  status = "not_checked",
  check_id = "bam_evidence_not_supplied",
  detail = "No BAM file is bundled with this vignette fixture.",
  stringsAsFactors = FALSE
)

library_preparation_evidence <- data.frame(
  schema_version = "flounder.library_preparation_qc.v1",
  status = "not_checked",
  evidence_kind = "porkchop_reserved",
  detail = paste(
    "Porkchop adapter, primer, barcode, kit, and cDNA evidence can populate",
    "this section when library-preparation screening evidence is supplied."
  ),
  stringsAsFactors = FALSE
)
```

## Create The Synoptikon Payload

[`as_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
assembles the sectioned evidence in memory.
[`write_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
serialises the same contract as JSON for handoff.

``` r

payload <- as_synoptikon_qc(
  run_summary = run_summary,
  report_cards = run_report_card,
  flowcell = list(channel_density = flowcell_density),
  barcode = list(composition = barcode_composition),
  pod5 = list(open_data_example = pod5_handoff),
  bam = list(report_card = bam_evidence),
  library_preparation = list(evidence = library_preparation_evidence),
  input_provenance = list(
    input_id = "example-sequencing-summary",
    kind = "sequencing_summary",
    uri = summary_file,
    size_bytes = unname(file.info(summary_file)$size),
    sha256 = NULL,
    metadata = list(package_fixture = TRUE)
  ),
  run_identity = list(
    run_id = "example-run",
    sample_id = "example-sample",
    project_id = "flounder-revival",
    flow_cell_id = NA_character_
  ),
  governance_context = list(
    tenant_id = "example-tenant",
    governance_domain_id = "example-domain",
    mneion_work_request_id = "example-work-request"
  ),
  limitations = list(
    limitation_id = "example-pod5-not-downloaded",
    severity = "info",
    section = "pod5",
    message = paste(
      "This vignette records ONT open-data POD5 provenance but does not",
      "download POD5 bytes by default."
    )
  ),
  payload_id = "example-synoptikon-qc-payload",
  generated_at_utc = "2026-06-14T06:39:36Z"
)

names(payload)
#>  [1] "schema_version"     "payload_id"         "generated_at_utc"  
#>  [4] "producer"           "run_identity"       "governance_context"
#>  [7] "input_provenance"   "qc_sections"        "report_cards"      
#> [10] "limitations"        "handoff"
names(payload$qc_sections)
#> [1] "sequencing_summary"  "flowcell"            "barcode"            
#> [4] "pod5"                "bam"                 "library_preparation"
payload$handoff
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

``` r

handoff_path <- file.path(tempdir(), "flounder-synoptikon-qc-payload.json")
write_synoptikon_qc(
  handoff_path,
  run_summary = run_summary,
  report_cards = run_report_card,
  flowcell = list(channel_density = flowcell_density),
  barcode = list(composition = barcode_composition),
  pod5 = list(open_data_example = pod5_handoff),
  bam = list(report_card = bam_evidence),
  library_preparation = list(evidence = library_preparation_evidence),
  input_provenance = payload$input_provenance,
  run_identity = payload$run_identity,
  governance_context = payload$governance_context,
  limitations = payload$limitations,
  payload_id = payload$payload_id,
  generated_at_utc = payload$generated_at_utc
)

jsonlite::fromJSON(handoff_path, simplifyVector = FALSE)$schema_version
#> [1] "flounder.synoptikon_qc_payload.v1"
```

The written JSON follows the installed schema at
`inst/schema/synoptikon-qc-payload-v1.schema.json`. Downstream
Synoptikon ingestion can use that schema to validate the handoff before
creating governed workbench records, report-build inputs, or audit
evidence.

## Where Grammateus Fits

The Synoptikon payload is QC evidence. It is not a trusted report
manifest and does not replace Grammateus lifecycle contracts. Once
Grammateus runtime support lands, the same sectioned QC data should feed
semantic report elements, governed figures, provenance tables, and
Mnemosyne Biosciences branded HTML/PDF reports.
