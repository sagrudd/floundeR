# Package index

## Run-Level QC Contracts

Tidy sequencing-summary QC tables and report-card checks.

- [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md)
  : Summarise a nanopore run for QC handoff
- [`qc_yield_over_time()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  [`qc_read_length_distribution()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  [`qc_quality_distribution()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  [`qc_channel_density()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  [`qc_barcode_composition()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  : Build tidy QC data from sequencing-summary inputs
- [`qc_report_card_thresholds()`](https://sagrudd.github.io/floundeR/reference/qc_report_card.md)
  [`qc_report_card()`](https://sagrudd.github.io/floundeR/reference/qc_report_card.md)
  : Build a run-level QC report card

## POD5 Evidence And ONT Open Data

Rust-backed POD5 discovery, integrity, manifest, and selected ONT Zymo
open-data helpers.

- [`pod5_compare()`](https://sagrudd.github.io/floundeR/reference/pod5_compare.md)
  : Compare two POD5 collections or manifests
- [`pod5_file_info()`](https://sagrudd.github.io/floundeR/reference/pod5_file_info.md)
  : Inspect local POD5 file metadata
- [`pod5_find()`](https://sagrudd.github.io/floundeR/reference/pod5_find.md)
  : Discover local POD5-containing folders
- [`pod5_folder_info()`](https://sagrudd.github.io/floundeR/reference/pod5_folder_info.md)
  : Summarise a local POD5 folder or run tree
- [`pod5_manifest()`](https://sagrudd.github.io/floundeR/reference/pod5_manifest.md)
  : Build a versioned POD5 collection manifest
- [`pod5_subdivide_plan()`](https://sagrudd.github.io/floundeR/reference/pod5_subdivide_plan.md)
  : Plan a read-only POD5 collection subdivision
- [`pod5_verify()`](https://sagrudd.github.io/floundeR/reference/pod5_verify.md)
  : Verify a local POD5 candidate file
- [`ont_open_data_fetch()`](https://sagrudd.github.io/floundeR/reference/ont_open_data_fetch.md)
  : Fetch one ONT open-data object into an explicit cache directory
- [`ont_open_data_list()`](https://sagrudd.github.io/floundeR/reference/ont_open_data_list.md)
  [`flounder_ont_zymo_pod5_prefix()`](https://sagrudd.github.io/floundeR/reference/ont_open_data_list.md)
  : List objects in an ONT open-data S3 prefix
- [`ont_zymo_pod5_dataset()`](https://sagrudd.github.io/floundeR/reference/ont_zymo_pod5_dataset.md)
  : Describe the canonical ONT Zymo fecal POD5 dataset
- [`ont_zymo_pod5_example_objects()`](https://sagrudd.github.io/floundeR/reference/ont_zymo_pod5_example_objects.md)
  : Return selected ONT Zymo fecal POD5 example objects

## BAM/Bamana QC Evidence

Curated alignment-level QC wrappers backed by Bamana semantics.

- [`bam_check_eof()`](https://sagrudd.github.io/floundeR/reference/bam_check_eof.md)
  : Check canonical BGZF EOF evidence
- [`bam_check_index()`](https://sagrudd.github.io/floundeR/reference/bam_check_index.md)
  : Check BAM index evidence for QC report cards
- [`bam_check_map()`](https://sagrudd.github.io/floundeR/reference/bam_check_map.md)
  : Check whether BAM mapping evidence is present
- [`bam_check_sort()`](https://sagrudd.github.io/floundeR/reference/bam_check_sort.md)
  : Check BAM sort declaration and observed sort evidence
- [`bam_check_tag()`](https://sagrudd.github.io/floundeR/reference/bam_check_tag.md)
  : Check BAM aux-tag evidence for library/QC review
- [`bam_qc_report_card()`](https://sagrudd.github.io/floundeR/reference/bam_qc_report_card.md)
  : Build a BAM QC report card
- [`bam_qc_report_card_thresholds()`](https://sagrudd.github.io/floundeR/reference/bam_qc_report_card_thresholds.md)
  : Default BAM report-card thresholds
- [`bam_summary()`](https://sagrudd.github.io/floundeR/reference/bam_summary.md)
  : Summarise a BAM file for nanopore QC
- [`bam_validate()`](https://sagrudd.github.io/floundeR/reference/bam_validate.md)
  : Validate BAM structure and selected consistency checks
- [`bam_verify()`](https://sagrudd.github.io/floundeR/reference/bam_verify.md)
  : Verify a BAM file header and container identity

## Library-Preparation Evidence

Curated Porkchop-backed kit, adapter, primer, barcode, and cDNA QC
evidence.

- [`library_kit_candidates()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  [`library_adapter_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  [`library_barcode_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  [`library_cdna_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  : Porkchop library-preparation QC evidence
- [`library_preparation_report_card()`](https://sagrudd.github.io/floundeR/reference/library_preparation_report_card.md)
  : Build a library-preparation QC report card
- [`library_preparation_report_card_thresholds()`](https://sagrudd.github.io/floundeR/reference/library_preparation_report_card_thresholds.md)
  : Default library-preparation report-card thresholds

## Synoptikon Handoff

Versioned JSON payloads for Mnemosyne/Synoptikon QC review.

- [`as_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
  [`write_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
  : Build and write Synoptikon QC payloads

## Grammateus Reporting

Semantic report elements, figure handoff, controlled R plot rendering,
runtime discovery, and Mnemosyne theme descriptors.

- [`grammateus_figure_from_file()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md)
  [`grammateus_figure_from_ggplot()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md)
  : Build Grammateus figure metadata from image artifacts
- [`grammateus_mnemosyne_theme()`](https://sagrudd.github.io/floundeR/reference/grammateus_mnemosyne_theme.md)
  [`grammateus_apply_theme()`](https://sagrudd.github.io/floundeR/reference/grammateus_mnemosyne_theme.md)
  : Build Mnemosyne Grammateus report theme metadata
- [`grammateus_plot_spec()`](https://sagrudd.github.io/floundeR/reference/grammateus_plot_spec.md)
  [`grammateus_render_plot()`](https://sagrudd.github.io/floundeR/reference/grammateus_plot_spec.md)
  : Build and render Grammateus semantic plot specifications
- [`grammateus_render_element()`](https://sagrudd.github.io/floundeR/reference/grammateus_render_element.md)
  [`grammateus_render_figure_html()`](https://sagrudd.github.io/floundeR/reference/grammateus_render_element.md)
  [`grammateus_render_figure_pdf()`](https://sagrudd.github.io/floundeR/reference/grammateus_render_element.md)
  : Render Grammateus report elements through the Rust binding
- [`grammateus_report_element()`](https://sagrudd.github.io/floundeR/reference/grammateus_report_element.md)
  [`grammateus_qc_report_elements()`](https://sagrudd.github.io/floundeR/reference/grammateus_report_element.md)
  : Build Grammateus semantic report elements
- [`grammateus_runtime_available()`](https://sagrudd.github.io/floundeR/reference/grammateus_runtime_available.md)
  [`grammateus_runtime_version()`](https://sagrudd.github.io/floundeR/reference/grammateus_runtime_available.md)
  [`grammateus_runtime_manifest()`](https://sagrudd.github.io/floundeR/reference/grammateus_runtime_available.md)
  [`grammateus_runtime_validate()`](https://sagrudd.github.io/floundeR/reference/grammateus_runtime_available.md)
  [`grammateus_runtime_install()`](https://sagrudd.github.io/floundeR/reference/grammateus_runtime_available.md)
  : Discover and validate a prebuilt Grammateus runtime
- [`qc_report()`](https://sagrudd.github.io/floundeR/reference/qc_report.md)
  : Assemble a Grammateus-backed QC report contract
- [`qc_plot_yield_over_time()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_quality_distribution()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_read_length_distribution()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_flowcell_density()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_barcode_balance()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_pod5_integrity()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_bam_mapping_summary()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_bam_mapq_distribution()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  [`qc_plot_bam_flag_summary()`](https://sagrudd.github.io/floundeR/reference/qc_plot_yield_over_time.md)
  : Build Grammateus plot specs for nanopore QC report families

## Rust Extension Availability

Helpers for detecting and skipping optional compiled Rust support.

- [`flounder_rust_capabilities()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md)
  [`flounder_rust_available()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md)
  [`skip_if_no_flounder_rust()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md)
  : Query compiled Rust support

## Legacy Sequence And Annotation Helpers

File wrappers and helpers retained outside the retired FAST5 API.

- [`Fasta`](https://sagrudd.github.io/floundeR/reference/Fasta.md) : R6
  Class for loading and analysing FASTA files

- [`Fastq`](https://sagrudd.github.io/floundeR/reference/Fastq.md) : R6
  Class for loading and analysing nanopore (and other) FASTQ files

- [`Blast`](https://sagrudd.github.io/floundeR/reference/Blast.md) :

  R6 Class for loading and analysing BLAST results in basic `Pairwise`
  format

- [`GenbankGenome`](https://sagrudd.github.io/floundeR/reference/GenbankGenome.md)
  : R6 Class for loading and analysing Genbank whole genome files

- [`SequencingSummary`](https://sagrudd.github.io/floundeR/reference/SequencingSummary.md)
  : R6 Class for loading and analysing nanopore sequencing_summary files

- [`SequencingSet`](https://sagrudd.github.io/floundeR/reference/SequencingSet.md)
  : R6 Class for loading and analysing sequence sets

- [`Flowcell`](https://sagrudd.github.io/floundeR/reference/Flowcell.md)
  : R6 Class for performing Flowcell centric analyses

- [`MultiplexSet`](https://sagrudd.github.io/floundeR/reference/MultiplexSet.md)
  : R6 Class for loading, visualising and analysing barcode information

- [`TemporalSet`](https://sagrudd.github.io/floundeR/reference/TemporalSet.md)
  : R6 Class for analysing sequence sets with accompanying temporal data

- [`BamFile`](https://sagrudd.github.io/floundeR/reference/BamFile.md) :
  R6 class for legacy BAM file references

- [`FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
  : R6 Class for floundeR based analyses

- [`flnDr()`](https://sagrudd.github.io/floundeR/reference/flnDr.md) :
  extract a floundeR packaged file

- [`phredmean()`](https://sagrudd.github.io/floundeR/reference/phredmean.md)
  : calculate mean Phred scores from list of Q values

- [`qualToMeanQ()`](https://sagrudd.github.io/floundeR/reference/qualToMeanQ.md)
  : calculate mean Phred score from an ASCII encoded phred string

- [`to_flowcell()`](https://sagrudd.github.io/floundeR/reference/to_flowcell.md)
  :

  Prepare a `flowcell` object from provided class

- [`to_sequencing_set()`](https://sagrudd.github.io/floundeR/reference/to_sequencing_set.md)
  :

  Prepare a `SequencingSet` object from provided class

## Visualisation Helpers

Legacy graphical classes and themes retained for compatibility.

- [`Angenieux`](https://sagrudd.github.io/floundeR/reference/Angenieux.md)
  : R6 Class for visualising floundeR based datasets
- [`AngenieuxDecoration`](https://sagrudd.github.io/floundeR/reference/AngenieuxDecoration.md)
  : R6 Class for describing additional Angenieux decorations
- [`Infographic`](https://sagrudd.github.io/floundeR/reference/Infographic.md)
  : R6 Class for loading and analysing sequence sets
- [`InfographicItem`](https://sagrudd.github.io/floundeR/reference/InfographicItem.md)
  : R6 Class for loading and analysing sequence sets
