# Real-Data QC Evidence

This document describes how to produce report-grade QC evidence for the
canonical ONT Zymo PAU85136 example without making floundeR a basecaller.

floundeR never provides direct basecalling functionality. The full report
workflow requires basecalled evidence, but that evidence is produced by
Mnematikon or another controlled upstream workflow. floundeR then ingests the
basecalled outputs and creates QC tables, Grammateus plot artifacts, report
contracts, manifests, and provenance.

## Canonical Inputs

The selected raw POD5 source is:

```text
s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
```

The selected fail object remains reserved for fail-state examples:

```text
s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5
```

Do not commit downloaded POD5, FASTQ, BAM, sequencing-summary, or derived
basecalling artifacts to this repository unless they are intentionally tiny test
fixtures with reviewed provenance.

## Workflow Boundary

The ownership boundary is:

- Mnematikon owns basecalling orchestration, model selection, container use,
  accelerator access, and basecalling provenance.
- floundeR owns raw-data metadata checks, basecalled-output ingestion, QC
  tables, plot artifacts, report contracts, report manifests, and Synoptikon
  handoff evidence.
- Grammateus rendering remains optional from the public package perspective.
  floundeR can still emit runtime-free report contracts and plot artifacts.

## Dry Run

The dry run writes source metadata, an analysis status file, a Mnematikon
handoff table, and a workflow manifest into a cache directory. It does not touch
S3 or any DGX host.

```sh
Rscript scripts/build-real-data-qc-evidence.R --dry-run
```

To choose a cache directory:

```sh
FLOUNDER_REAL_DATA_QC_DIR=/path/to/cache \
Rscript scripts/build-real-data-qc-evidence.R --dry-run
```

## Ingest Basecalled Evidence

After Mnematikon has produced a Dorado sequencing summary for the selected
POD5, point floundeR at that file:

```sh
FLOUNDER_REAL_DATA_QC=true \
FLOUNDER_REAL_DATA_QC_SEQUENCE_SUMMARY=/path/to/sequencing_summary.tsv \
Rscript scripts/build-real-data-qc-evidence.R
```

Alternatively, point floundeR at a basecalled output directory and the script
will search for a sequencing-summary file:

```sh
FLOUNDER_REAL_DATA_QC=true \
FLOUNDER_REAL_DATA_QC_BASECALLED_DIR=/path/to/mnematikon/basecalled-output \
Rscript scripts/build-real-data-qc-evidence.R
```

When a sequencing summary is available, the workflow writes:

- `qc-tables/*.tsv`
- `plots/<plot_id>/*.png`
- `plots/<plot_id>/*.svg`
- `report/report_ont_zymo_pau85136_qc-contract.json`
- `report/report_ont_zymo_pau85136_qc-manifest.json`
- `workflow-manifest.json`

The report contract includes the generated figures and tables, but records that
basecalling provenance must come from the upstream Mnematikon artifact manifest.

## DGX And Mnematikon Notes

Hostnames, IP addresses, SSH keys, and container availability are
environment-specific. The handoff table records optional hints through
environment variables such as `FLOUNDER_DGX_HOST`, but those values are not
package contracts and should not be hard-coded into floundeR.

The current development environment may refer to `../mnematikon` and use
available Mnematikon containers for upstream basecalling. That remains an
operational dependency for producing complete example evidence, not a floundeR
runtime feature.

