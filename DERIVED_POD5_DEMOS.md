# Derived POD5 Demonstration Artifacts

This document defines how floundeR should use lightweight POD5 demonstration
artifacts derived from the canonical ONT Zymo fecal POD5 source.

The goal is to make examples, reports, and Synoptikon demos practical without
committing downloaded or generated POD5 files to this repository.

## Source Objects

Derived demonstration artifacts should start from the selected ONT objects in
`DATASETS.md`:

- routine pass source:
  `s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5`
- fail-state source:
  `s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5`

Use the pass source by default. Use the fail source only when demonstrating
fail-state QC evidence.

## Maintained Workflow

Use `scripts/derive-pod5-demo-workflow.R` as the maintained preparation entry
point. The workflow is opt-in and writes metadata outside the repository by
default.

Metadata-only dry run:

```sh
Rscript scripts/derive-pod5-demo-workflow.R --dry-run
```

Local-source planning run:

```sh
FLOUNDER_DERIVE_POD5_DEMO=true \
FLOUNDER_DERIVED_POD5_SOURCE=/path/to/PAU85136_pass_279c9095_68316534_8289.pod5 \
Rscript scripts/derive-pod5-demo-workflow.R
```

Explicit opt-in download plus planning run:

```sh
FLOUNDER_DERIVE_POD5_DEMO=true \
FLOUNDER_RUN_NETWORK_TESTS=true \
FLOUNDER_DERIVE_POD5_DEMO_DOWNLOAD=true \
Rscript scripts/derive-pod5-demo-workflow.R
```

The default output directory is
`tools::R_user_dir("floundeR", "cache")/derived-pod5-demo`. Override it with
`FLOUNDER_DERIVED_POD5_DEMO_DIR`.

## Produced Metadata

The maintained workflow writes:

- `source-metadata.tsv`: selected ONT source object metadata;
- `demo-workflow.json`: floundeR version, `pod5-tools` source, source object,
  cache path, and workflow policy;
- `README.md`: human-readable cache-directory summary;
- `subdivide-plan.tsv`: read-only `pod5_subdivide_plan()` output when a local
  POD5 source is supplied and Rust support is available;
- `pod5-manifest.tsv`: `pod5_manifest()` output for the local source when
  available.

By design, the workflow does not write derived POD5 files. If a future
controlled `pod5-tools` library path creates a small derived artifact, record at
least:

- source bucket, key, size, timestamp, and checksum or manifest reference;
- source cache path and whether it was downloaded by this workflow;
- `pod5-tools` version or path dependency source;
- subdivision strategy and parameters;
- selected source read/file identifiers or relative paths;
- output artifact path, size, checksum, and read count;
- floundeR version and generated timestamp.

## Repository Policy

Do not commit:

- downloaded ONT POD5 files;
- large derived POD5 files;
- local cache directories;
- transient workflow outputs from a developer machine.

Commit only:

- documentation describing the workflow;
- small text manifests or fixtures when they are intentionally curated for
  package tests;
- tiny POD5 fixtures only when they are deliberately added under the repository
  fixture policy and reviewed with visible provenance.

## Reporting And Synoptikon Use

Reports and Synoptikon handoff payloads should reference derived demonstrations
through their metadata and checksums. They should not assume the large source or
derived POD5 bytes are present inside the package.

For routine examples, cite the pass source and the derived manifest. For
failure-state examples, cite the fail source and make the fail-state purpose
explicit in the provenance.
