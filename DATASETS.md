# External Datasets

`floundeR` should use real nanopore data for examples and integration tests,
but large raw-data files must not be committed to the package repository.

## ONT Zymo Fecal POD5 Dataset

Canonical contemporary POD5 example dataset:

```text
s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/
```

Verified on 2026-06-13 using the public S3 listing API and `aws.s3` anonymous
bucket access.

Bucket details:

- bucket: `ont-open-data`
- region: `eu-west-1`
- prefix: `zymo_fecal_2025.05/raw/PAU85136/pod5/`
- total POD5 objects observed: `9107`
- pass POD5 objects observed: `8290`
- fail POD5 objects observed: `817`
- total listed bytes: `2676973535744`

Selected example objects:

```text
zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5
```

Decision:

- Use `PAU85136_pass_279c9095_68316534_8289.pod5` as the primary routine
  example file. It is the smallest pass POD5 observed in the prefix at
  `47077200` bytes, last modified `2025-05-19T23:24:03.000Z`.
- Use `PAU85136_fail_279c9095_68316534_0.pod5` only when an example needs a
  fail-state POD5 file. It is the smallest fail POD5 observed in the prefix at
  `163007608` bytes, last modified `2025-05-19T14:31:37.000Z`.

The combined selected example set is about 210 MB. That is acceptable for
explicit opt-in integration examples, but too large for package fixtures or
default tests. Many other objects are hundreds of megabytes, so examples should
download only these named objects unless a user explicitly requests a broader
sample.

## Derived POD5 Demonstration Subsets

The selected ONT files are source inputs, not necessarily the final files used
in every demonstration. `../pod5-tools` may be used to create smaller POD5
subsets or split outputs for documentation, tests, report examples, and
synoptikon demos.

Rules for derived subsets:

- Derive subsets from the selected pass/fail source files by default.
- Keep generated demonstration POD5 files as small as practical.
- Plan candidate splits first with `pod5_subdivide_plan()` and record the
  resulting plan table before creating any derived artifact.
- Use `scripts/derive-pod5-demo-workflow.R` as the maintained opt-in workflow
  for example preparation. By default it writes source metadata, a workflow
  manifest, and, when a local source POD5 is supplied, read-only plan/manifest
  tables outside the repository. It does not write derived POD5 files.
- Record the source bucket, source key, source object size, source timestamp,
  `pod5-tools` version, subdivision strategy, target values, selected relative
  paths, output file size, and output checksum.
- Store large derived files in an explicit cache or artifact location, not in
  the package repository.
- Commit only tiny derived fixtures if they are intentionally created for unit
  tests and are small enough for normal package distribution.

Example metadata-only dry run:

```sh
Rscript scripts/derive-pod5-demo-workflow.R --dry-run
```

Example local-source planning run:

```sh
FLOUNDER_DERIVE_POD5_DEMO=true \
FLOUNDER_DERIVED_POD5_SOURCE=/path/to/PAU85136_pass_279c9095_68316534_8289.pod5 \
Rscript scripts/derive-pod5-demo-workflow.R
```

Example explicit source download and planning run:

```sh
FLOUNDER_DERIVE_POD5_DEMO=true \
FLOUNDER_RUN_NETWORK_TESTS=true \
FLOUNDER_DERIVE_POD5_DEMO_DOWNLOAD=true \
Rscript scripts/derive-pod5-demo-workflow.R
```

## Access From R

`floundeR` wraps the canonical ONT Zymo POD5 source in small helpers so S3
bucket details are not repeated across examples:

```r
ont_zymo_pod5_dataset()
ont_zymo_pod5_example_objects(role = "pass")
ont_zymo_pod5_example_objects(role = "fail")
```

The pass object is the default real-data demonstration input. The fail object is
only for fail-state examples. Both remain opt-in downloads and must be cached
outside the repository.

Anonymous S3 listing is still available for controlled inspection:

```r
objects <- ont_open_data_list(
  prefix = paste0(
    flounder_ont_zymo_pod5_prefix(),
    "PAU85136_pass"
  ),
  max = 5,
  region = "eu-west-1",
  anonymous = TRUE
)
```

Downloads must select exactly one object key:

```r
if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  object <- ont_zymo_pod5_example_objects(role = "pass")
  ont_open_data_fetch(
    key = object$key,
    cache_dir = tools::R_user_dir("floundeR", "cache")
  )
}
```

## Repository Policy

- Do not commit downloaded POD5 files.
- The repository `.gitignore` ignores `*.pod5`, common local ONT cache
  directories, and temporary POD5 download names by default. Intentional tiny
  POD5 fixtures must use the `tests/testthat/fixtures/tiny-*.pod5` naming
  convention and be added deliberately after reviewer-visible provenance is
  recorded.
- Do not make package tests depend on network access.
- Use tiny generated fixtures for unit tests.
- Use the selected pass file for routine opt-in integration tests, vignettes,
  benchmarking, and report examples.
- Use the selected fail file only for failure-state QC/report examples.
- Prefer small `pod5-tools`-derived subsets/splits when full selected objects
  are larger than the workflow requires.
- Cache downloaded files outside the repository by default.
- Record object keys, sizes, timestamps, and checksums/manifests in examples
  when they are used.
