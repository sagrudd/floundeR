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
- Record the source bucket, source key, source object size, source timestamp,
  `pod5-tools` version, subdivision strategy, target values, selected relative
  paths, output file size, and output checksum.
- Store large derived files in an explicit cache or artifact location, not in
  the package repository.
- Commit only tiny derived fixtures if they are intentionally created for unit
  tests and are small enough for normal package distribution.

## Access From R

`aws.s3` works for anonymous listing when the bucket region is supplied:

```r
library(aws.s3)

Sys.setenv(AWS_NO_SIGN_REQUEST = "true")

objects <- aws.s3::get_bucket(
  bucket = "ont-open-data",
  prefix = "zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass",
  max = 5,
  region = "eu-west-1"
)
```

Future package helpers should wrap this pattern in a small, explicit API such
as `ont_open_data_list()` and `ont_open_data_fetch()`.

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
