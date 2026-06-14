# GitHub Installation Paths

This document records the GitHub installation paths for open-source floundeR
users and for authorized users who also receive private prebuilt Grammateus
runtime assets.

floundeR, pod5-tools, bamana, and porkchop are open-source. Grammateus remains
private and optional from the public package perspective. Core nanopore QC,
POD5 discovery and integrity checks, BAM/BGZF/FASTQ QC, library-preparation
evidence, report contracts, and Synoptikon payload construction must remain
usable without private Grammateus assets.

## Public GitHub Source Install

Public users install floundeR from GitHub as an R source package.

System prerequisites:

- R with package development tools.
- Cargo and `rustc` satisfying the floor in `DESCRIPTION`.
- Platform build tools suitable for compiling R packages with native code.

Example using `pak`:

```r
install.packages("pak")
pak::pak("sagrudd/floundeR")
```

Example using `remotes`:

```r
install.packages("remotes")
remotes::install_github("sagrudd/floundeR")
```

This path must not require:

- Grammateus source code;
- Grammateus runtime archives;
- GitHub private-release credentials;
- ONT open-data downloads;
- network-dependent tests or examples during package checks.

After install, public users can confirm the optional Grammateus runtime is
absent without treating absence as a package failure:

```r
library(floundeR)

grammateus_runtime_available()
grammateus_runtime_version()
```

`grammateus_runtime_available()` returning `FALSE` is expected for public core
installs that have not been configured with private runtime assets.

## Authorized Grammateus Runtime Install

Authorized users may receive private prebuilt Grammateus runtime assets through
a GitHub release or a private Mnemosyne release channel. Those assets are not
part of the public floundeR source package.

Expected release asset set:

```text
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}.sha256
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}.sig
grammateus-runtime-manifest-{runtime_version}.json
grammateus-runtime-manifest-{runtime_version}.json.sig
grammateus-runtime-index-{runtime_version}.json
grammateus-runtime-index-{runtime_version}.json.sig
flounder-runtime-compatibility-{runtime_version}.json
flounder-runtime-compatibility-{runtime_version}.json.sig
```

The normative asset naming and manifest contract is maintained in
`GRAMMATEUS_RELEASE_ASSETS.md`.

Authorized users should:

1. Download only the platform-specific runtime assets for their system.
2. Verify the archive checksum before unpacking.
3. Preserve the runtime manifest, compatibility manifest, and signature files
   beside the extracted runtime.
4. Keep all runtime archives and extracted private assets outside the floundeR
   source repository.
5. Configure floundeR with an explicit runtime root.

Example shell outline:

```sh
mkdir -p "${HOME}/.cache/floundeR/grammateus-runtime"
cd "${HOME}/.cache/floundeR/grammateus-runtime"

# Download authorized private assets through the approved GitHub or Mnemosyne
# release channel before running these checks.
shasum -a 256 -c grammateus-runtime-0.6.0-aarch64-apple-darwin.tar.gz.sha256
tar -xzf grammateus-runtime-0.6.0-aarch64-apple-darwin.tar.gz
```

Example R validation:

```r
library(floundeR)

runtime_root <- path.expand("~/.cache/floundeR/grammateus-runtime/0.6.0/aarch64-apple-darwin")
validation <- grammateus_runtime_validate(runtime_root)

stopifnot(isTRUE(validation$valid))
```

Users may configure discovery through `GRAMMATEUS_HOME`:

```sh
export GRAMMATEUS_HOME="${HOME}/.cache/floundeR/grammateus-runtime/0.6.0/aarch64-apple-darwin"
```

or through an R option:

```r
options(floundeR.grammateus_home = runtime_root)
```

When a runtime has already been extracted and validated, it may be copied into
the floundeR-managed cache explicitly:

```r
installed <- grammateus_runtime_install(
  source_runtime_root = runtime_root,
  overwrite = FALSE
)

stopifnot(isTRUE(installed$valid))
```

No helper should silently download private runtime assets during package load,
examples, vignettes, tests, or ordinary QC operations.

## GitHub Actions Posture

GitHub Actions are temporarily disabled during the active reboot. Disabled
workflow templates are stored under `.github/disabled-workflows/`.

The public check template is runtime-free and suitable for ordinary open-source
checks when restored.

The private Grammateus runtime template is manual-only and requires authorized
runtime asset secrets. It downloads private assets, verifies the archive
checksum, validates the runtime through floundeR, and runs report-rendering
tests with `GRAMMATEUS_HOME` set.

Do not move disabled workflows back to `.github/workflows/` until the project
is ready to resume Actions.

## Failure Expectations

Core QC should continue to work when Grammateus is absent. Only governed
HTML/PDF rendering paths should fail, and they should fail with typed,
actionable runtime conditions such as:

- runtime not installed;
- runtime version incompatible;
- manifest, compatibility, checksum, or signature validation failed;
- unsupported platform or ABI;
- private release credentials missing.
