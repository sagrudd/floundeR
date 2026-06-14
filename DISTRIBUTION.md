# Distribution Plan

`floundeR` will remain fully open-source. `../pod5-tools`, `../bamana`, and
`../porkchop` will also remain open and available. Grammateus will remain
private, so floundeR must not require Grammateus source code to install, check,
or use the core QC toolbox.

## Package Shape

The public floundeR package should have three layers:

1. Open-source R package core.
2. Open-source Rust-backed QC engines from `pod5-tools`, `bamana`, and
   `porkchop`.
3. Optional private Grammateus runtime for governed HTML/PDF report rendering.

Core QC, data import, POD5 metadata checks, BAM QC summaries, adapter/primer
and kit evidence, report-card data, and synoptikon payload generation should
work without private Grammateus source. When Grammateus is unavailable, report
APIs should fail with a clear typed R condition that explains how to install or
configure the runtime.

## Grammateus Distribution

Do not vendor Grammateus source into floundeR.

Initial distribution should use GitHub release assets attached to floundeR or a
private Mnemosyne/Grammateus release channel. The assets should be prebuilt for
the supported platforms and referenced by an explicit manifest.

Recommended artifact set per release:

- `grammateus-runtime-{version}-{platform}.{tar.gz|zip}`
- `grammateus-runtime-{version}-{platform}.{tar.gz|zip}.sha256`
- `grammateus-runtime-{version}-{platform}.{tar.gz|zip}.sig`
- `grammateus-runtime-manifest-{version}.json`
- `grammateus-runtime-manifest-{version}.json.sig`
- `grammateus-runtime-index-{version}.json`
- `grammateus-runtime-index-{version}.json.sig`
- `flounder-runtime-compatibility-{version}.json`
- `flounder-runtime-compatibility-{version}.json.sig`
- release notes identifying the compatible floundeR version range

The normative GitHub/private-release asset names, platform identifiers,
release-level runtime index, compatibility manifest, and publication rules are
defined in `GRAMMATEUS_RELEASE_ASSETS.md`.

The runtime archive should contain only redistributable runtime content needed
by floundeR reporting:

- dynamic/static library or executable adapter required by the R extension;
- report templates and themes;
- Mnemosyne Biosciences branding assets approved for report output;
- schema/contracts needed at runtime;
- license and provenance metadata;
- SBOM or dependency manifest where available.

The archive must not contain private source code.

## Runtime Discovery

floundeR should discover Grammateus in this order:

1. explicit `GRAMMATEUS_HOME`;
2. package option such as `options(floundeR.grammateus_home = "...")`;
3. user cache managed by a helper such as `grammateus_runtime_install()`;
4. system installation paths, if later supported.

Helpers should include:

- `grammateus_runtime_available()`;
- `grammateus_runtime_version()`;
- `grammateus_runtime_install(version = NULL, platform = NULL)`;
- `grammateus_runtime_validate()`;
- `grammateus_runtime_manifest()`.

Installation must be explicit. Do not silently download private runtime assets
during package load, examples, tests, or ordinary QC operations.

The normative runtime-bundle, manifest, validation, discovery, and Rust response
contracts are defined in `GRAMMATEUS_RUNTIME_INTERFACE.md`. In short, public
floundeR builds consume a signed manifest plus a platform-specific runtime
archive rather than private Grammateus source. The manifest must declare:

- `schema_version = "flounder.grammateus_runtime_manifest.v1"`;
- the Grammateus runtime and release versions;
- the compatible floundeR version range;
- the target platform and ABI;
- every runtime library, template, theme, schema, branding, license, and SBOM
  artifact with relative path, byte length, and SHA-256 checksum;
- the rendering and trusted-report capabilities available in the runtime;
- signature/governance-envelope metadata;
- build provenance including source commit and SBOM path.

Runtime validation must fail closed for missing manifests, invalid signatures,
checksum mismatches, platform/ABI mismatches, incompatible versions, unsafe
paths, missing artifacts, and missing capabilities. Runtime-aware report
rendering must still use registered Rust entry points inside the R extension;
it must not call Grammateus CLI binaries as the normal package API.

## GitHub Distribution

Initial GitHub distribution should provide:

- source package install from GitHub;
- binary package builds where practical;
- private Grammateus runtime assets for authorized users;
- a public floundeR release manifest that records which open-source
  `pod5-tools`, `bamana`, and `porkchop` versions were tested;
- a private Grammateus runtime manifest that records compatible runtime builds.

The concrete public GitHub source-install path and authorized private
Grammateus runtime setup path are documented in `GITHUB_INSTALLATION.md`.

GitHub Actions should build and check floundeR without private Grammateus
access. A separate private or credentialed workflow may verify Grammateus report
rendering using private runtime assets.

During the active reboot, workflow files remain disabled under
`.github/disabled-workflows/`. The disabled public workflow is explicitly
runtime-free and suitable for open-source checks once restored. The disabled
private workflow is manual-only, credentialed, downloads authorized prebuilt
Grammateus runtime assets, verifies the archive checksum, validates the runtime
through floundeR, and then runs the Grammateus report-rendering tests with
`GRAMMATEUS_HOME` set.

## Future Bioconductor Posture

If floundeR is submitted to Bioconductor in future, the package must remain
installable, checkable, and useful without access to private Grammateus source
or private downloads.

Bioconductor-facing posture:

- keep Grammateus support optional;
- do not require private runtime downloads during install or check;
- skip Grammateus integration tests unless an explicit environment variable and
  runtime path are supplied;
- keep examples runnable without Grammateus or mark them as non-running where
  policy allows;
- document that governed PDF/HTML rendering requires an optional private
  Grammateus runtime distributed outside the open-source package.

Before any Bioconductor submission, re-check current Bioconductor policy around
external binaries, downloads, optional system requirements, and non-open
components.

## Open Rust Dependencies

`pod5-tools`, `bamana`, and `porkchop` should stay open-source and can be
integrated more directly than Grammateus.

The preferred model is to expose only curated QC/review library APIs through the
floundeR Rust extension. Do not expose their entire command surfaces just
because the libraries are available.

## Security And Verification

Private Grammateus runtime artifacts should be verified before use:

- checksum verification;
- manifest signature file verification;
- artifact signature file verification when declared by the runtime manifest;
- version compatibility checks against floundeR in both the runtime manifest
  and release-level compatibility manifest;
- platform and ABI checks;
- clear provenance in rendered reports.

Report manifests should record:

- floundeR version;
- Grammateus runtime version and artifact hash;
- `pod5-tools` and `bamana` versions where used;
- `porkchop` version where adapter, primer, barcode, kit, or cDNA evidence is
  used;
- R version and relevant R package versions;
- input data source hashes and external object keys;
- generated plot/report artifact hashes.

## Failure Modes

Core QC APIs must not fail merely because Grammateus is absent.

Only report-rendering calls that require Grammateus should fail, and they should
fail with actionable messages:

- runtime not installed;
- runtime version incompatible;
- manifest or signature verification failed;
- platform unsupported;
- private artifact credentials missing;
- governed container/R backend unavailable.

This boundary is enforced by package tests that deliberately point Grammateus
discovery at missing runtime paths while exercising core sequencing-summary QC,
tidy QC tables, report cards, POD5 discovery, and Synoptikon JSON export.
Release-style checks must continue to pass without any private Grammateus
runtime assets present.
