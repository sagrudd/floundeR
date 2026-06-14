# Grammateus Runtime Interface

This document defines the prebuilt Grammateus runtime interface that public
floundeR builds will use when private Grammateus source is unavailable.
Grammateus remains private and internal-only. floundeR remains open-source and
must install, check, and provide core QC functionality without this runtime.

The interface is intentionally narrow. It exists only to support governed
HTML/PDF QC report rendering, report manifests, provenance, Mnemosyne
Biosciences styling, and Synoptikon trusted-report handoff.

## Runtime Bundle Contract

A Grammateus runtime bundle is a versioned archive distributed outside the
public package source. It may be attached to an authorized GitHub release or a
private Mnemosyne release channel.

Bundle archive name:

```text
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}.sha256
grammateus-runtime-{runtime_version}-{platform}.{tar.gz|zip}.sig
```

Manifest name:

```text
grammateus-runtime-manifest-{runtime_version}.json
grammateus-runtime-manifest-{runtime_version}.json.sig
grammateus-runtime-index-{runtime_version}.json
grammateus-runtime-index-{runtime_version}.json.sig
flounder-runtime-compatibility-{runtime_version}.json
flounder-runtime-compatibility-{runtime_version}.json.sig
```

The detailed GitHub/private-release asset contract, platform matrix, runtime
index shape, compatibility manifest shape, and publication rules are maintained
in `GRAMMATEUS_RELEASE_ASSETS.md`.

The bundle must not contain private Grammateus source code. It may contain:

- compiled runtime library or adapter required by the floundeR Rust extension;
- approved report templates, themes, fonts, and branding assets;
- rendering schemas and runtime contract files;
- license, provenance, and SBOM metadata;
- examples or fixtures needed to validate the installed runtime.

## Manifest Schema

The runtime manifest is the authoritative machine-readable contract. floundeR
helpers must validate it before using the runtime.

Required fields:

```json
{
  "schema_version": "flounder.grammateus_runtime_manifest.v1",
  "runtime_name": "grammateus-runtime",
  "runtime_version": "0.6.0",
  "grammateus_release": "0.6.0",
  "flounder_version_min": "0.16.0",
  "flounder_version_max_exclusive": "0.17.0",
  "platform": "aarch64-apple-darwin",
  "abi": {
    "r_version_min": "4.4.0",
    "rust_toolchain": "1.85",
    "library_kind": "cdylib"
  },
  "artifacts": [
    {
      "path": "lib/libgrammateus_runtime.dylib",
      "kind": "runtime_library",
      "sha256": "sha256:0000000000000000000000000000000000000000000000000000000000000000",
      "byte_len": 0
    },
    {
      "path": "templates/mnemosyne_qc_report_v1.json",
      "kind": "template",
      "sha256": "sha256:0000000000000000000000000000000000000000000000000000000000000000",
      "byte_len": 0
    }
  ],
  "capabilities": {
    "render_figure_html": true,
    "render_figure_pdf": true,
    "render_table_html": true,
    "render_table_pdf": true,
    "render_report_html": true,
    "render_report_pdf": true,
    "trusted_report_manifest": true,
    "mnemosyne_biosciences_theme": true
  },
  "signing": {
    "signature_file": "grammateus-runtime-manifest-0.6.0.json.sig",
    "signature_type": "mnemosyne_governance_envelope_v1",
    "authority_id": "auth_mnemosyne_biosciences",
    "key_id": "authority_key_2026_001"
  },
  "provenance": {
    "built_at_utc": "2026-06-14T00:00:00Z",
    "source_commit": "0fba9ab3f23018133e4b52af47378e5a69519eff",
    "builder": "mnemosyne-release",
    "sbom_path": "sbom/grammateus-runtime.spdx.json"
  }
}
```

Rules:

- `schema_version` must be exactly recognised by floundeR.
- `runtime_name` must be `grammateus-runtime`.
- `runtime_version` must satisfy the floundeR compatibility window.
- `platform` must match the current host or an explicitly selected override.
- Every artifact entry must have `path`, `kind`, `sha256`, and `byte_len`.
- Artifact paths must be relative to the runtime root and must not contain
  `..`, absolute paths, or symlinks escaping the runtime root.
- Runtime libraries, templates, themes, schemas, fonts, and branding assets
  must be listed explicitly.
- Signature verification is required before governed rendering is enabled.

## Discovery Order

floundeR runtime discovery must be deterministic and explicit:

1. `GRAMMATEUS_HOME`
2. `options(floundeR.grammateus_home = "...")`
3. floundeR-managed user cache
4. system installation paths, if later supported

Package load must not download, install, or validate private assets. Runtime
discovery helpers may inspect paths, but installation and network access must
require an explicit user call.

## Install Layout

A valid runtime root should use this layout:

```text
{runtime_root}/
  manifest.json
  manifest.json.sig
  lib/
  templates/
  themes/
  assets/
  schemas/
  licenses/
  sbom/
```

The floundeR-managed cache root should be outside the package repository, for
example:

```text
tools::R_user_dir("floundeR", "cache")/grammateus/{runtime_version}/{platform}
```

Downloaded archives and extracted runtimes must not be committed to the
floundeR repository.

## R Helper Contract

The next implementation slice should expose:

- `grammateus_runtime_available()`: logical scalar, never errors for absence.
- `grammateus_runtime_version()`: character scalar or `NA_character_`.
- `grammateus_runtime_manifest()`: parsed manifest list after basic shape
  validation.
- `grammateus_runtime_validate()`: detailed validation result with typed
  failure categories.
- `grammateus_runtime_install()`: explicit installer for authorized users.

Validation result shape:

```r
list(
  schema_version = "flounder.grammateus_runtime_validation.v1",
  available = TRUE,
  valid = TRUE,
  runtime_root = "/path/to/runtime",
  runtime_version = "0.6.0",
  platform = "aarch64-apple-darwin",
  capabilities = list(render_report_html = TRUE, render_report_pdf = TRUE),
  artifact_count = 12L,
  artifact_sha256 = "sha256:...",
  failures = data.frame(
    category = character(),
    path = character(),
    message = character()
  )
)
```

Typed failure categories:

- `runtime_not_found`
- `manifest_missing`
- `manifest_schema`
- `manifest_signature`
- `artifact_signature`
- `checksum_mismatch`
- `platform_unsupported`
- `version_incompatible`
- `artifact_missing`
- `path_escape`
- `capability_missing`
- `credentials_missing`

## Rust Boundary

Normal report-rendering APIs must call Rust entry points registered in the
floundeR R extension. They must not invoke Grammateus CLI binaries.

The public build may return `runtime_unavailable` from the Rust boundary. An
authorized runtime-enabled build may replace that path with dynamic library or
adapter calls, provided the R-facing response envelope remains stable:

```r
list(
  ok = TRUE,
  data = "<html>...</html>",
  error = NULL,
  category = NULL,
  operation = "figure_html"
)
```

Failure envelopes must include:

```r
list(
  ok = FALSE,
  data = NULL,
  error = "actionable message",
  category = "runtime_unavailable",
  operation = "figure_html"
)
```

## Security And Trust

Before governed rendering is enabled, floundeR must verify:

- runtime root path safety;
- manifest schema;
- manifest signature or governance envelope declaration;
- manifest and artifact signature file presence, path safety, non-empty content,
  and declared SHA-256 checksums where provided;
- every declared artifact checksum and byte length;
- platform and ABI compatibility;
- floundeR/runtime version compatibility;
- required rendering and trusted-report capabilities.

Rendered report manifests must record the Grammateus runtime version, runtime
manifest hash, selected template/theme identifiers, and relevant artifact
hashes. Grammateus must remain rendering infrastructure; Synoptikon/Mneion owns
identity, roles, authority keys, approval policy, and governance signing.

## Public Package Behavior

The open-source package must continue to pass checks without a runtime. In that
state:

- core QC functions work normally;
- plot and figure handoff helpers work normally;
- render calls fail with `flounder_grammateus_runtime_unavailable`;
- tests that require a real private runtime are skipped unless an explicit
  runtime path and opt-in environment variable are supplied;
- examples must not download private runtime assets.
