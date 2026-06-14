# Grammateus Runtime Release Assets

This document defines the GitHub release asset names and manifest set for
prebuilt private Grammateus runtime bundles consumed by open-source floundeR
builds.

The assets are private distribution artifacts. They must not be committed to
the floundeR source repository and must not contain private Grammateus source
code.

## Release Channels

floundeR supports two release channels for Grammateus runtime assets:

- `github-private`: authorized GitHub release assets attached to a floundeR
  release or a private Mnemosyne/Grammateus release repository.
- `mnemosyne-private`: an internal Mnemosyne release channel with the same
  file names, manifests, checksums, and signatures.

The public floundeR package may know the asset naming convention and manifest
schema, but it must not require credentials or private assets during install,
package load, examples, or default tests.

## Platform Identifiers

Use Rust target triples as the platform component:

- `aarch64-apple-darwin`
- `x86_64-apple-darwin`
- `x86_64-unknown-linux-gnu`
- `aarch64-unknown-linux-gnu`
- `x86_64-pc-windows-msvc`

Additional platforms may be added only when the runtime manifest declares the
ABI and the floundeR validator recognises the platform.

## Asset Names

For a Grammateus runtime version `{runtime_version}`, floundeR compatibility
window `{flounder_min}` to `{flounder_max_exclusive}`, and platform
`{platform}`, publish exactly these assets:

```text
grammateus-runtime-{runtime_version}-{platform}.tar.gz
grammateus-runtime-{runtime_version}-{platform}.tar.gz.sha256
grammateus-runtime-{runtime_version}-{platform}.tar.gz.sig
```

Windows builds may use `.zip` instead of `.tar.gz`:

```text
grammateus-runtime-{runtime_version}-{platform}.zip
grammateus-runtime-{runtime_version}-{platform}.zip.sha256
grammateus-runtime-{runtime_version}-{platform}.zip.sig
```

Each release must also publish release-level manifests:

```text
grammateus-runtime-manifest-{runtime_version}.json
grammateus-runtime-manifest-{runtime_version}.json.sig
grammateus-runtime-index-{runtime_version}.json
grammateus-runtime-index-{runtime_version}.json.sig
flounder-runtime-compatibility-{runtime_version}.json
flounder-runtime-compatibility-{runtime_version}.json.sig
```

## Asset Roles

- runtime archive: extracted into the floundeR-managed Grammateus cache or an
  explicitly configured `GRAMMATEUS_HOME`.
- archive `.sha256`: single-line checksum file containing
  `sha256:<hex>  <archive-name>`.
- archive `.sig`: governance signature or signature envelope for the archive.
- runtime manifest: per-runtime contract consumed by
  `grammateus_runtime_validate()`.
- runtime index: release-level index listing every platform archive, archive
  checksum, archive signature, manifest checksum, and expected install root.
- compatibility manifest: floundeR-facing compatibility window and tested
  open-source engine versions for `pod5-tools`, `bamana`, and `porkchop`.
  Installed runtimes should either place this file at
  `flounder-runtime-compatibility-{runtime_version}.json` or declare its
  relative path in `manifest.json` as `compatibility_file`.

## Runtime Index Shape

The runtime index is the public shape floundeR can inspect after an authorized
download. It does not contain private source paths.

```json
{
  "schema_version": "flounder.grammateus_runtime_index.v1",
  "runtime_name": "grammateus-runtime",
  "runtime_version": "0.6.0",
  "release_channel": "github-private",
  "published_at_utc": "2026-06-14T00:00:00Z",
  "manifest_asset": {
    "name": "grammateus-runtime-manifest-0.6.0.json",
    "sha256": "sha256:0000000000000000000000000000000000000000000000000000000000000000",
    "signature_asset": "grammateus-runtime-manifest-0.6.0.json.sig"
  },
  "platforms": [
    {
      "platform": "aarch64-apple-darwin",
      "archive_asset": "grammateus-runtime-0.6.0-aarch64-apple-darwin.tar.gz",
      "archive_sha256": "sha256:0000000000000000000000000000000000000000000000000000000000000000",
      "archive_signature_asset": "grammateus-runtime-0.6.0-aarch64-apple-darwin.tar.gz.sig",
      "archive_bytes": 0,
      "abi": {
        "r_version_min": "4.4.0",
        "rust_toolchain": "1.85",
        "library_kind": "cdylib"
      }
    }
  ],
  "signing": {
    "signature_asset": "grammateus-runtime-index-0.6.0.json.sig",
    "signature_type": "mnemosyne_governance_envelope_v1",
    "authority_id": "auth_mnemosyne_biosciences",
    "key_id": "authority_key_2026_001"
  }
}
```

## Compatibility Manifest Shape

```json
{
  "schema_version": "flounder.grammateus_runtime_compatibility.v1",
  "runtime_version": "0.6.0",
  "flounder_version_min": "0.21.0",
  "flounder_version_max_exclusive": "0.22.0",
  "required_runtime_capabilities": [
    "render_report_html",
    "render_report_pdf",
    "trusted_report_manifest",
    "mnemosyne_biosciences_theme"
  ],
  "tested_open_engines": {
    "pod5_tools": {
      "version": "recorded-at-release",
      "source": "open-source"
    },
    "bamana": {
      "version": "recorded-at-release",
      "source": "open-source"
    },
    "porkchop": {
      "version": "recorded-at-release",
      "source": "open-source"
    }
  },
  "notes": [
    "Core floundeR QC must remain usable without this runtime."
  ]
}
```

## Publication Rules

- Publish all platform archives before publishing the runtime index.
- Publish checksums and signatures for every archive and manifest.
- The runtime manifest must list every file inside the extracted runtime root
  with relative path, kind, byte length, and SHA-256 checksum.
- The runtime index must list every platform archive expected for the release.
- The compatibility manifest must declare the supported floundeR version range.
- Release notes must identify whether the runtime supports HTML, PDF, trusted
  report manifests, Mnemosyne templates/themes, and controlled R plot
  execution.
- Do not publish private Grammateus source, build worktrees, credentials, or
  unredacted internal paths.

## Install Expectations

`grammateus_runtime_install()` should eventually accept either an explicit
runtime root or a downloaded asset set that matches this document. Installation
must remain explicit. The helper must not silently download assets during
package load, examples, vignettes, or ordinary QC operations.
