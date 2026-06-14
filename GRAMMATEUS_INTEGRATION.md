# Grammateus Integration Notes

This document records floundeR's Grammateus reporting integration preflight and
the boundaries that keep the public package open-source while Grammateus remains
private reporting infrastructure.

## Checkout Preflight

Checked on 2026-06-14 before adding any Rust dependency or R reporting API:

- local checkout: `../grammateus`
- branch: `main`
- package version declared by Grammateus: `0.6.0`
- local commit: `0fba9ab3f23018133e4b52af47378e5a69519eff`
- remote: `origin` at `https://github.com/sagrudd/grammateus`
- remote comparison after `git fetch --prune`: `0` commits ahead and `0`
  commits behind `origin/main`
- worktree state: clean

This confirms the local private source checkout is current and clean enough for
the next Slice 13 steps: canonical contract verification and public library API
identification.

No Grammateus source code is vendored into floundeR by this preflight. No
floundeR build path, package check, or core QC API depends on private
Grammateus source or runtime assets at this stage.

## Distribution Boundary

The public floundeR package must continue to install, check, and provide core
QC functionality without Grammateus source. Development may use an authorized
private checkout to design and test in-process bindings, but public
distribution must use optional prebuilt Grammateus runtime/reporting artifacts
with explicit discovery, version validation, manifests, and checksums.

The next integration steps must preserve the canonical Grammateus contracts
listed in `GOVERNANCE.md`, especially rendering elements, trusted-report
lifecycle, and trusted-report Synoptikon contracts.
