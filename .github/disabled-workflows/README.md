# Disabled GitHub Actions

GitHub Actions are disabled temporarily while the floundeR revival is producing
heavy code churn.

The workflow definitions are kept here for restoration, but they are not active
while they remain outside `.github/workflows/`.

## Public Check

The public workflow must not require private Grammateus source, private runtime
assets, network-only ONT data, or private credentials. It explicitly clears
`GRAMMATEUS_HOME` and checks that core runtime discovery reports absence
cleanly.

To re-enable the public R package check, move:

```text
.github/disabled-workflows/r-check.yml.disabled
```

back to:

```text
.github/workflows/r-check.yml
```

## Private Grammateus Runtime Check

The private workflow is credentialed and manual-only. It downloads authorized
prebuilt Grammateus runtime assets, verifies the archive checksum, unpacks the
runtime outside the repository, validates the runtime through floundeR, and
runs report-rendering tests with `GRAMMATEUS_HOME` set.

To re-enable it, move:

```text
.github/disabled-workflows/grammateus-runtime-check.yml.disabled
```

to:

```text
.github/workflows/grammateus-runtime-check.yml
```

Required repository secrets:

- `GRAMMATEUS_RUNTIME_ARCHIVE_URL`
- `GRAMMATEUS_RUNTIME_SHA256_URL`
- `GRAMMATEUS_RUNTIME_SIG_URL`
- `GRAMMATEUS_RUNTIME_TOKEN`
