# R6 class for legacy BAM file references

`BamFile` is a minimal legacy holder for a BAM file path. The rebooted
BAM QC surface is expected to move to curated in-process Rust bindings
through Bamana, but this class remains exported until that replacement
API is defined.

## Format

An `R6Class` generator object.

## Methods

### Method `new()`

#### Usage

    BamFile$new(bamfile)

#### Arguments

- `bamfile`:

  Path to a BAM file.

#### Returns

A new `BamFile` object.
