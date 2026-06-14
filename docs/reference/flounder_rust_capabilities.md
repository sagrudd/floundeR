# Query compiled Rust support

`flounder_rust_capabilities()` reports whether the compiled Rust
extension is available and which curated Rust-backed engines have been
linked into the package. `flounder_rust_available()` is a predicate for
package code and downstream tests. `skip_if_no_flounder_rust()` skips a
`testthat` test when a Rust-backed feature is unavailable.

## Usage

``` r
flounder_rust_capabilities(required = TRUE)

flounder_rust_available()

skip_if_no_flounder_rust(feature = "compiled Rust support")
```

## Arguments

- required:

  Logical scalar. If `TRUE`, missing compiled Rust support is reported
  as a typed `floundeR_rust_unavailable` condition. If `FALSE`, a
  capability list with `compiled_support = FALSE` is returned.

- feature:

  Short feature label used in skip messages.

## Value

`flounder_rust_capabilities()` returns a named list with schema version,
package/crate metadata, `compiled_support`, and per-engine link status.
`flounder_rust_available()` returns `TRUE` or `FALSE`.
`skip_if_no_flounder_rust()` returns the capability list invisibly when
Rust is available, otherwise calls
[`testthat::skip()`](https://testthat.r-lib.org/reference/skip.html).
