library(floundeR)

local_fake_grammateus_runtime <- function(
    capabilities = list(
      render_report_html = TRUE,
      render_report_pdf = TRUE,
      render_element_html = TRUE,
      render_element_pdf = TRUE
    ),
    corrupt_artifact_hash = FALSE,
    include_artifact_signatures = FALSE,
    corrupt_signature_hash = FALSE,
    empty_manifest_signature = FALSE) {
  runtime <- tempfile("grammateus-runtime-")
  dir.create(file.path(runtime, "lib"), recursive = TRUE)
  dir.create(file.path(runtime, "templates"), recursive = TRUE)
  artifacts <- c(
    "lib/libgrammateus_runtime.dylib",
    "templates/mnemosyne-qc-theme.json"
  )
  writeLines("runtime-binary-placeholder", file.path(runtime, artifacts[[1L]]))
  writeLines("{}", file.path(runtime, artifacts[[2L]]))
  artifact_manifest <- lapply(artifacts, function(path) {
    full_path <- file.path(runtime, path)
    hash <- paste0("sha256:", unname(tools::sha256sum(full_path)))
    if (isTRUE(corrupt_artifact_hash) && identical(path, artifacts[[1L]])) {
      hash <- paste0("sha256:", strrep("0", 64L))
    }
    list(
      path = path,
      sha256 = hash,
      byte_len = unname(file.info(full_path)$size)
    )
  })
  if (isTRUE(include_artifact_signatures)) {
    artifact_manifest <- Map(function(artifact, path) {
      signature_file <- paste0(path, ".sig")
      writeLines(paste("fixture-signature-for", path),
                 file.path(runtime, signature_file))
      artifact$signature_file <- signature_file
      artifact$signature_sha256 <- paste0(
        "sha256:",
        unname(tools::sha256sum(file.path(runtime, signature_file)))
      )
      artifact
    }, artifact_manifest, artifacts)
  }
  manifest <- list(
    schema_version = "flounder.grammateus_runtime_manifest.v1",
    runtime_name = "grammateus-runtime",
    runtime_version = "0.6.0",
    grammateus_release = "grammateus-0.6.0",
    flounder_version_min = "0.0.0",
    flounder_version_max_exclusive = "99.0.0",
    platform = R.version$platform,
    abi = list(r = as.character(getRversion()), rust = "1.85"),
    artifacts = artifact_manifest,
    capabilities = capabilities,
    signing = list(
      signature_file = "manifest.json.sig",
      signature_algorithm = "fixture"
    ),
    provenance = list(
      source = "test fixture",
      built_at = "2026-06-14T00:00:00Z"
    )
  )
  jsonlite::write_json(
    manifest,
    path = file.path(runtime, "manifest.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  if (isTRUE(empty_manifest_signature)) {
    file.create(file.path(runtime, "manifest.json.sig"))
  } else {
    writeLines("fixture-signature", file.path(runtime, "manifest.json.sig"))
  }
  manifest$signing$signature_sha256 <- paste0(
    "sha256:",
    unname(tools::sha256sum(file.path(runtime, "manifest.json.sig")))
  )
  if (isTRUE(corrupt_signature_hash)) {
    manifest$signing$signature_sha256 <- paste0("sha256:", strrep("1", 64L))
  }
  jsonlite::write_json(
    manifest,
    path = file.path(runtime, "manifest.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  runtime
}

test_that("Grammateus runtime validation reports absent runtimes cleanly", {
  old_env <- Sys.getenv("GRAMMATEUS_HOME", unset = NA_character_)
  old_option <- getOption("floundeR.grammateus_home", default = NULL)
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("GRAMMATEUS_HOME")
    } else {
      Sys.setenv(GRAMMATEUS_HOME = old_env)
    }
    options(floundeR.grammateus_home = old_option)
  }, add = TRUE)
  Sys.unsetenv("GRAMMATEUS_HOME")
  options(floundeR.grammateus_home = NULL)

  validation <- grammateus_runtime_validate(runtime_root = tempfile("missing-"))

  expect_s3_class(validation, "flounder_grammateus_runtime_validation")
  expect_false(validation$available)
  expect_false(validation$valid)
  expect_equal(validation$failures$category, "runtime_not_found")
  expect_false(grammateus_runtime_available(runtime_root = tempfile("missing-")))
  expect_identical(grammateus_runtime_version(runtime_root = tempfile("missing-")),
                   NA_character_)
})

test_that("Grammateus runtime validation accepts a complete local runtime", {
  runtime <- local_fake_grammateus_runtime(include_artifact_signatures = TRUE)

  validation <- grammateus_runtime_validate(runtime)

  expect_true(validation$available)
  expect_true(validation$valid)
  expect_equal(validation$runtime_version, "0.6.0")
  expect_equal(validation$platform, R.version$platform)
  expect_equal(validation$artifact_count, 2L)
  expect_match(validation$artifact_sha256, "^sha256:[0-9a-f]{64}$")
  expect_equal(nrow(validation$failures), 0L)
  expect_true(grammateus_runtime_available(runtime))
  expect_equal(grammateus_runtime_version(runtime), "0.6.0")
  expect_equal(grammateus_runtime_manifest(runtime)$runtime_name,
               "grammateus-runtime")
})

test_that("Grammateus runtime validation reports manifest and artifact failures", {
  runtime <- local_fake_grammateus_runtime(
    capabilities = list(render_report_html = TRUE),
    corrupt_artifact_hash = TRUE
  )

  validation <- grammateus_runtime_validate(runtime)

  expect_true(validation$available)
  expect_false(validation$valid)
  expect_true("capability_missing" %in% validation$failures$category)
  expect_true("checksum_mismatch" %in% validation$failures$category)
})

test_that("Grammateus runtime validation verifies signature artifacts", {
  runtime <- local_fake_grammateus_runtime(
    include_artifact_signatures = TRUE,
    corrupt_signature_hash = TRUE
  )

  validation <- grammateus_runtime_validate(runtime)

  expect_true(validation$available)
  expect_false(validation$valid)
  expect_true("checksum_mismatch" %in% validation$failures$category)

  runtime <- local_fake_grammateus_runtime(
    include_artifact_signatures = TRUE,
    empty_manifest_signature = TRUE
  )
  validation <- grammateus_runtime_validate(runtime)
  expect_false(validation$valid)
  expect_true("manifest_signature" %in% validation$failures$category)

  runtime <- local_fake_grammateus_runtime(include_artifact_signatures = TRUE)
  unlink(file.path(runtime, "lib/libgrammateus_runtime.dylib.sig"))
  validation <- grammateus_runtime_validate(runtime)
  expect_false(validation$valid)
  expect_true("artifact_signature" %in% validation$failures$category)
})

test_that("Grammateus runtime install copies a validated runtime into cache", {
  source <- local_fake_grammateus_runtime()
  cache_root <- tempfile("grammateus-cache-")

  installed <- grammateus_runtime_install(source, cache_root = cache_root)

  expect_true(installed$valid)
  expect_true(dir.exists(installed$runtime_root))
  expect_match(installed$runtime_root,
               "grammateus/0[.]6[.]0/.+$")
  expect_error(
    grammateus_runtime_install(source, cache_root = cache_root),
    class = "flounder_grammateus_runtime_install_error"
  )
})
