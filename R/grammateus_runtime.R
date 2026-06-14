#' Discover and validate a prebuilt Grammateus runtime
#'
#' These helpers inspect optional private Grammateus runtime bundles without
#' making Grammateus a required dependency of the public floundeR package.
#' Discovery checks `GRAMMATEUS_HOME`, then
#' `options(floundeR.grammateus_home = "...")`, then the floundeR user cache.
#' No helper downloads private assets during package load or ordinary QC work.
#'
#' @param runtime_root Optional explicit runtime root. When omitted, the
#'   standard discovery order is used.
#' @param required_capabilities Character vector of manifest capabilities that
#'   must be `TRUE`.
#' @param source_runtime_root Existing extracted runtime root to copy into the
#'   floundeR-managed cache.
#' @param cache_root Optional cache root. Defaults to
#'   `tools::R_user_dir("floundeR", "cache")`.
#' @param overwrite Whether `grammateus_runtime_install()` may replace an
#'   existing cached runtime directory.
#'
#' @return
#' `grammateus_runtime_available()` returns a logical scalar.
#' `grammateus_runtime_version()` returns a character scalar or `NA_character_`.
#' `grammateus_runtime_manifest()` returns the parsed manifest list or `NULL`.
#' `grammateus_runtime_validate()` returns a validation list with schema version,
#' availability, validity, runtime path, version, platform, capabilities,
#' artifact counts, aggregate artifact hash, and failure records.
#' `grammateus_runtime_install()` returns the validation list for the installed
#' runtime.
#'
#' @export
grammateus_runtime_available <- function(runtime_root = NULL) {
  isTRUE(grammateus_runtime_validate(runtime_root = runtime_root)$valid)
}

#' @rdname grammateus_runtime_available
#' @export
grammateus_runtime_version <- function(runtime_root = NULL) {
  manifest <- grammateus_runtime_manifest(runtime_root = runtime_root)
  version <- manifest$runtime_version
  if (!is.character(version) || length(version) != 1L || is.na(version)) {
    return(NA_character_)
  }
  version
}

#' @rdname grammateus_runtime_available
#' @export
grammateus_runtime_manifest <- function(runtime_root = NULL) {
  runtime_root <- .grammateus_runtime_root(runtime_root)
  if (is.na(runtime_root)) {
    return(NULL)
  }
  manifest_path <- file.path(runtime_root, "manifest.json")
  if (!file.exists(manifest_path)) {
    return(NULL)
  }
  tryCatch(
    jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
    error = function(error) NULL
  )
}

#' @rdname grammateus_runtime_available
#' @export
grammateus_runtime_validate <- function(
    runtime_root = NULL,
    required_capabilities = c("render_report_html", "render_report_pdf")) {
  runtime_root <- .grammateus_runtime_root(runtime_root)
  failures <- .grammateus_runtime_failures()
  manifest <- NULL
  capabilities <- list()
  runtime_version <- NA_character_
  platform <- NA_character_
  artifact_hashes <- character()
  artifact_count <- 0L

  if (is.na(runtime_root)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "runtime_not_found",
      NA_character_,
      "No Grammateus runtime root was discovered."
    )
    return(.grammateus_runtime_validation(
      available = FALSE,
      valid = FALSE,
      runtime_root = NA_character_,
      runtime_version = runtime_version,
      platform = platform,
      capabilities = capabilities,
      artifact_count = artifact_count,
      artifact_sha256 = NA_character_,
      failures = failures
    ))
  }

  manifest_path <- file.path(runtime_root, "manifest.json")
  if (!file.exists(manifest_path)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "manifest_missing",
      manifest_path,
      "Runtime manifest.json is missing."
    )
  } else {
    manifest <- tryCatch(
      jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
      error = function(error) {
        failures <<- .grammateus_runtime_add_failure(
          failures,
          "manifest_schema",
          manifest_path,
          paste("Runtime manifest is not valid JSON:", conditionMessage(error))
        )
        NULL
      }
    )
  }

  if (is.list(manifest)) {
    failures <- .grammateus_runtime_validate_manifest_shape(
      manifest,
      runtime_root,
      failures
    )
    runtime_version <- manifest$runtime_version %||% NA_character_
    platform <- manifest$platform %||% NA_character_
    capabilities <- manifest$capabilities %||% list()
    failures <- .grammateus_runtime_validate_compatibility(
      manifest,
      failures
    )
    failures <- .grammateus_runtime_validate_compatibility_manifest(
      manifest,
      runtime_root,
      failures
    )
    failures <- .grammateus_runtime_validate_capabilities(
      capabilities,
      required_capabilities,
      failures
    )
    artifacts <- .grammateus_runtime_validate_artifacts(
      manifest$artifacts,
      runtime_root,
      failures
    )
    failures <- artifacts$failures
    artifact_count <- artifacts$count
    artifact_hashes <- artifacts$hashes
  }

  .grammateus_runtime_validation(
    available = TRUE,
    valid = nrow(failures) == 0L,
    runtime_root = runtime_root,
    runtime_version = runtime_version,
    platform = platform,
    capabilities = capabilities,
    artifact_count = artifact_count,
    artifact_sha256 = .grammateus_runtime_aggregate_hash(artifact_hashes),
    failures = failures
  )
}

#' @rdname grammateus_runtime_available
#' @export
grammateus_runtime_install <- function(
    source_runtime_root,
    cache_root = tools::R_user_dir("floundeR", "cache"),
    overwrite = FALSE) {
  source_runtime_root <- .grammateus_required_runtime_root(
    source_runtime_root,
    "source_runtime_root"
  )
  validation <- grammateus_runtime_validate(runtime_root = source_runtime_root)
  if (!isTRUE(validation$valid)) {
    .grammateus_runtime_abort(
      "Cannot install an invalid Grammateus runtime.",
      "flounder_grammateus_runtime_validation_error",
      validation = validation
    )
  }
  cache_root <- .grammateus_required_runtime_root(
    cache_root,
    "cache_root",
    must_exist = FALSE
  )
  platform <- validation$platform
  runtime_version <- validation$runtime_version
  destination <- file.path(cache_root, "grammateus", runtime_version, platform)
  if (dir.exists(destination)) {
    if (!isTRUE(overwrite)) {
      .grammateus_runtime_abort(
        paste("Grammateus runtime cache already exists:", destination),
        "flounder_grammateus_runtime_install_error"
      )
    }
    unlink(destination, recursive = TRUE, force = TRUE)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(
    from = source_runtime_root,
    to = dirname(destination),
    recursive = TRUE,
    copy.mode = FALSE,
    copy.date = TRUE
  )
  if (!isTRUE(ok)) {
    .grammateus_runtime_abort(
      paste("Failed to copy Grammateus runtime into cache:", destination),
      "flounder_grammateus_runtime_install_error"
    )
  }
  if (!identical(basename(source_runtime_root), basename(destination))) {
    copied <- file.path(dirname(destination), basename(source_runtime_root))
    if (dir.exists(copied)) {
      renamed <- file.rename(copied, destination)
      if (!isTRUE(renamed)) {
        .grammateus_runtime_abort(
          paste("Failed to move copied Grammateus runtime into cache:", destination),
          "flounder_grammateus_runtime_install_error"
        )
      }
    }
  }
  grammateus_runtime_validate(runtime_root = destination)
}

.grammateus_runtime_root <- function(runtime_root = NULL) {
  if (!is.null(runtime_root)) {
    if (dir.exists(runtime_root)) {
      return(normalizePath(runtime_root, winslash = "/", mustWork = TRUE))
    }
    return(NA_character_)
  }
  candidates <- character()
  env_home <- Sys.getenv("GRAMMATEUS_HOME", unset = "")
  if (nzchar(env_home)) {
    candidates <- c(candidates, env_home)
  }
  option_home <- getOption("floundeR.grammateus_home", default = NULL)
  if (!is.null(option_home)) {
    candidates <- c(candidates, option_home)
  }
  candidates <- c(candidates, .grammateus_runtime_cache_candidates())
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

.grammateus_runtime_cache_candidates <- function() {
  cache_root <- tools::R_user_dir("floundeR", "cache")
  root <- file.path(cache_root, "grammateus")
  if (!dir.exists(root)) {
    return(character())
  }
  manifests <- list.files(
    root,
    pattern = "^manifest[.]json$",
    recursive = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )
  dirname(manifests)
}

.grammateus_runtime_validate_manifest_shape <- function(
    manifest,
    runtime_root,
    failures) {
  required <- c(
    "schema_version",
    "runtime_name",
    "runtime_version",
    "grammateus_release",
    "flounder_version_min",
    "flounder_version_max_exclusive",
    "platform",
    "abi",
    "artifacts",
    "capabilities",
    "signing",
    "provenance"
  )
  for (field in required) {
    if (is.null(manifest[[field]])) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "manifest_schema",
        "manifest.json",
        paste("Runtime manifest is missing required field:", field)
      )
    }
  }
  if (!identical(manifest$schema_version,
                 "flounder.grammateus_runtime_manifest.v1")) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "manifest_schema",
      "manifest.json",
      "Runtime manifest schema_version is unsupported."
    )
  }
  if (!identical(manifest$runtime_name, "grammateus-runtime")) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "manifest_schema",
      "manifest.json",
      "Runtime manifest runtime_name must be grammateus-runtime."
    )
  }
  signature_file <- manifest$signing$signature_file
  if (!is.character(signature_file) || length(signature_file) != 1L ||
      is.na(signature_file) || !nzchar(signature_file)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "manifest_signature",
      "manifest.json",
      "Runtime manifest does not declare a signature_file."
    )
  } else {
    failures <- .grammateus_runtime_validate_signature_file(
      runtime_root = runtime_root,
      signature_file = signature_file,
      expected_hash = manifest$signing$signature_sha256 %||% NA_character_,
      category = "manifest_signature",
      failures = failures
    )
  }
  failures
}

.grammateus_runtime_validate_compatibility <- function(manifest, failures) {
  platform <- manifest$platform %||% NA_character_
  current_platform <- R.version$platform
  if (!identical(platform, current_platform)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "platform_unsupported",
      "manifest.json",
      paste("Runtime platform", platform, "does not match", current_platform)
    )
  }
  version <- tryCatch(
    utils::packageVersion("floundeR"),
    error = function(error) numeric_version("0.0.0")
  )
  min_version <- .grammateus_runtime_numeric_version(
    manifest$flounder_version_min
  )
  max_version <- .grammateus_runtime_numeric_version(
    manifest$flounder_version_max_exclusive
  )
  if (is.null(min_version) || is.null(max_version) ||
      version < min_version || version >= max_version) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "version_incompatible",
      "manifest.json",
      "Runtime floundeR compatibility window does not include this package."
    )
  }
  failures
}

.grammateus_runtime_validate_compatibility_manifest <- function(
    manifest,
    runtime_root,
    failures) {
  runtime_version <- manifest$runtime_version %||% NA_character_
  compatibility_file <- manifest$compatibility_file %||%
    paste0("flounder-runtime-compatibility-", runtime_version, ".json")
  if (!.grammateus_runtime_safe_relative_path(compatibility_file)) {
    return(.grammateus_runtime_add_failure(
      failures,
      "compatibility_manifest",
      as.character(compatibility_file),
      "Runtime compatibility manifest path must be relative and must not escape root."
    ))
  }
  compatibility_path <- file.path(runtime_root, compatibility_file)
  if (!file.exists(compatibility_path)) {
    return(.grammateus_runtime_add_failure(
      failures,
      "compatibility_manifest",
      compatibility_file,
      "Runtime compatibility manifest is missing."
    ))
  }
  compatibility <- tryCatch(
    jsonlite::fromJSON(compatibility_path, simplifyVector = FALSE),
    error = function(error) {
      failures <<- .grammateus_runtime_add_failure(
        failures,
        "compatibility_manifest",
        compatibility_file,
        paste(
          "Runtime compatibility manifest is not valid JSON:",
          conditionMessage(error)
        )
      )
      NULL
    }
  )
  if (!is.list(compatibility)) {
    return(failures)
  }
  if (!identical(
    compatibility$schema_version,
    "flounder.grammateus_runtime_compatibility.v1"
  )) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "compatibility_manifest",
      compatibility_file,
      "Runtime compatibility manifest schema_version is unsupported."
    )
  }
  if (!identical(compatibility$runtime_version, runtime_version)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "version_incompatible",
      compatibility_file,
      "Runtime compatibility manifest version does not match manifest.json."
    )
  }
  manifest_min <- manifest$flounder_version_min %||% NA_character_
  manifest_max <- manifest$flounder_version_max_exclusive %||% NA_character_
  if (!identical(compatibility$flounder_version_min, manifest_min) ||
      !identical(compatibility$flounder_version_max_exclusive, manifest_max)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "version_incompatible",
      compatibility_file,
      "Runtime compatibility window disagrees with manifest.json."
    )
  }
  failures <- .grammateus_runtime_validate_compatibility_window(
    compatibility$flounder_version_min,
    compatibility$flounder_version_max_exclusive,
    compatibility_file,
    failures
  )
  required <- .grammateus_runtime_character_vector(
    compatibility$required_runtime_capabilities %||% character()
  )
  if (!is.character(required)) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "compatibility_manifest",
      compatibility_file,
      "Runtime compatibility manifest required_runtime_capabilities must be a character vector."
    )
  } else {
    failures <- .grammateus_runtime_validate_capabilities(
      manifest$capabilities %||% list(),
      required,
      failures
    )
  }
  signature_file <- compatibility$signing$signature_file %||%
    paste0(compatibility_file, ".sig")
  failures <- .grammateus_runtime_validate_signature_file(
    runtime_root = runtime_root,
    signature_file = signature_file,
    expected_hash = compatibility$signing$signature_sha256 %||% NA_character_,
    category = "compatibility_signature",
    failures = failures
  )
  failures
}

.grammateus_runtime_character_vector <- function(value) {
  if (is.character(value)) {
    return(value)
  }
  if (is.list(value) && all(vapply(value, function(item) {
    is.character(item) && length(item) == 1L && !is.na(item)
  }, logical(1L)))) {
    return(unlist(value, use.names = FALSE))
  }
  value
}

.grammateus_runtime_validate_compatibility_window <- function(
    min_version,
    max_version,
    path,
    failures) {
  version <- tryCatch(
    utils::packageVersion("floundeR"),
    error = function(error) numeric_version("0.0.0")
  )
  min_version <- .grammateus_runtime_numeric_version(min_version)
  max_version <- .grammateus_runtime_numeric_version(max_version)
  if (is.null(min_version) || is.null(max_version) ||
      version < min_version || version >= max_version) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "version_incompatible",
      path,
      "Runtime floundeR compatibility window does not include this package."
    )
  }
  failures
}

.grammateus_runtime_validate_capabilities <- function(
    capabilities,
    required_capabilities,
    failures) {
  for (capability in required_capabilities) {
    if (!isTRUE(capabilities[[capability]])) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "capability_missing",
        "manifest.json",
        paste("Runtime capability is missing or false:", capability)
      )
    }
  }
  failures
}

.grammateus_runtime_validate_artifacts <- function(
    artifacts,
    runtime_root,
    failures) {
  hashes <- character()
  if (!is.list(artifacts) || length(artifacts) == 0L) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      "manifest_schema",
      "manifest.json",
      "Runtime manifest must contain at least one artifact."
    )
    return(list(failures = failures, count = 0L, hashes = hashes))
  }

  for (artifact in artifacts) {
    path <- artifact$path %||% NA_character_
    expected_hash <- artifact$sha256 %||% NA_character_
    expected_bytes <- artifact$byte_len %||% NA_real_
    if (!.grammateus_runtime_safe_relative_path(path)) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "path_escape",
        as.character(path),
        "Runtime artifact path must be relative and must not escape root."
      )
      next
    }
    full_path <- file.path(runtime_root, path)
    if (!file.exists(full_path)) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "artifact_missing",
        path,
        "Runtime artifact is missing."
      )
      next
    }
    normalized_root <- normalizePath(runtime_root, winslash = "/", mustWork = TRUE)
    normalized_artifact <- normalizePath(full_path, winslash = "/", mustWork = TRUE)
    if (!startsWith(normalized_artifact, paste0(normalized_root, "/"))) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "path_escape",
        path,
        "Runtime artifact resolves outside the runtime root."
      )
      next
    }
    actual_hash <- unname(tools::sha256sum(full_path))
    actual_hash <- paste0("sha256:", tolower(actual_hash))
    hashes <- c(hashes, actual_hash)
    if (!identical(actual_hash, tolower(expected_hash))) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "checksum_mismatch",
        path,
        "Runtime artifact SHA-256 does not match manifest."
      )
    }
    actual_bytes <- unname(file.info(full_path)$size)
    if (is.na(expected_bytes) || !identical(as.numeric(expected_bytes),
                                           as.numeric(actual_bytes))) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "checksum_mismatch",
        path,
        "Runtime artifact byte length does not match manifest."
      )
    }
    signature_file <- artifact$signature_file %||% artifact$signature_asset %||%
      NA_character_
    if (is.character(signature_file) && length(signature_file) == 1L &&
        !is.na(signature_file) && nzchar(signature_file)) {
      failures <- .grammateus_runtime_validate_signature_file(
        runtime_root = runtime_root,
        signature_file = signature_file,
        expected_hash = artifact$signature_sha256 %||% NA_character_,
        category = "artifact_signature",
        failures = failures
      )
    }
  }
  list(failures = failures, count = length(artifacts), hashes = hashes)
}

.grammateus_runtime_validate_signature_file <- function(
    runtime_root,
    signature_file,
    expected_hash,
    category,
    failures) {
  if (!.grammateus_runtime_safe_relative_path(signature_file)) {
    return(.grammateus_runtime_add_failure(
      failures,
      category,
      as.character(signature_file),
      "Runtime signature path must be relative and must not escape root."
    ))
  }
  signature_path <- file.path(runtime_root, signature_file)
  if (!file.exists(signature_path)) {
    return(.grammateus_runtime_add_failure(
      failures,
      category,
      signature_file,
      "Runtime signature file is missing."
    ))
  }
  normalized_root <- normalizePath(runtime_root, winslash = "/", mustWork = TRUE)
  normalized_signature <- normalizePath(
    signature_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (!startsWith(normalized_signature, paste0(normalized_root, "/"))) {
    return(.grammateus_runtime_add_failure(
      failures,
      category,
      signature_file,
      "Runtime signature resolves outside the runtime root."
    ))
  }
  signature_bytes <- unname(file.info(signature_path)$size)
  if (is.na(signature_bytes) || signature_bytes <= 0) {
    failures <- .grammateus_runtime_add_failure(
      failures,
      category,
      signature_file,
      "Runtime signature file is empty."
    )
  }
  if (is.character(expected_hash) && length(expected_hash) == 1L &&
      !is.na(expected_hash) && nzchar(expected_hash)) {
    actual_hash <- paste0(
      "sha256:",
      tolower(unname(tools::sha256sum(signature_path)))
    )
    expected_hash <- tolower(expected_hash)
    if (!grepl("^sha256:", expected_hash)) {
      expected_hash <- paste0("sha256:", expected_hash)
    }
    if (!identical(actual_hash, expected_hash)) {
      failures <- .grammateus_runtime_add_failure(
        failures,
        "checksum_mismatch",
        signature_file,
        "Runtime signature SHA-256 does not match manifest."
      )
    }
  }
  failures
}

.grammateus_runtime_safe_relative_path <- function(path) {
  is.character(path) &&
    length(path) == 1L &&
    !is.na(path) &&
    nzchar(path) &&
    !grepl("^(/|[A-Za-z]:)", path) &&
    !any(strsplit(path, "/", fixed = TRUE)[[1L]] %in% c("..", ""))
}

.grammateus_runtime_validation <- function(
    available,
    valid,
    runtime_root,
    runtime_version,
    platform,
    capabilities,
    artifact_count,
    artifact_sha256,
    failures) {
  structure(
    list(
      schema_version = "flounder.grammateus_runtime_validation.v1",
      available = isTRUE(available),
      valid = isTRUE(valid),
      runtime_root = runtime_root,
      runtime_version = runtime_version,
      platform = platform,
      capabilities = capabilities,
      artifact_count = as.integer(artifact_count),
      artifact_sha256 = artifact_sha256,
      failures = failures
    ),
    class = c("flounder_grammateus_runtime_validation", "list")
  )
}

.grammateus_runtime_failures <- function() {
  data.frame(
    category = character(),
    path = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

.grammateus_runtime_add_failure <- function(failures, category, path, message) {
  rbind(
    failures,
    data.frame(
      category = category,
      path = path,
      message = message,
      stringsAsFactors = FALSE
    )
  )
}

.grammateus_runtime_aggregate_hash <- function(hashes) {
  if (length(hashes) == 0L) {
    return(NA_character_)
  }
  path <- tempfile("flounder-grammateus-runtime-hashes-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(sort(hashes), path, useBytes = TRUE)
  paste0("sha256:", unname(tools::sha256sum(path)))
}

.grammateus_runtime_numeric_version <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    return(NULL)
  }
  tryCatch(numeric_version(value), error = function(error) NULL)
}

.grammateus_required_runtime_root <- function(
    path,
    field,
    must_exist = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    .grammateus_runtime_abort(
      paste(field, "must be a non-empty character scalar."),
      "flounder_grammateus_runtime_argument_error"
    )
  }
  if (must_exist && !dir.exists(path)) {
    .grammateus_runtime_abort(
      paste(field, "does not exist:", path),
      "flounder_grammateus_runtime_argument_error"
    )
  }
  normalizePath(path, winslash = "/", mustWork = must_exist)
}

.grammateus_runtime_abort <- function(message, class, validation = NULL) {
  stop(
    structure(
      list(message = message, validation = validation, call = NULL),
      class = c(class, "flounder_grammateus_error", "error", "condition")
    )
  )
}
