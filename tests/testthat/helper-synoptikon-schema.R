synoptikon_payload_schema_path <- function() {
  path <- file.path(
    synoptikon_flounder_root(),
    "inst", "schema", "synoptikon-qc-payload-v1.schema.json")
  if (!file.exists(path)) {
    path <- system.file(
      "schema", "synoptikon-qc-payload-v1.schema.json",
      package = "floundeR", mustWork = FALSE)
  }
  path
}

synoptikon_mnemosyne_schema_candidates <- function() {
  env <- Sys.getenv("FLOUNDER_MNEMOSYNE_SYNOPTIKON_SCHEMA", unset = NA_character_)
  root <- normalizePath(file.path(synoptikon_flounder_root(), ".."), mustWork = FALSE)
  unique(stats::na.omit(c(
    env,
    file.path(root, "mnemosyne", "schemas", "flounder", "synoptikon-qc-payload-v1.schema.json"),
    file.path(root, "mnemosyne", "mneion-api-types", "schemas", "flounder", "synoptikon-qc-payload-v1.schema.json"),
    file.path(root, "mnemosyne-docs", "products", "mneion", "schemas", "flounder", "synoptikon-qc-payload-v1.schema.json")
  )))
}

synoptikon_flounder_root <- function() {
  path <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    description <- file.path(path, "DESCRIPTION")
    if (file.exists(description)) {
      package <- tryCatch(
        read.dcf(description, fields = "Package")[[1]],
        error = function(e) NA_character_)
      if (identical(package, "floundeR")) {
        return(path)
      }
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      break
    }
    path <- parent
  }
  normalizePath(getwd(), mustWork = TRUE)
}

synoptikon_mnemosyne_schema_path <- function() {
  found <- synoptikon_mnemosyne_schema_candidates()
  found <- found[file.exists(found)]
  if (length(found) == 0) {
    return(NULL)
  }
  found[[1]]
}

expect_synoptikon_schema_valid <- function(instance, schema_path = synoptikon_payload_schema_path()) {
  schema <- jsonlite::fromJSON(schema_path, simplifyVector = FALSE)
  errors <- synoptikon_validate_schema_node(instance, schema, schema, "$")
  if (length(errors) > 0) {
    fail(paste(
      "Synoptikon payload did not validate against JSON schema:",
      paste(utils::head(errors, 20), collapse = "\n"),
      sep = "\n"
    ))
  }
  succeed()
}

synoptikon_validate_schema_node <- function(value, schema, root_schema, path) {
  if (!is.null(schema[["$ref"]])) {
    schema <- synoptikon_resolve_schema_ref(schema[["$ref"]], root_schema)
  }

  errors <- character()
  if (!is.null(schema$type) &&
      !synoptikon_schema_type_matches(value, schema$type)) {
    errors <- c(
      errors,
      sprintf(
        "%s expected type %s, got %s",
        path,
        paste(schema$type, collapse = " or "),
        synoptikon_json_type(value)
      )
    )
    return(errors)
  }

  if (!is.null(schema$const) && !identical(value, schema$const)) {
    errors <- c(errors, sprintf("%s expected const %s", path, schema$const))
  }
  if (!is.null(schema$enum) &&
      !any(vapply(schema$enum, identical, logical(1), y = value))) {
    errors <- c(errors, sprintf("%s value is not in enum", path))
  }
  if (!is.null(schema$minLength) && is.character(value) &&
      nchar(value, type = "chars") < schema$minLength) {
    errors <- c(errors, sprintf("%s is shorter than minLength", path))
  }
  if (!is.null(schema$minimum) && is.numeric(value) && value < schema$minimum) {
    errors <- c(errors, sprintf("%s is less than minimum", path))
  }
  if (!is.null(schema$pattern) && is.character(value) &&
      !grepl(schema$pattern, value, perl = TRUE)) {
    errors <- c(errors, sprintf("%s does not match pattern", path))
  }
  if (!is.null(schema$format) && identical(schema$format, "date-time") &&
      is.character(value) && !synoptikon_is_rfc3339_datetime(value)) {
    errors <- c(errors, sprintf("%s is not an RFC3339 date-time", path))
  }

  if (synoptikon_json_type(value) == "object") {
    names_value <- names(value) %||% character()
    required <- schema$required %||% character()
    missing <- setdiff(required, names_value)
    if (length(missing) > 0) {
      errors <- c(errors, sprintf("%s missing required property %s", path, missing))
    }

    properties <- schema$properties %||% list()
    for (name in intersect(names(properties), names_value)) {
      errors <- c(
        errors,
        synoptikon_validate_schema_node(
          value[[name]], properties[[name]], root_schema,
          paste0(path, ".", name)
        )
      )
    }

    extra <- setdiff(names_value, names(properties))
    if (length(extra) > 0) {
      additional <- schema$additionalProperties
      if (identical(additional, FALSE)) {
        errors <- c(errors, sprintf("%s has undeclared property %s", path, extra))
      } else if (is.list(additional)) {
        for (name in extra) {
          errors <- c(
            errors,
            synoptikon_validate_schema_node(
              value[[name]], additional, root_schema,
              paste0(path, ".", name)
            )
          )
        }
      }
    }
  }

  if (synoptikon_json_type(value) == "array" && !is.null(schema$items)) {
    for (i in seq_along(value)) {
      errors <- c(
        errors,
        synoptikon_validate_schema_node(
          value[[i]], schema$items, root_schema,
          paste0(path, "[", i, "]")
        )
      )
    }
    if (!is.null(schema$minItems) && length(value) < schema$minItems) {
      errors <- c(errors, sprintf("%s has fewer items than minItems", path))
    }
  }

  errors
}

synoptikon_resolve_schema_ref <- function(ref, root_schema) {
  if (!startsWith(ref, "#/")) {
    stop("Only local JSON schema references are supported in tests.", call. = FALSE)
  }
  parts <- strsplit(sub("^#/", "", ref), "/", fixed = TRUE)[[1]]
  node <- root_schema
  for (part in parts) {
    node <- node[[part]]
  }
  node
}

synoptikon_schema_type_matches <- function(value, types) {
  any(vapply(types, function(type) {
    switch(
      type,
      "null" = is.null(value),
      "boolean" = is.logical(value) && length(value) == 1 && !is.na(value),
      "string" = is.character(value) && length(value) == 1 && !is.na(value),
      "integer" = is.numeric(value) && length(value) == 1 &&
        !is.na(value) && value == floor(value),
      "number" = is.numeric(value) && length(value) == 1 && !is.na(value),
      "object" = is.list(value) && !is.null(names(value)),
      "array" = is.list(value) && is.null(names(value)),
      FALSE
    )
  }, logical(1)))
}

synoptikon_json_type <- function(value) {
  if (is.null(value)) {
    return("null")
  }
  if (is.logical(value)) {
    return("boolean")
  }
  if (is.character(value)) {
    return("string")
  }
  if (is.numeric(value) && value == floor(value)) {
    return("integer")
  }
  if (is.numeric(value)) {
    return("number")
  }
  if (is.list(value) && !is.null(names(value))) {
    return("object")
  }
  if (is.list(value)) {
    return("array")
  }
  typeof(value)
}

synoptikon_is_rfc3339_datetime <- function(value) {
  grepl(
    paste0(
      "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}",
      "(\\.\\d+)?(Z|[+-]\\d{2}:\\d{2})$"
    ),
    value,
    perl = TRUE
  )
}
