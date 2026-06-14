#' List objects in an ONT open-data S3 prefix
#'
#' `ont_open_data_list()` performs an explicit anonymous listing of a public ONT
#' open-data S3 prefix. The default prefix is the selected Zymo fecal POD5
#' source used by floundeR integration examples. The helper lists metadata only;
#' it does not download POD5 files.
#'
#' @param prefix Character scalar. S3 object prefix to list.
#' @param bucket Character scalar. Public S3 bucket name.
#' @param region Character scalar. AWS region for the bucket.
#' @param max Integer scalar. Maximum number of objects to return.
#' @param anonymous Logical scalar. If `TRUE`, temporarily sets
#'   `AWS_NO_SIGN_REQUEST=true` while listing.
#'
#' @return A data frame with one row per object and the columns `bucket`, `key`,
#'   `size`, `last_modified_utc`, `etag`, and `storage_class`.
#'
#' @examples
#' \dontrun{
#' if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
#'   ont_open_data_list(max = 5)
#' }
#' }
#'
#' @export
ont_open_data_list <- function(
    prefix = flounder_ont_zymo_pod5_prefix(),
    bucket = "ont-open-data",
    region = "eu-west-1",
    max = 1000L,
    anonymous = TRUE) {
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    .flounder_open_data_error("`prefix` must be a non-missing character scalar.")
  }
  if (!is.character(bucket) || length(bucket) != 1L || is.na(bucket)) {
    .flounder_open_data_error("`bucket` must be a non-missing character scalar.")
  }
  if (!is.character(region) || length(region) != 1L || is.na(region)) {
    .flounder_open_data_error("`region` must be a non-missing character scalar.")
  }
  if (!is.numeric(max) || length(max) != 1L || is.na(max) || max < 1) {
    .flounder_open_data_error("`max` must be a positive numeric scalar.")
  }
  if (!is.logical(anonymous) || length(anonymous) != 1L || is.na(anonymous)) {
    .flounder_open_data_error("`anonymous` must be TRUE or FALSE.")
  }

  old_no_sign <- Sys.getenv("AWS_NO_SIGN_REQUEST", unset = NA_character_)
  if (isTRUE(anonymous)) {
    Sys.setenv(AWS_NO_SIGN_REQUEST = "true")
    on.exit(.flounder_restore_env("AWS_NO_SIGN_REQUEST", old_no_sign), add = TRUE)
  }

  objects <- tryCatch(
    .flounder_ont_open_data_get_bucket(
      bucket = bucket,
      prefix = prefix,
      max = as.integer(max),
      region = region
    ),
    error = function(error) {
      .flounder_open_data_error(
        paste("Failed to list ONT open-data objects:", conditionMessage(error)),
        parent = error
      )
    }
  )

  .flounder_ont_open_data_listing(objects, bucket = bucket)
}

#' @rdname ont_open_data_list
#' @export
flounder_ont_zymo_pod5_prefix <- function() {
  "zymo_fecal_2025.05/raw/PAU85136/pod5/"
}

.flounder_ont_open_data_get_bucket <- function(bucket, prefix, max, region) {
  aws.s3::get_bucket(
    bucket = bucket,
    prefix = prefix,
    max = max,
    region = region
  )
}

.flounder_ont_open_data_listing <- function(objects, bucket) {
  if (length(objects) == 0L) {
    return(data.frame(
      bucket = character(),
      key = character(),
      size = numeric(),
      last_modified_utc = character(),
      etag = character(),
      storage_class = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(objects, function(object) {
    data.frame(
      bucket = bucket,
      key = .flounder_aws_field(object, "Key"),
      size = as.numeric(.flounder_aws_field(object, "Size")),
      last_modified_utc = .flounder_aws_field(object, "LastModified"),
      etag = .flounder_aws_field(object, "ETag"),
      storage_class = .flounder_aws_field(object, "StorageClass"),
      stringsAsFactors = FALSE
    )
  })

  listing <- do.call(rbind, rows)
  row.names(listing) <- NULL
  listing
}

.flounder_aws_field <- function(object, name) {
  value <- object[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_character_)
  }
  as.character(value[[1L]])
}

.flounder_restore_env <- function(name, value) {
  if (is.na(value)) {
    Sys.unsetenv(name)
  } else {
    args <- list(value)
    names(args) <- name
    do.call(Sys.setenv, args)
  }
}

.flounder_open_data_error <- function(message, parent = NULL, call = NULL) {
  condition <- structure(
    list(message = message, parent = parent, call = call),
    class = c(
      "floundeR_open_data_error",
      "floundeR_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}
