library(floundeR)

local_fixture_path <- function(...) {
  candidates <- c(
    file.path("fixtures", ...),
    file.path("tests", "testthat", "fixtures", ...)
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(existing[[1L]])
  }
  candidates[[1L]]
}

local_without_grammateus_runtime <- function(code) {
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
  Sys.setenv(GRAMMATEUS_HOME = tempfile("missing-grammateus-env-"))
  options(floundeR.grammateus_home = tempfile("missing-grammateus-option-"))
  force(code)
}

test_that("core QC and Synoptikon APIs do not require Grammateus assets", {
  local_without_grammateus_runtime({
    expect_false(grammateus_runtime_available(tempfile("missing-runtime-")))

    summary_path <- local_fixture_path("sequencing_summary_dorado.tsv")
    summary <- qc_run_summary(summary_path, source_id = "no-grammateus")
    yield <- qc_yield_over_time(summary_path, resolution_minutes = 1)
    read_lengths <- qc_read_length_distribution(summary_path, bins = 3)
    qualities <- qc_quality_distribution(summary_path, bins = 3)
    channels <- qc_channel_density(summary_path)
    barcodes <- qc_barcode_composition(summary_path)
    card <- qc_report_card(summary, barcode_composition = barcodes)
    output <- tempfile(fileext = ".json")

    payload_path <- write_synoptikon_qc(
      output,
      run_summary = summary,
      flowcell = list(channel_density = channels),
      barcode = list(composition = barcodes),
      report_cards = card,
      payload_id = "no-grammateus-payload",
      generated_at_utc = "2026-06-14T10:57:08Z"
    )
    payload <- jsonlite::fromJSON(payload_path, simplifyVector = FALSE)

    expect_equal(summary$schema_version, "flounder.qc_run_summary.v1")
    expect_equal(unique(yield$schema_version),
                 "flounder.qc_yield_over_time.v1")
    expect_equal(unique(read_lengths$schema_version),
                 "flounder.qc_read_length_distribution.v1")
    expect_equal(unique(qualities$schema_version),
                 "flounder.qc_quality_distribution.v1")
    expect_equal(unique(channels$schema_version),
                 "flounder.qc_channel_density.v1")
    expect_equal(unique(barcodes$schema_version),
                 "flounder.qc_barcode_composition.v1")
    expect_equal(unique(card$schema_version), "flounder.qc_report_card.v1")
    expect_equal(payload$schema_version, "flounder.synoptikon_qc_payload.v1")
    expect_true(payload$qc_sections$sequencing_summary$present)
    expect_true(file.exists(payload_path))
  })
})

test_that("POD5 discovery does not require Grammateus assets", {
  skip_if_no_flounder_rust("POD5 discovery without Grammateus runtime")

  local_without_grammateus_runtime({
    root <- tempfile("no-grammateus-pod5-find-")
    dir.create(root)

    discovered <- pod5_find(root)

    expect_s3_class(discovered, "data.frame")
    expect_equal(nrow(discovered), 0L)
  })
})

test_that("report rendering remains runtime-gated when Grammateus is absent", {
  skip_if_no_flounder_rust("Grammateus render binding without runtime")

  local_without_grammateus_runtime({
    image_path <- tempfile("no-grammateus-render-", fileext = ".png")
    grDevices::png(image_path, width = 24, height = 24)
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    grDevices::dev.off()
    on.exit(unlink(image_path), add = TRUE)

    figure <- grammateus_figure_from_file(
      path = image_path,
      figure_id = "figure_no_grammateus_runtime",
      caption = "No-runtime rendering boundary smoke test.",
      alt_text = "Blank image for no-runtime rendering boundary test."
    )

    expect_error(
      grammateus_render_figure_html(figure),
      class = "flounder_grammateus_runtime_unavailable"
    )
  })
})
