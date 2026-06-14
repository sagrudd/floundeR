test_that("Grammateus render bindings report unavailable runtime in public builds", {
  skip_if_no_flounder_rust("Grammateus Rust render binding")

  image_path <- tempfile("flounder-grammateus-render-", fileext = ".png")
  grDevices::png(image_path, width = 24, height = 24)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  grDevices::dev.off()
  on.exit(unlink(image_path), add = TRUE)

  figure <- floundeR::grammateus_figure_from_file(
    path = image_path,
    figure_id = "figure_render_binding",
    caption = "Render binding smoke test.",
    alt_text = "Blank render binding smoke-test image."
  )

  capabilities <- floundeR::flounder_rust_capabilities()
  expect_identical(capabilities$grammateus, "not_linked")
  expect_error(
    floundeR::grammateus_render_figure_html(figure),
    class = "flounder_grammateus_runtime_unavailable"
  )
  expect_error(
    floundeR::grammateus_render_figure_pdf(figure),
    class = "flounder_grammateus_runtime_unavailable"
  )
  expect_error(
    floundeR::grammateus_render_element(figure, format = "html"),
    class = "flounder_grammateus_runtime_unavailable"
  )
})

test_that("Grammateus render dispatcher rejects unsupported elements", {
  expect_error(
    floundeR::grammateus_render_element(list(kind = "table"), format = "html"),
    class = "flounder_grammateus_unsupported_element"
  )
})
