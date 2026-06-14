# Render Grammateus report elements through the Rust binding

`grammateus_render_element()` sends a prepared floundeR Grammateus
semantic element to the compiled Rust report-rendering boundary.
`format = "html"` returns rendered HTML as a character scalar, while
`format = "pdf"` returns PDF bytes as a raw vector. Public floundeR
builds do not bundle the private Grammateus runtime; when it is not
linked these functions fail with a typed
`flounder_grammateus_runtime_unavailable` condition.

## Usage

``` r
grammateus_render_element(element, format = c("html", "pdf"))

grammateus_render_figure_html(element)

grammateus_render_figure_pdf(element)
```

## Arguments

- element:

  A Grammateus semantic element prepared by floundeR, currently a
  `flounder_grammateus_figure`.

- format:

  Output format: `html` or `pdf`.

## Value

Rendered HTML as a character scalar for `format = "html"`, or PDF bytes
as a raw vector for `format = "pdf"`.
