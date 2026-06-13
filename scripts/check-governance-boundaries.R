#!/usr/bin/env Rscript

read_text <- function(path) {
  if (!file.exists(path)) {
    stop("Required governance file is missing: ", path, call. = FALSE)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

require_pattern <- function(text, pattern, label, fixed = FALSE) {
  if (fixed) {
    found <- grepl(tolower(pattern), tolower(text), fixed = TRUE)
  } else {
    found <- grepl(pattern, text, perl = TRUE, ignore.case = TRUE)
  }
  if (!found) {
    stop("Governance boundary check failed: ", label, call. = FALSE)
  }
}

distribution <- read_text("DISTRIBUTION.md")
governance <- read_text("GOVERNANCE.md")
agents <- read_text("AGENTS.md")
roadmap <- read_text("ROADMAP.md")

for (open_project in c("floundeR", "pod5-tools", "bamana", "porkchop")) {
  require_pattern(
    distribution,
    paste0(open_project, ".*open|open.*", open_project),
    paste0(open_project, " must be documented as open-source/open")
  )
}

require_pattern(
  distribution,
  "Grammateus.*private|private.*Grammateus",
  "Grammateus must be documented as private"
)
require_pattern(
  distribution,
  "Grammateus.*optional|optional.*Grammateus",
  "Grammateus must be documented as optional from the public package perspective"
)
require_pattern(
  distribution,
  "prebuilt.*Grammateus|Grammateus.*prebuilt",
  "Grammateus distribution must use prebuilt runtime/reporting artifacts"
)

for (porkchop_doc in c(
  "../porkchop/AGENTS.md",
  "../porkchop/ROADMAP.md",
  "../porkchop/README.md",
  "../porkchop/docs/output/json.rst",
  "../porkchop/docs/kits/provenance.rst",
  "../porkchop/docs/validation/index.rst"
)) {
  require_pattern(
    governance,
    porkchop_doc,
    paste0("Porkchop governing document must be listed: ", porkchop_doc),
    fixed = TRUE
  )
}

require_pattern(
  agents,
  "Do not mirror every `pod5-tools`, `bamana`, or `porkchop` capability",
  "AGENTS.md must preserve the curated integration rule"
)
require_pattern(
  roadmap,
  "Distribution boundary: `floundeR`, `pod5-tools`, `bamana`, and `porkchop` are\\s+open-source",
  "ROADMAP.md must preserve the open/private distribution boundary"
)

message("Governance boundary checks passed.")
