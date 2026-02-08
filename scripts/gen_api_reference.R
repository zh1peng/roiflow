#!/usr/bin/env Rscript

# Generate a single markdown file listing Usage + Arguments for all exported
# functions, based on the package's Rd files under man/.
#
# Output: inst/examples/api_reference.md

out_path <- file.path("inst", "examples", "api_reference.md")
rd_files <- sort(list.files("man", pattern = "\\.Rd$", full.names = TRUE))

`%||%` <- function(x, y) if (is.null(x)) y else x

.rd_tag <- function(x) attr(x, "Rd_tag")

.flatten <- function(node) {
  if (is.null(node)) return("")
  if (is.character(node)) return(paste(node, collapse = ""))
  if (!is.list(node)) return("")

  tag <- .rd_tag(node)
  txt <- paste(vapply(node, .flatten, character(1)), collapse = "")

  if (!is.null(tag) && tag %in% c("\\code", "\\verb")) {
    return(txt)
  }
  txt
}

.extract_section <- function(rd, tag) {
  tags <- vapply(rd, .rd_tag, character(1))
  idx <- which(tags == tag)
  if (length(idx) == 0) return(NULL)
  rd[[idx[1]]]
}

.extract_first_text <- function(rd, tag) {
  sec <- .extract_section(rd, tag)
  if (is.null(sec)) return(NULL)
  s <- trimws(.flatten(sec))
  if (!nzchar(s)) return(NULL)
  s
}

.extract_usage <- function(rd) {
  sec <- .extract_section(rd, "\\usage")
  if (is.null(sec)) return(NULL)
  s <- trimws(paste(as.character(sec), collapse = ""))
  if (!nzchar(s)) return(NULL)
  s
}

.extract_arguments <- function(rd) {
  sec <- .extract_section(rd, "\\arguments")
  if (is.null(sec)) return(list())

  items <- sec[vapply(sec, function(x) identical(.rd_tag(x), "\\item"), logical(1))]
  out <- list()
  for (it in items) {
    if (!is.list(it) || length(it) < 2) next
    nm <- trimws(.flatten(it[[1]]))
    desc <- trimws(.flatten(it[[2]]))
    desc <- gsub("[[:space:]]+", " ", desc)
    if (!nzchar(nm)) next
    out[[nm]] <- desc
  }
  out
}

lines <- character(0)
lines <- c(lines, "# roiflow API: Usage + Arguments")
lines <- c(lines, "")
lines <- c(lines, "This reference is generated from `man/*.Rd` and shows the exact function signatures and documented parameters.")
lines <- c(lines, "")

for (f in rd_files) {
  rd <- tools::parse_Rd(f)
  name <- .extract_first_text(rd, "\\name") %||% tools::file_path_sans_ext(basename(f))
  title <- .extract_first_text(rd, "\\title")
  usage <- .extract_usage(rd)
  args <- .extract_arguments(rd)

  lines <- c(lines, sprintf("## `%s`", name))
  if (!is.null(title)) lines <- c(lines, "", title)

  if (!is.null(usage)) {
    lines <- c(lines, "", "**Usage**", "")
    lines <- c(lines, "```r", usage, "```")
  }

  if (length(args) > 0) {
    lines <- c(lines, "", "**Arguments**", "")
    for (nm in names(args)) {
      lines <- c(lines, sprintf("- `%s`: %s", nm, args[[nm]]))
    }
  } else {
    lines <- c(lines, "", "**Arguments**", "", "- (No documented arguments.)")
  }

  lines <- c(lines, "")
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, out_path, useBytes = TRUE)
cat(sprintf("Wrote: %s\n", out_path))

