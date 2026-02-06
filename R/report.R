# Reporting helpers

#' Save Results, Plots, and Manifest
#'
#' @description
#' Saves tables and plots from a \code{profile_result} object and writes a
#' manifest with spec and session info.
#'
#' @param result A \code{profile_result} object.
#' @param out_dir Output directory.
#' @param report_spec List with options: \code{table_format} (csv/tsv),
#'   \code{plot_format} (png/pdf), \code{write_manifest} (TRUE/FALSE).
#'
#' @return Invisibly returns a manifest list.
#' @export
make_report <- function(result, out_dir, report_spec = list()) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  table_format <- report_spec$table_format %||% "csv"
  plot_format <- report_spec$plot_format %||% "png"
  write_manifest <- isTRUE(report_spec$write_manifest %||% TRUE)

  if (!is.null(result$data)) {
    fn <- file.path(out_dir, paste0("results.", table_format))
    if (table_format == "tsv") {
      utils::write.table(result$data, fn, sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      utils::write.csv(result$data, fn, row.names = FALSE)
    }
  }

  if (!is.null(result$diagnostics)) {
    fn <- file.path(out_dir, paste0("diagnostics.", table_format))
    if (table_format == "tsv") {
      utils::write.table(result$diagnostics, fn, sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      utils::write.csv(result$diagnostics, fn, row.names = FALSE)
    }
  }

  if (!is.null(result$summaries)) {
    fn <- file.path(out_dir, paste0("summaries.", table_format))
    if (table_format == "tsv") {
      utils::write.table(result$summaries, fn, sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      utils::write.csv(result$summaries, fn, row.names = FALSE)
    }
  }

  if (!is.null(result$plots)) {
    .check_pkg("ggplot2", "plot saving")
    for (nm in names(result$plots)) {
      plot_obj <- result$plots[[nm]]
      fn <- file.path(out_dir, paste0(nm, ".", plot_format))
      ggplot2::ggsave(fn, plot = plot_obj, bg = "white")
    }
  }

  manifest <- list(
    timestamp = Sys.time(),
    spec = result$meta$spec %||% list(),
    session = sessionInfo(),
    log = result$log
  )

  if (write_manifest) {
    fn <- file.path(out_dir, "manifest.json")
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::write_json(manifest, fn, auto_unbox = TRUE, pretty = TRUE)
    } else {
      utils::write.table(capture.output(str(manifest)), fn, row.names = FALSE, col.names = FALSE)
    }
  }
  invisible(manifest)
}
