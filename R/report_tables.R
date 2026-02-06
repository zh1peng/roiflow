# Results table formatting + export helpers

#' Format a Results Table
#'
#' @description
#' Formats a standardized results table (e.g., from \code{\link{compare_groups}})
#' into a stable, publication-friendly layout. This function is deterministic:
#' no sorting, filtering, or p-adjustment is performed here.
#'
#' @param res_tbl Data frame. Standard results table.
#' @param style Character. One of \code{"minimal"} or \code{"publication"}.
#' @param roi_drop_suffix Optional character. If provided, remove this suffix from
#'   ROI names (e.g., \code{"_thickavg"}).
#' @param roi_drop_prefix Optional character. Regex prefix to remove from ROI
#'   names (e.g., \code{"^SV_"}).
#' @param digits Integer. Rounding digits for effect/estimate columns.
#' @param p_digits Integer. Significant digits for p-values.
#'
#' @return A formatted data frame.
#' @export
format_results_table <- function(res_tbl,
                                 style = c("minimal", "publication"),
                                 roi_drop_suffix = "_thickavg",
                                 roi_drop_prefix = "^SV_",
                                 digits = 3,
                                 p_digits = 3) {
  style <- match.arg(style)
  stopifnot(is.data.frame(res_tbl))

  required <- c("var", "estimate", "std_error", "statistic", "p_value", "p_adj",
                "ref_level", "target_level", "es_value", "es_method")
  missing <- setdiff(required, names(res_tbl))
  if (length(missing) > 0) {
    stop(sprintf("E_BAD_RESULTS_SCHEMA: Missing required columns in res_tbl: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  .sig_star <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
  }

  out <- res_tbl

  out$ROI <- out$var
  if (!is.null(roi_drop_suffix) && nzchar(roi_drop_suffix)) {
    out$ROI <- gsub(paste0(roi_drop_suffix, "$"), "", out$ROI)
  }
  if (!is.null(roi_drop_prefix) && nzchar(roi_drop_prefix)) {
    out$ROI <- gsub(roi_drop_prefix, "", out$ROI)
  }

  out$Sig <- vapply(out$p_adj, .sig_star, character(1))

  # stable rounding (keep raw columns intact in minimal style)
  if (style == "minimal") {
    out$estimate_round <- round(out$estimate, digits)
    out$std_error_round <- round(out$std_error, digits)
    out$statistic_round <- round(out$statistic, digits)
    out$p_value_round <- signif(out$p_value, p_digits)
    out$p_adj_round <- signif(out$p_adj, p_digits)
    out$es_value_round <- round(out$es_value, digits)
    if ("es_ci_low" %in% names(out)) out$es_ci_low_round <- round(out$es_ci_low, digits)
    if ("es_ci_high" %in% names(out)) out$es_ci_high_round <- round(out$es_ci_high, digits)
    return(out)
  }

  # publication layout: compact, consistent columns
  pub <- data.frame(
    ROI = out$ROI,
    term = if ("term" %in% names(out)) out$term else NA_character_,
    beta = round(out$estimate, digits),
    se = round(out$std_error, digits),
    statistic = round(out$statistic, digits),
    p = signif(out$p_value, p_digits),
    FDR_p = signif(out$p_adj, p_digits),
    Sig = out$Sig,
    EffectSize = round(out$es_value, digits),
    ES_method = out$es_method,
    ES_CI_low = if ("es_ci_low" %in% names(out)) round(out$es_ci_low, digits) else NA_real_,
    ES_CI_high = if ("es_ci_high" %in% names(out)) round(out$es_ci_high, digits) else NA_real_,
    Ref = out$ref_level,
    Target = out$target_level,
    contrast_label = if ("contrast_label" %in% names(out)) out$contrast_label else NA_character_,
    stringsAsFactors = FALSE
  )

  pub
}

.write_table <- function(x, path, format = c("csv", "tsv")) {
  format <- match.arg(format)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (format == "tsv") {
    utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

#' Export Raw + Formatted Results Tables
#'
#' @description
#' Writes a raw results table and a formatted (publication-style) table to disk.
#'
#' @param res_tbl Data frame. Raw results table.
#' @param out_dir Character. Output directory for tables.
#' @param base_name Character. File prefix (without extension).
#' @param formats Character vector. Any of \code{"csv"} or \code{"tsv"}.
#' @param style Character. Formatting style for the formatted table.
#' @param overwrite Logical. If \code{FALSE} (default), error if outputs exist.
#' @param ... Passed to \code{\link{format_results_table}}.
#'
#' @return Invisibly returns a list of written file paths.
#' @export
export_results_tables <- function(res_tbl, out_dir, base_name = "results",
                                 formats = c("csv"),
                                 style = c("publication", "minimal"),
                                 overwrite = FALSE,
                                 ...) {
  style <- match.arg(style)
  stopifnot(is.data.frame(res_tbl), is.character(out_dir), length(out_dir) == 1)
  stopifnot(is.character(base_name), length(base_name) == 1)
  formats <- unique(formats)
  bad <- setdiff(formats, c("csv", "tsv"))
  if (length(bad) > 0) stop("Unsupported formats: ", paste(bad, collapse = ", "), call. = FALSE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  raw_tbl <- res_tbl
  fmt_tbl <- format_results_table(res_tbl, style = style, ...)

  written <- list(raw = character(0), formatted = character(0))
  for (fmt in formats) {
    raw_path <- file.path(out_dir, paste0(base_name, "_raw.", fmt))
    fmt_path <- file.path(out_dir, paste0(base_name, "_formatted.", fmt))
    if (!overwrite && (file.exists(raw_path) || file.exists(fmt_path))) {
      stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s",
                   paste(c(raw_path, fmt_path), collapse = ", ")), call. = FALSE)
    }
    .write_table(raw_tbl, raw_path, format = fmt)
    .write_table(fmt_tbl, fmt_path, format = fmt)
    written$raw <- c(written$raw, raw_path)
    written$formatted <- c(written$formatted, fmt_path)
  }

  invisible(written)
}

