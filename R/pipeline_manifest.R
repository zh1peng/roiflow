# Manifest + logging helpers for config-first pipeline

.pkg_version <- function(pkg = "roiflow") {
  v <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NULL)
  if (!is.null(v)) return(v)
  dcf <- tryCatch(read.dcf("DESCRIPTION"), error = function(e) NULL)
  if (!is.null(dcf) && "Version" %in% colnames(dcf)) return(as.character(dcf[1, "Version"]))
  NA_character_
}

.git_commit <- function() {
  sha <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) character(0))
  if (length(sha) == 0) return(NA_character_)
  sha <- sha[1]
  if (!is.character(sha) || length(sha) != 1) return(NA_character_)
  sha
}

.file_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)[1])
}

.sanitize_for_json <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = ""))
  if (inherits(x, "POSIXt")) return(as.character(x))
  if (is.factor(x)) return(as.character(x))
  if (is.function(x)) return("<function>")
  if (is.environment(x)) return("<environment>")
  if (is.list(x)) {
    out <- lapply(x, .sanitize_for_json)
    # preserve names, but drop empty placeholders
    if (!is.null(names(out))) {
      empty <- vapply(out, function(z) is.null(z) || (is.character(z) && length(z) == 1 && z == ""), logical(1))
      out[!empty]
    } else {
      out
    }
  } else {
    x
  }
}

#' Build a Run Manifest From Pipeline Context
#'
#' @description
#' Creates a manifest list derived from the resolved config, actual spec objects,
#' step timings, and results returned by the pipeline. This function does not
#' write to disk.
#'
#' @param ctx Pipeline context returned by \code{\link{run_pipeline}}.
#'
#' @return A nested list suitable for JSON serialization.
#' @export
build_manifest <- function(ctx) {
  if (!is.list(ctx)) stop("ctx must be a list.", call. = FALSE)

  cfg <- ctx$cfg %||% list()
  input <- cfg$io$input %||% NA_character_
  out_root <- cfg$io$out_root %||% NA_character_

  df_in <- ctx$input$data %||% NULL
  df_final <- ctx$data %||% NULL

  analysis_tbl <- ctx$results$compare_groups$data %||% ctx$results$analysis$data %||% NULL

  # small derived summary from results, if available
  summary <- list()
  if (is.data.frame(analysis_tbl) && nrow(analysis_tbl) > 0) {
    if ("p_adj" %in% names(analysis_tbl)) {
      summary$min_fdr_p <- suppressWarnings(min(analysis_tbl$p_adj, na.rm = TRUE))
    }
    if ("es_value" %in% names(analysis_tbl)) {
      ord <- order(abs(analysis_tbl$es_value), decreasing = TRUE)
      top <- head(analysis_tbl[ord, c(intersect(c("var", "es_value", "p_adj"), names(analysis_tbl)))], 5)
      summary$top5_by_abs_es <- top
    }
    ft <- analysis_tbl$formula_text %||% NULL
    if (!is.null(ft) && length(ft) > 0) {
      f1 <- as.character(ft[1])
      summary$formula_template <- sub("^\\s*[^~]+~", "<y> ~", f1)
    } else {
      summary$formula_template <- NA_character_
    }
    summary$ref_level <- unique(analysis_tbl$ref_level %||% NA_character_)
    summary$target_level <- unique(analysis_tbl$target_level %||% NA_character_)
    summary$es_method <- unique(analysis_tbl$es_method %||% NA_character_)
    summary$p_adjust_method <- unique(analysis_tbl$p_adj_method %||% NA_character_)
  }

  list(
    timestamp = as.character(Sys.time()),
    package = list(
      name = "roiflow",
      version = .pkg_version("roiflow"),
      git_commit = .git_commit()
    ),
    io = list(
      input = input,
      input_md5 = .file_md5(input),
      out_root = out_root
    ),
    resolved = cfg$resolved %||% list(),
    specs = .sanitize_for_json(ctx$specs %||% list()),
    counts = list(
      n_rows_input = if (is.data.frame(df_in)) nrow(df_in) else NA_integer_,
      n_rows_final = if (is.data.frame(df_final)) nrow(df_final) else NA_integer_,
      roi_count = length(cfg$resolved$roi_cols %||% character(0)),
      roi_regex = cfg$roi$regex %||% NA_character_
    ),
    group = list(
      group_col = cfg$resolved$group_col %||% NA_character_,
      counts = ctx$meta$group_counts %||% NULL
    ),
    batch = list(
      site_col = cfg$columns$site %||% NA_character_,
      levels = ctx$meta$batch_levels %||% NULL
    ),
    timings = ctx$timings %||% list(),
    log = ctx$log %||% list(),
    results_summary = summary,
    session = list(
      r_version = paste0(R.version$major, ".", R.version$minor),
      platform = R.version$platform
    )
  )
}

#' Write a Manifest to JSON
#'
#' @param manifest List from \code{\link{build_manifest}}.
#' @param path Output file path.
#' @param overwrite Logical. If \code{FALSE}, error if file exists.
#'
#' @return Invisibly returns the output path.
#' @export
write_manifest <- function(manifest, path, overwrite = FALSE) {
  path <- .as_scalar_chr(path, "path")
  overwrite <- isTRUE(overwrite)
  if (!overwrite && file.exists(path)) {
    stop(sprintf("E_EXISTS: Manifest exists (set overwrite=TRUE): %s", path), call. = FALSE)
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  .check_pkg("jsonlite", "writing manifest")
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}
