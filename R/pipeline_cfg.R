# Config-first pipeline helpers (JSON -> validated/resolved config)

.cfg_stop <- function(code, message) {
  stop(sprintf("%s: %s", code, message), call. = FALSE)
}

.as_scalar_chr <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1) {
    .cfg_stop("E_CFG_TYPE", sprintf("%s must be a single string.", name))
  }
  x
}

.as_scalar_lgl <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.logical(x) || length(x) != 1) {
    .cfg_stop("E_CFG_TYPE", sprintf("%s must be TRUE/FALSE.", name))
  }
  x
}

.as_chr_vec <- function(x, name) {
  if (is.null(x)) return(character(0))
  if (!is.character(x)) .cfg_stop("E_CFG_TYPE", sprintf("%s must be a string array.", name))
  x
}

.default_cfg <- function() {
  list(
    project = list(name = "two_group", version = NA_character_),
    io = list(
      input = file.path("data", "sample.csv"),
      out_root = file.path("outputs", "two_group_com"),
      overwrite = FALSE,
      save_data = TRUE,
      save_tables = TRUE,
      save_plots = TRUE,
      plot_formats = c("png")
    ),
    columns = list(
      site = "Site",
      group = "auto",
      age = "Age",
      sex = "Sex",
      icv = "ICV"
    ),
    roi = list(
      regex = "^(L_|R_)",
      drop_suffix = "_thickavg"
    ),
    pipeline = list(
      list(step = "prep", params = list()),
      list(step = "combat", params = list()),
      list(step = "compare_groups", params = list()),
      list(step = "plots", params = list()),
      list(step = "export", params = list())
    ),
    manifest = list(enabled = TRUE, path = file.path("logs", "run_manifest.json"))
  )
}

#' Read a Pipeline Config (JSON)
#'
#' @description
#' Reads a JSON config file used by \code{\link{run_pipeline}}.
#'
#' @param path Character. Path to a JSON file.
#'
#' @return A nested list.
#' @export
read_cfg <- function(path) {
  path <- .as_scalar_chr(path, "path")
  if (!file.exists(path)) .cfg_stop("E_CFG_NOT_FOUND", sprintf("Config file not found: %s", path))
  .check_pkg("jsonlite", "reading JSON config")
  jsonlite::read_json(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
}

.cfg_merge <- function(x, y) {
  # shallow-ish recursive merge: lists are merged by name; scalars overwrite
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  if (!is.list(x) || !is.list(y)) return(y)
  nms <- names(y)
  # JSON arrays often arrive as unnamed lists: treat as replacement.
  if (is.null(nms) || all(!nzchar(nms))) return(y)
  out <- x
  for (nm in nms) {
    out[[nm]] <- .cfg_merge(out[[nm]], y[[nm]])
  }
  out
}

.validate_step <- function(step, i) {
  if (!is.list(step)) .cfg_stop("E_CFG_TYPE", sprintf("pipeline[%d] must be an object.", i))
  if (is.null(step$step)) .cfg_stop("E_CFG_MISSING", sprintf("pipeline[%d].step is required.", i))
  if (!is.character(step$step) || length(step$step) != 1) {
    .cfg_stop("E_CFG_TYPE", sprintf("pipeline[%d].step must be a string.", i))
  }
  allowed <- c("prep", "combat", "compare_groups", "pairwise_groups", "assoc_scan", "corr_scan", "plots", "export")
  if (!step$step %in% allowed) {
    .cfg_stop("E_CFG_STEP", sprintf("Unsupported step '%s'. Allowed: %s",
                                    step$step, paste(allowed, collapse = ", ")))
  }
  if (is.null(step$params)) step$params <- list()
  if (!is.list(step$params)) .cfg_stop("E_CFG_TYPE", sprintf("pipeline[%d].params must be an object.", i))
  step
}

#' Validate and Normalize a Pipeline Config
#'
#' @description
#' Validates a config list (typically from \code{\link{read_cfg}}) and fills
#' missing fields with defaults. This function does not read input data; use
#' \code{\link{resolve_cfg}} for data-dependent auto-detection.
#'
#' @param cfg A nested list.
#'
#' @return A normalized config list.
#' @export
validate_cfg <- function(cfg) {
  if (!is.list(cfg)) .cfg_stop("E_CFG_TYPE", "cfg must be a list/object.")
  cfg <- .cfg_merge(.default_cfg(), cfg)

  # io
  cfg$io$input <- .as_scalar_chr(cfg$io$input, "io.input")
  cfg$io$out_root <- .as_scalar_chr(cfg$io$out_root, "io.out_root")
  cfg$io$overwrite <- .as_scalar_lgl(cfg$io$overwrite, "io.overwrite") %||% FALSE
  cfg$io$save_data <- .as_scalar_lgl(cfg$io$save_data, "io.save_data") %||% TRUE
  cfg$io$save_tables <- .as_scalar_lgl(cfg$io$save_tables, "io.save_tables") %||% TRUE
  cfg$io$save_plots <- .as_scalar_lgl(cfg$io$save_plots, "io.save_plots") %||% TRUE
  cfg$io$plot_formats <- .as_chr_vec(cfg$io$plot_formats, "io.plot_formats")
  if (length(cfg$io$plot_formats) == 0) cfg$io$plot_formats <- c("png")

  # columns
  cfg$columns$site <- .as_scalar_chr(cfg$columns$site, "columns.site")
  cfg$columns$group <- .as_scalar_chr(cfg$columns$group, "columns.group")
  cfg$columns$age <- .as_scalar_chr(cfg$columns$age, "columns.age")
  cfg$columns$sex <- .as_scalar_chr(cfg$columns$sex, "columns.sex")
  cfg$columns$icv <- .as_scalar_chr(cfg$columns$icv, "columns.icv")

  # roi
  cfg$roi$regex <- .as_scalar_chr(cfg$roi$regex, "roi.regex")
  cfg$roi$drop_suffix <- .as_scalar_chr(cfg$roi$drop_suffix, "roi.drop_suffix") %||% ""

  # pipeline steps
  if (is.null(cfg$pipeline) || !is.list(cfg$pipeline) || length(cfg$pipeline) == 0) {
    .cfg_stop("E_CFG_MISSING", "pipeline must be a non-empty array of steps.")
  }
  cfg$pipeline <- lapply(seq_along(cfg$pipeline), function(i) .validate_step(cfg$pipeline[[i]], i))

  # manifest
  cfg$manifest$enabled <- .as_scalar_lgl(cfg$manifest$enabled, "manifest.enabled") %||% TRUE
  cfg$manifest$path <- .as_scalar_chr(cfg$manifest$path, "manifest.path")

  cfg
}

.detect_group_col <- function(names_vec) {
  candidates <- c("Group", "group", "groups", "CaseControl", "Diagnosis")
  hit <- candidates[candidates %in% names_vec]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

#' Resolve a Config Against a Data Frame
#'
#' @description
#' Performs data-dependent resolution (e.g., auto-detect group column, derive ROI
#' columns). Used by \code{\link{run_pipeline}}.
#'
#' @param cfg Config list (validated).
#' @param df Optional data frame. If provided, enables auto-detection.
#'
#' @return Resolved config list with \code{resolved} fields filled.
#' @export
resolve_cfg <- function(cfg, df = NULL) {
  cfg <- validate_cfg(cfg)

  cfg$resolved <- cfg$resolved %||% list()

  if (!is.null(df)) {
    stopifnot(is.data.frame(df))

    group_col <- cfg$columns$group
    if (identical(group_col, "auto")) {
      group_col <- .detect_group_col(names(df))
      if (is.null(group_col)) {
        .cfg_stop("E_GROUP_DETECT",
                  "Could not auto-detect group column. Set columns.group explicitly.")
      }
    }
    cfg$resolved$group_col <- group_col

    roi_cols <- cfg$roi$cols %||% NULL
    if (is.null(roi_cols)) {
      roi_cols <- grep(cfg$roi$regex, names(df), value = TRUE)
    }
    if (length(roi_cols) == 0) {
      .cfg_stop("E_ROI_EMPTY", sprintf("No ROI columns found for roi.regex='%s'.", cfg$roi$regex))
    }
    cfg$resolved$roi_cols <- roi_cols
  }

  cfg
}
