#' Create Preprocessing Specification
#'
#' @description
#' Creates a configuration object that controls how \code{\link{prep}} preprocesses
#' ROI data. This allows you to customize column detection patterns, factor levels,
#' outlier handling, and missing data policies.
#'
#' @param roi_regex Character. Regular expression to identify ROI columns.
#'   Default: \code{"^(L_|R_)"} matches columns starting with L_ or R_.
#' @param numeric_regex Character. Regular expression to identify columns that
#'   should be converted to numeric. Default matches ROI columns, Age, ICV, and mod_PDS.
#' @param zero_to_na_regex Character. Regular expression to identify columns where
#'   zeros should be converted to NA. Default matches ROI columns and ICV.
#' @param factor_levels Named list. Specifies factor variables and their reference
#'   levels. Default: \code{list(Sex = "0", Group = "0")}.
#' @param site_col Character. Name of the site/scanner column. Default: \code{"Site"}.
#' @param site_ref Character or vector. How to set the site reference level:
#'   \itemize{
#'     \item \code{"largest"}: Use the site with most observations (default)
#'     \item \code{"keep"}: Keep existing factor order
#'     \item A specific site name: Use that site as reference
#'   }
#' @param icv_col Character. Name of the intracranial volume column. Default: \code{"ICV"}.
#' @param icv_scale Numeric. Factor to divide ICV by. Default: \code{1e6} (converts to millions).
#' @param na_action Character. How to handle missing data:
#'   \itemize{
#'     \item \code{"drop_required"}: Drop rows with NA in required columns (default)
#'     \item \code{"drop_all"}: Drop rows with any NA
#'     \item \code{"keep"}: Keep all rows
#'   }
#' @param required_cols Character vector. Columns that must be complete for
#'   \code{na_action = "drop_required"}. If NULL (default), inferred as key covariates
#'   (Site, Sex, Group, Age, ICV if present).
#' @param outlier_action Character. How to handle outliers:
#'   \itemize{
#'     \item \code{"clip"}: Clip outliers to fence boundaries (default)
#'     \item \code{"flag"}: Detect but don't modify
#'     \item \code{"none"}: Skip outlier detection
#'   }
#' @param outlier_method Character. Method for outlier detection:
#'   \itemize{
#'     \item \code{"iqr"}: Interquartile range method (default)
#'     \item \code{"sd"}: Standard deviation method
#'   }
#' @param outlier_k Numeric. Multiplier for outlier fences. Default: 3.
#'   For IQR: fences at Q1 - k*IQR and Q3 + k*IQR.
#'   For SD: fences at mean ± k*SD.
#'
#' @return A list with class attributes containing all preprocessing specifications.
#'
#' @details
#' The specification object controls all aspects of the preprocessing pipeline.
#' Common customizations:
#' \itemize{
#'   \item Change ROI pattern: \code{roi_regex = "^ROI_"}
#'   \item More conservative outliers: \code{outlier_k = 2.5}
#'   \item Keep all data: \code{na_action = "keep", outlier_action = "none"}
#'   \item Custom factor levels: \code{factor_levels = list(Sex = "M", Diagnosis = "HC")}
#' }
#'
#' @examples
#' # Default specification
#' spec <- prep_spec()
#'
#' # Custom ROI pattern and stricter outlier detection
#' spec <- prep_spec(
#'   roi_regex = "^ROI_",
#'   outlier_k = 2.5,
#'   outlier_method = "sd"
#' )
#'
#' # Keep all data, no outlier clipping
#' spec <- prep_spec(
#'   na_action = "keep",
#'   outlier_action = "flag"
#' )
#'
#' # Use specific site as reference
#' spec <- prep_spec(site_ref = "SITE_A")
#'
#' @seealso \code{\link{prep}} for the main preprocessing function
#' @export
prep_spec <- function(
  roi_regex        = "^(L_|R_)",
  numeric_regex    = "^(L_|R_)|Age$|ICV$|mod_PDS$",
  zero_to_na_regex = "^(L_|R_)|^ICV$",
  factor_levels    = list(
    Sex   = "0",
    Group = "0"
  ),
  site_col         = "Site",
  site_ref         = c("largest", "keep"),  # or a specific level like "SITE_A"
  icv_col          = "ICV",
  icv_scale        = 1e6,
  na_action        = c("drop_required", "drop_all", "keep"),
  required_cols    = NULL,                  # if NULL, inferred from key covariates
  outlier_action   = c("clip", "flag", "none"),
  outlier_method   = c("iqr", "sd"),
  outlier_k        = 3
) {
  na_action <- match.arg(na_action)
  outlier_action <- match.arg(outlier_action)
  outlier_method <- match.arg(outlier_method)

  # normalize site_ref without blocking explicit level
  site_ref_norm <- site_ref
  if (length(site_ref_norm) == 1 && site_ref_norm %in% c("largest", "keep")) {
    # ok
  } else if (length(site_ref_norm) == 1) {
    site_ref_norm <- as.character(site_ref_norm)
  } else {
    site_ref_norm <- "largest"
  }
  site_ref <- site_ref_norm

  list(
    roi_regex = roi_regex,
    numeric_regex = numeric_regex,
    zero_to_na_regex = zero_to_na_regex,
    factor_levels = factor_levels,
    site_col = site_col,
    site_ref = site_ref,
    icv_col = icv_col,
    icv_scale = icv_scale,
    na_action = na_action,
    required_cols = required_cols,
    outlier_action = outlier_action,
    outlier_method = outlier_method,
    outlier_k = outlier_k
  )
}

# ---------- Helpers ----------

#' Safely Convert to Numeric
#'
#' @description
#' Internal helper that converts values to numeric while tracking how many
#' new NAs are introduced by coercion.
#'
#' @param x Vector to convert to numeric.
#'
#' @return List with two elements:
#'   \itemize{
#'     \item \code{value}: Numeric vector
#'     \item \code{n_na_new}: Count of new NAs introduced by coercion
#'   }
#'
#' @keywords internal
#' @noRd
.safe_as_numeric <- function(x) {
  # preserve existing NAs, warn only via returned stats
  x_chr <- as.character(x)
  if (length(x_chr) > 0) {
    x_chr[trimws(x_chr) == ""] <- NA_character_
  }
  x_num <- suppressWarnings(as.numeric(x_chr))
  list(value = x_num, n_na_new = sum(is.na(x_num) & !is.na(x_chr)))
}

#' Append Log Entry
#'
#' @description
#' Internal helper to append structured log entries during preprocessing.
#'
#' @param log List of existing log entries.
#' @param code Character. Log code (e.g., "W_NO_ROI_COLS", "I_ICV_RESCALE").
#' @param message Character. Human-readable log message.
#' @param context List. Additional context information.
#'
#' @return Updated log list.
#'
#' @keywords internal
#' @noRd
._append_log <- function(log, code, message, context = list()) {
  log[[length(log) + 1]] <- list(code = code, message = message, context = context)
  log
}

# ---------- Main ----------

#' Preprocess ROI Data
#'
#' @description
#' Main preprocessing function for ROI-based neuroimaging data. Performs a
#' comprehensive, configurable pipeline including type conversion, rescaling,
#' factor handling, missing data management, and outlier detection/clipping.
#'
#' This function is designed to be both human and LLM-agent friendly with:
#' \itemize{
#'   \item Flexible configuration via \code{\link{prep_spec}}
#'   \item Detailed logging of all operations
#'   \item Quality control metrics
#'   \item Informative error messages
#'   \item Option to return full result object or just cleaned data
#' }
#'
#' @param df Data frame containing ROI measurements and covariates.
#' @param spec Preprocessing specification created by \code{\link{prep_spec}}.
#'   Default uses standard settings for L_/R_ prefixed ROI columns.
#' @param return Character. What to return:
#'   \itemize{
#'     \item \code{"data"}: Return only the cleaned data frame (default)
#'     \item \code{"result"}: Return full result object with data, QC metrics, and logs
#'   }
#' @param verbose Logical. Print progress summary. Default: \code{TRUE} in interactive
#'   sessions, \code{FALSE} otherwise.
#' @param strict Logical. If \code{TRUE}, stop on warnings (e.g., missing reference
#'   levels, no ROI columns). Default: \code{FALSE}.
#'
#' @return
#' If \code{return = "data"}: A cleaned data frame.
#'
#' If \code{return = "result"}: An object of class \code{roiflow_prep_result} containing:
#' \itemize{
#'   \item \code{data}: Cleaned data frame
#'   \item \code{qc}: List of quality control metrics:
#'     \itemize{
#'       \item \code{n_rows_in}, \code{n_rows_out}: Input/output row counts
#'       \item \code{n_cols_in}, \code{n_cols_out}: Input/output column counts
#'       \item \code{n_roi_cols}: Number of ROI columns detected
#'       \item \code{n_rows_dropped_na}: Rows removed due to missing data
#'       \item \code{n_zero_to_na}: Number of zeros converted to NA
#'       \item \code{n_outlier_points}: Total outlier values detected/clipped
#'       \item \code{outliers_per_roi}: Named vector of outlier counts per ROI
#'       \item \code{outlier_index}: List of row indices per ROI (after NA filtering)
#'       \item \code{required_cols}: Columns used for completeness check
#'     }
#'   \item \code{log}: List of log entries with codes, messages, and context
#'   \item \code{spec}: The preprocessing specification used
#' }
#'
#' @details
#' The preprocessing pipeline executes these steps in order:
#'
#' \strong{1. Column Identification}
#' \itemize{
#'   \item Identifies ROI columns using \code{roi_regex}
#'   \item Identifies numeric columns using \code{numeric_regex}
#'   \item Identifies zero-to-NA columns using \code{zero_to_na_regex}
#' }
#'
#' \strong{2. Type Conversion}
#' \itemize{
#'   \item Converts matched columns to numeric
#'   \item Tracks new NAs introduced by coercion
#'   \item Logs warnings for problematic conversions
#' }
#'
#' \strong{3. ICV Rescaling}
#' \itemize{
#'   \item Divides ICV by \code{icv_scale} (default: 1,000,000)
#'   \item Useful for numerical stability in models
#' }
#'
#' \strong{4. Factor Handling}
#' \itemize{
#'   \item Converts specified columns to factors
#'   \item Sets reference levels (e.g., Sex = "0", Group = "0")
#'   \item For Site: optionally sets largest site as reference
#'   \item Warns if specified reference levels don't exist
#' }
#'
#' \strong{5. Zero to NA Conversion}
#' \itemize{
#'   \item Converts zeros to NA in ROI and ICV columns
#'   \item Common in neuroimaging when 0 indicates missing/failed measurement
#' }
#'
#' \strong{6. Missing Data Handling}
#' \itemize{
#'   \item \code{na_action = "drop_required"}: Removes rows with NA in required columns
#'   \item \code{na_action = "drop_all"}: Removes rows with any NA
#'   \item \code{na_action = "keep"}: Keeps all rows
#'   \item Required columns auto-detected from key covariates if not specified
#' }
#'
#' \strong{7. Outlier Detection/Clipping}
#' \itemize{
#'   \item IQR method: Outliers beyond Q1 - k*IQR or Q3 + k*IQR
#'   \item SD method: Outliers beyond mean ± k*SD
#'   \item \code{outlier_action = "clip"}: Clips to fence boundaries
#'   \item \code{outlier_action = "flag"}: Detects but doesn't modify
#'   \item Applied only to ROI columns
#' }
#'
#' @section Log Codes:
#' The log uses structured codes for programmatic parsing:
#' \itemize{
#'   \item \code{W_*}: Warnings (e.g., \code{W_NO_ROI_COLS}, \code{W_NUMERIC_COERCE_NA})
#'   \item \code{I_*}: Informational (e.g., \code{I_ICV_RESCALE}, \code{I_FACTOR_RELEVEL})
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - returns cleaned data
#' clean_data <- prep(raw_data)
#'
#' # Get full result with QC metrics and logs
#' result <- prep(raw_data, return = "result")
#' print(result)
#' summary(result$qc)
#' View(result$log)
#'
#' # Custom specification
#' spec <- prep_spec(
#'   outlier_k = 2.5,
#'   na_action = "drop_all",
#'   site_ref = "SITE_A"
#' )
#' clean_data <- prep(raw_data, spec = spec)
#'
#' # Strict mode - stop on warnings
#' clean_data <- prep(raw_data, strict = TRUE)
#'
#' # Silent mode
#' clean_data <- prep(raw_data, verbose = FALSE)
#' }
#'
#' @seealso
#' \code{\link{prep_spec}} for creating preprocessing specifications
#'
#' @export
#' @importFrom stats quantile complete.cases relevel sd
prep <- function(df,
                 spec = prep_spec(),
                 return = c("data", "result"),
                 verbose = interactive(),
                 strict = FALSE) {
  stopifnot(is.data.frame(df))
  return <- match.arg(return)

  log <- list()
  qc <- list()
  qc$outlier_index <- list()

  # ---- identify columns
  cols <- names(df)
  roi_cols <- grep(spec$roi_regex, cols, value = TRUE)
  num_cols <- grep(spec$numeric_regex, cols, value = TRUE)
  zero_cols <- grep(spec$zero_to_na_regex, cols, value = TRUE)

  qc$n_rows_in <- nrow(df)
  qc$n_cols_in <- ncol(df)
  qc$n_roi_cols <- length(roi_cols)

  if (length(roi_cols) == 0) {
    log <- ._append_log(log, "W_NO_ROI_COLS",
                        "No ROI columns detected (roi_regex matched nothing).",
                        list(roi_regex = spec$roi_regex))
    if (strict) stop("No ROI columns detected.")
  }

  # ---- coerce numeric
  n_na_new_total <- 0
  for (col in intersect(num_cols, cols)) {
    res <- .safe_as_numeric(df[[col]])
    df[[col]] <- res$value
    if (res$n_na_new > 0) {
      n_na_new_total <- n_na_new_total + res$n_na_new
      log <- ._append_log(
        log, "W_NUMERIC_COERCE_NA",
        sprintf("Coercion introduced %d new NA(s) in %s.", res$n_na_new, col),
        list(column = col, n_na_new = res$n_na_new)
      )
    }
  }
  qc$n_na_new_from_numeric <- n_na_new_total

  # ---- ICV rescale
  if (!is.null(spec$icv_col) && spec$icv_col %in% cols) {
    df[[spec$icv_col]] <- df[[spec$icv_col]] / spec$icv_scale
    log <- ._append_log(log, "I_ICV_RESCALE",
                        sprintf("Rescaled %s by /%g.", spec$icv_col, spec$icv_scale),
                        list(icv_col = spec$icv_col, icv_scale = spec$icv_scale))
  }

  # ---- factor handling (Sex/Group etc.)
  if (length(spec$factor_levels) > 0) {
    for (nm in names(spec$factor_levels)) {
      if (nm %in% cols) {
        ref <- spec$factor_levels[[nm]]
        df[[nm]] <- factor(as.character(df[[nm]]))
        if (ref %in% levels(df[[nm]])) {
          df[[nm]] <- stats::relevel(df[[nm]], ref = ref)
          log <- ._append_log(log, "I_FACTOR_RELEVEL",
                              sprintf("Releveled %s with reference %s.", nm, ref),
                              list(column = nm, ref = ref))
        } else {
          log <- ._append_log(log, "W_FACTOR_REF_MISSING",
                              sprintf("Reference level %s not found in %s.", ref, nm),
                              list(column = nm, ref = ref, levels = levels(df[[nm]])))
          if (strict) stop(sprintf("Reference level %s not found in %s.", ref, nm))
        }
      }
    }
  }

  # ---- Site factor + reference
  if (!is.null(spec$site_col) && spec$site_col %in% cols) {
    df[[spec$site_col]] <- factor(as.character(df[[spec$site_col]]))
    if (spec$site_ref == "largest") {
      tab <- sort(table(df[[spec$site_col]]), decreasing = TRUE)
      if (length(tab) == 0) {
        log <- ._append_log(log, "W_SITE_ALL_NA",
                            "Site column has no non-missing values.",
                            list(site_col = spec$site_col))
      } else {
        ref <- names(tab)[1]
        df[[spec$site_col]] <- stats::relevel(df[[spec$site_col]], ref = ref)
        log <- ._append_log(log, "I_SITE_REF_LARGEST",
                            sprintf("Set %s reference to largest site: %s.", spec$site_col, ref),
                            list(site_col = spec$site_col, ref = ref,
                                 counts = stats::setNames(as.integer(tab), names(tab))))
      }
    } else if (spec$site_ref != "keep") {
      if (spec$site_ref %in% levels(df[[spec$site_col]])) {
        df[[spec$site_col]] <- stats::relevel(df[[spec$site_col]], ref = spec$site_ref)
        log <- ._append_log(log, "I_SITE_REF_EXPLICIT",
                            sprintf("Set %s reference to %s.", spec$site_col, spec$site_ref),
                            list(site_col = spec$site_col, ref = spec$site_ref))
      } else {
        log <- ._append_log(log, "W_SITE_REF_MISSING",
                            "Explicit site reference not found in levels.",
                            list(site_col = spec$site_col, ref = spec$site_ref,
                                 levels = levels(df[[spec$site_col]])))
        if (strict) stop("Explicit site reference not found in levels.")
      }
    }
  }

  # ---- zeros to NA (ROI/ICV)
  n_zero_to_na <- 0
  for (col in intersect(zero_cols, cols)) {
    if (is.numeric(df[[col]])) {
      idx <- which(df[[col]] == 0)
      if (length(idx) > 0) {
        df[[col]][idx] <- NA_real_
        n_zero_to_na <- n_zero_to_na + length(idx)
        log <- ._append_log(log, "I_ZERO_TO_NA",
                            sprintf("Set %d zero(s) to NA in %s.", length(idx), col),
                            list(column = col, n = length(idx)))
      }
    }
  }
  qc$n_zero_to_na <- n_zero_to_na

  # ---- NA policy
  if (is.null(spec$required_cols)) {
    # infer a conservative required set: key covariates if present
    req <- intersect(c(spec$site_col, names(spec$factor_levels), "Age", spec$icv_col), cols)
    spec$required_cols <- req
  } else {
    spec$required_cols <- intersect(spec$required_cols, cols)
  }

  qc$required_cols <- spec$required_cols
  n_before_na <- nrow(df)

  if (spec$na_action == "drop_required" && length(spec$required_cols) > 0) {
    keep <- stats::complete.cases(df[, spec$required_cols, drop = FALSE])
    df <- df[keep, , drop = FALSE]
  } else if (spec$na_action == "drop_all") {
    df <- df[stats::complete.cases(df), , drop = FALSE]
  } else {
    # keep
  }

  qc$n_rows_after_na <- nrow(df)
  qc$n_rows_dropped_na <- n_before_na - nrow(df)

  # ---- outlier handling on ROI columns
  outlier_counts <- integer(0)
  if (spec$outlier_action != "none" && length(roi_cols) > 0) {
    roi_cols2 <- intersect(roi_cols, names(df))
    for (region in roi_cols2) {
      x <- df[[region]]
      if (!is.numeric(x)) next

      if (spec$outlier_method == "iqr") {
        q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
        q1 <- as.numeric(q[1]); q3 <- as.numeric(q[2])
        iqrk <- spec$outlier_k * (q3 - q1)
        low <- q1 - iqrk
        up  <- q3 + iqrk
      } else {
        mu <- mean(x, na.rm = TRUE)
        sdv <- stats::sd(x, na.rm = TRUE)
        low <- mu - spec$outlier_k * sdv
        up  <- mu + spec$outlier_k * sdv
      }

      idx_low <- which(x < low)
      idx_up  <- which(x > up)
      idx <- sort(unique(c(idx_low, idx_up)))
      if (length(idx) > 0) {
        outlier_counts[region] <- length(idx)
        qc$outlier_index[[region]] <- idx

        if (spec$outlier_action == "clip") {
          x[idx_low] <- low
          x[idx_up]  <- up
          df[[region]] <- x
        }
      }
    }
  }
  qc$outliers_per_roi <- outlier_counts
  qc$n_outlier_points <- sum(outlier_counts)

  qc$n_rows_out <- nrow(df)
  qc$n_cols_out <- ncol(df)

  if (verbose) {
    message(sprintf("prep(): %d -> %d rows; NA-dropped=%d; outliers=%d points",
                    qc$n_rows_in, qc$n_rows_out, qc$n_rows_dropped_na, qc$n_outlier_points))
  }

  result <- structure(
    list(data = df, qc = qc, log = log, spec = spec),
    class = "roiflow_prep_result"
  )

  if (return == "data") return(result$data)
  result
}

# ---------- Print Method ----------

#' Print Method for Preprocessing Results
#'
#' @description
#' Prints a summary of preprocessing results including row changes,
#' ROI column count, outlier statistics, and data quality metrics.
#'
#' @param x Object of class \code{roiflow_prep_result}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun{
#' result <- prep(raw_data, return = "result")
#' print(result)
#' # <roiflow_prep_result>
#' # Rows: 1000 -> 987 (dropped NA: 13)
#' # ROI cols: 68 | Outlier points: 142 | New NA from numeric: 0
#' }
#'
#' @export
print.roiflow_prep_result <- function(x, ...) {
  qc <- x$qc
  cat("<roiflow_prep_result>\n")
  cat(sprintf("Rows: %d -> %d (dropped NA: %d)\n",
              qc$n_rows_in, qc$n_rows_out, qc$n_rows_dropped_na))
  cat(sprintf("ROI cols: %d | Outlier points: %d | New NA from numeric: %d\n",
              qc$n_roi_cols, qc$n_outlier_points, qc$n_na_new_from_numeric))
  invisible(x)
}
