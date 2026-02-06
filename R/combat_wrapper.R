################################################################################
# INTEGRATED COMBAT WRAPPERS (LEVEL-SAFE VERSION)
################################################################################

#' Define ComBat Wrapper Specifications
#'
#' @description
#' Creates a specification list for the high-level ComBat data frame wrappers.
#' These settings control ROI selection, batch column, covariate formula, and
#' safety behaviors for missingness and single-batch data.
#'
#' @param roi_cols Character vector. Explicit ROI column names. If \code{NULL},
#'   columns are inferred via \code{roi_regex}.
#' @param roi_regex Character. Regex to identify ROI columns when \code{roi_cols}
#'   is \code{NULL}. Default matches L_/R_ prefixed columns.
#' @param batch_col Character. Batch/site column name.
#' @param mod_formula Model formula for covariates (e.g., \code{~ Age + Sex}).
#'   Use \code{NULL} for no covariates.
#' @param eb Logical. Use empirical Bayes shrinkage in ComBat.
#' @param missing_fit Character. How to handle missing values during fit:
#'   \code{"impute"} (default) or \code{"error"} to stop if any NA is present
#'   in ROI data.
#' @param single_batch Character. If only one batch level is present:
#'   \code{"error"} (default) or \code{"noop"} to return a pass-through model.
#' @param keep_inputs Logical. If \code{TRUE}, return all columns on apply;
#'   otherwise return only ROI columns.
#'
#' @return A list of specification values consumed by the wrapper functions.
#' @export
combat_spec <- function(
  roi_cols      = NULL,                # explicit ROI columns
  roi_regex     = "^(L_|R_)",          # infer ROI columns if roi_cols is NULL
  batch_col     = "Site",
  mod_formula   = NULL,                # e.g., ~ Age + Sex
  eb            = TRUE,
  missing_fit   = c("impute", "error"),# if any NA present during fit
  single_batch  = c("error", "noop"),  # if only 1 batch level
  keep_inputs   = TRUE                 # keep non-ROI columns in output
) {
  list(
    roi_cols = roi_cols,
    roi_regex = roi_regex,
    batch_col = batch_col,
    mod_formula = mod_formula,
    eb = eb,
    missing_fit = match.arg(missing_fit),
    single_batch = match.arg(single_batch),
    keep_inputs = keep_inputs
  )
}

# -------- Internal Helpers --------

.append_log <- function(log, code, message, context = list()) {
  log[[length(log) + 1]] <- list(code = code, message = message, context = context)
  log
}

.infer_roi_cols <- function(df, spec) {
  if (!is.null(spec$roi_cols)) return(intersect(spec$roi_cols, names(df)))
  grep(spec$roi_regex, names(df), value = TRUE)
}

# Internal: level-safe model matrix builder.
# Ensures apply-time data has the exact same factor levels as fit-time data.
.build_mod <- function(df, spec, fit_levels = NULL) {
  if (is.null(spec$mod_formula)) return(NULL)
  
  # Sync factor levels to prevent design matrix column mismatch
  if (!is.null(fit_levels)) {
    for (nm in names(fit_levels)) {
      if (nm %in% names(df) && !is.null(fit_levels[[nm]])) {
        x0 <- as.character(df[[nm]])
        x1 <- factor(x0, levels = fit_levels[[nm]])
        if (any(!is.na(x0) & is.na(x1))) {
          stop(sprintf("New level(s) in %s at apply()", nm))
        }
        df[[nm]] <- x1
      }
    }
  }

  # Use model.frame to strictly pull from df environment
  mf <- stats::model.frame(spec$mod_formula, data = df, na.action = na.pass)
  mm <- stats::model.matrix(spec$mod_formula, data = mf)

  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  if (anyNA(mm)) {
    stop("Design matrix contains NA values (check covariate levels and missingness).")
  }
  mm
}

# -------- Main Functions --------

#' Fit ComBat on a Data Frame
#'
#' @description
#' Fits a ComBat model on a data frame containing ROI columns, a batch column,
#' and optional covariates. This wrapper handles ROI detection, model-matrix
#' construction with level safety, and returns a structured fit object.
#'
#' @param df Data frame with ROI columns and batch column.
#' @param spec Specification list from \code{\link{combat_spec}}.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class \code{combat_fit_obj} with elements:
#' \itemize{
#'   \item \code{model}: fitted model and metadata
#'   \item \code{qc}: basic QC metrics (rows, batches, ROI count)
#'   \item \code{log}: warning/info log entries
#'   \item \code{spec}: the specification used
#' }
#' @export
combat_fit_df <- function(df, spec = combat_spec(), verbose = interactive()) {
  stopifnot(is.data.frame(df))
  log <- list()
  qc  <- list()

  if (!spec$batch_col %in% names(df)) stop(paste("batch_col not found:", spec$batch_col))
  batch <- factor(as.character(df[[spec$batch_col]]))

  qc$n_rows_in    <- nrow(df)
  qc$n_batches    <- nlevels(batch)
  qc$batch_counts <- as.integer(table(batch))
  names(qc$batch_counts) <- levels(batch)

  # Handle Single Batch Case
  if (qc$n_batches < 2) {
    msg <- "ComBat requires >= 2 batches."
    if (spec$single_batch == "noop") {
      log <- .append_log(log, "W_SINGLE_BATCH_NOOP", msg, list(batch_col = spec$batch_col))
      return(structure(list(model = list(noop = TRUE, levels_batch = levels(batch)),
                            qc = qc, log = log, spec = spec), class = "combat_fit_obj"))
    }
    stop(msg)
  }

  roi_cols <- .infer_roi_cols(df, spec)
  if (length(roi_cols) == 0) stop("No ROI columns found.")
  qc$n_roi_cols <- length(roi_cols)

  # Prepare Data: Subjects as ROWS for this wrapper, 
  # but combat_fit expects ROIs as rows (usually), so we transpose for the internal call.
  dat <- as.matrix(df[, roi_cols, drop = FALSE])
  
  # Build Model and Capture Factor Levels for 'Apply' safety
  mod <- .build_mod(df, spec)
  if (!is.null(spec$mod_formula)) {
    mod_vars <- all.vars(spec$mod_formula)
    fit_levels <- lapply(df[, intersect(mod_vars, names(df)), drop=FALSE], function(x) {
      if (is.factor(x) || is.character(x)) levels(as.factor(x)) else NULL
    })
    fit_levels <- fit_levels[!vapply(fit_levels, is.null, logical(1))]
  } else {
    fit_levels <- list()
  }

  # Call patched internal fit (subjects are columns for the math)
  if (spec$missing_fit == "error" && anyNA(dat)) {
    stop("Missing values present in ROI data during fit and missing_fit = \"error\".")
  }
  tmp <- combat_fit(dat = t(dat), batch = batch, mod = mod, eb = spec$eb, verbose = verbose)

  model <- list(
    tmp = tmp,
    roi_cols = roi_cols,
    batch_col = spec$batch_col,
    mod_formula = spec$mod_formula,
    fit_levels = fit_levels,
    noop = FALSE
  )

  structure(list(model = model, qc = qc, log = log, spec = spec), class = "combat_fit_obj")
}

#' Apply ComBat to a Data Frame
#'
#' @description
#' Applies a fitted ComBat model to a new data frame. Validates batch levels and
#' ROI columns, aligns factor levels to the fit-time design, and returns a
#' harmonized data frame.
#'
#' @param fit Object returned by \code{\link{combat_fit_df}}.
#' @param df Data frame to harmonize.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class \code{combat_result} containing:
#' \itemize{
#'   \item \code{data}: harmonized data frame
#'   \item \code{qc}, \code{log}, \code{spec}, \code{model}: from fit
#' }
#' @export
combat_apply_df <- function(fit, df, verbose = interactive()) {
  stopifnot(inherits(fit, "combat_fit_obj"), is.data.frame(df))

  spec  <- fit$spec
  model <- fit$model
  
  if (isTRUE(model$noop)) {
    out <- if (spec$keep_inputs) df else df[, .infer_roi_cols(df, spec), drop = FALSE]
    return(structure(list(data = out, qc = fit$qc, log = fit$log, spec = spec, model = model), 
                     class = "combat_result"))
  }

  # Validate Batch
  batch_raw <- as.character(df[[model$batch_col]])
  batch <- factor(batch_raw, levels = model$tmp$levels_batch)
  if (any(is.na(batch))) stop("New/unknown batch levels in apply().")

  # Validate ROIs
  roi_cols <- intersect(model$roi_cols, names(df))
  if (length(roi_cols) != length(model$roi_cols)) stop("ROI column mismatch between Fit and Apply.")

  # Prepare Data (Transposed for internal ComBat math)
  dat <- as.matrix(df[, roi_cols, drop = FALSE])
  mod <- .build_mod(df, spec, fit_levels = model$fit_levels)

  # Internal Apply (t(dat) makes subjects columns)
  applied <- combat_apply(tmp = model$tmp, dat = t(dat), batch = batch, mod = mod, verbose = verbose)
  
  # Transpose back to subject-rows
  dat_h <- t(applied$dat.combat)

  df_out <- if (spec$keep_inputs) df else df[, roi_cols, drop = FALSE]
  df_out[, roi_cols] <- dat_h

  structure(list(data = df_out, qc = fit$qc, log = fit$log, spec = spec, model = model), 
            class = "combat_result")
}

#' One-Liner for Fit + Apply
#'
#' @description
#' Convenience helper that fits and immediately applies ComBat on the same data.
#'
#' @param df Data frame with ROI columns and batch column.
#' @param spec Specification list from \code{\link{combat_spec}}.
#' @param return Character. \code{"data"} returns only the harmonized data frame;
#'   \code{"result"} returns the full \code{combat_result} object.
#' @param verbose Logical. Print progress messages.
#'
#' @return Harmonized data frame or a \code{combat_result} object.
#' @export
combat_one_liner <- function(df, spec = combat_spec(), return = c("data", "result"), verbose = interactive()) {
  return <- match.arg(return)
  fit <- combat_fit_df(df, spec = spec, verbose = verbose)
  res <- combat_apply_df(fit, df, verbose = verbose)
  if (return == "data") return(res$data)
  res
}

#' Print Method for \code{combat_result}
#'
#' @param x Object of class \code{combat_result}.
#' @param ... Unused.
#' @export
print.combat_result <- function(x, ...) {
  cat("<combat_result>\n")
  cat(sprintf("Rows: %d | Batches: %d | ROI cols: %d\n", 
              nrow(x$data), x$qc$n_batches, length(x$model$roi_cols)))
  invisible(x)
}
