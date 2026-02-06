# Utilities for profiling analyses (logging, checks, helpers)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.log_add <- function(log, code, message, context = list()) {
  log[[length(log) + 1]] <- list(code = code, message = message, context = context)
  log
}

.check_pkg <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- if (is.null(reason)) {
      sprintf("Package '%s' is required but not installed.", pkg)
    } else {
      sprintf("Package '%s' is required for %s but not installed.", pkg, reason)
    }
    stop(msg, call. = FALSE)
  }
  TRUE
}

.safe_factor_levels <- function(x) {
  if (is.factor(x)) return(levels(x))
  if (is.character(x)) return(levels(factor(x)))
  NULL
}

.set_ref_factor <- function(x, ref, var_name = "factor") {
  if (is.null(ref) || is.na(ref) || !nzchar(as.character(ref))) {
    stop(sprintf("E_REF_NOT_SET: Reference level for %s is not set.", var_name), call. = FALSE)
  }
  f <- factor(as.character(x))
  if (!ref %in% levels(f)) {
    stop(sprintf("E_REF_NOT_FOUND: Reference level '%s' not found in %s.", ref, var_name), call. = FALSE)
  }
  stats::relevel(f, ref = ref)
}

.infer_ref_level <- function(levels_vec, preferred = c("Control", "0")) {
  for (p in preferred) {
    if (p %in% levels_vec) return(p)
  }
  NULL
}

.apply_ref_levels <- function(df, ref_levels) {
  if (is.null(ref_levels) || length(ref_levels) == 0) return(df)
  for (nm in names(ref_levels)) {
    if (nm %in% names(df)) {
      df[[nm]] <- .set_ref_factor(df[[nm]], ref_levels[[nm]], nm)
    }
  }
  df
}

.quote_var <- function(x) {
  ifelse(grepl("[^A-Za-z0-9_.]", x), sprintf("`%s`", x), x)
}

.ensure_cols <- function(df, cols, context = "") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("E_MISSING_COL: Missing required columns%s: %s",
                 if (nzchar(context)) paste0(" in ", context) else "",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  TRUE
}

.capture_warnings <- function(expr) {
  warns <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warns)
}

.classify_warning <- function(msg) {
  if (grepl("converge|convergence", msg, ignore.case = TRUE)) return("W_CONVERGENCE")
  if (grepl("singular", msg, ignore.case = TRUE)) return("W_SINGULAR_FIT")
  if (grepl("separat", msg, ignore.case = TRUE)) return("W_SEPARATION")
  if (grepl("rank deficient", msg, ignore.case = TRUE)) return("W_RANK_DEFICIENT")
  "W_MODEL_WARNING"
}

.result_object <- function(data, diagnostics = NULL, summaries = NULL,
                           plots = NULL, meta = list(), log = list()) {
  structure(
    list(
      data = data,
      diagnostics = diagnostics,
      summaries = summaries,
      plots = plots,
      meta = meta,
      log = log
    ),
    class = "profile_result"
  )
}

.complete_cases <- function(df) {
  stats::complete.cases(df)
}
