# Spin-correlation I/O helpers.

#' Tidy Spin Correlation Results
#'
#' @description
#' Converts spin correlation outputs into a long table with one row per map pair.
#'
#' @param result A \code{spin_correlation_result} object, or a list containing
#' \code{observed} and \code{p_value} matrices.
#' @param upper_only Logical. If \code{TRUE}, returns only upper triangle when
#' the matrix is square and names match.
#' @param drop_diag Logical. Drop diagonal entries.
#'
#' @return A data frame with columns:
#' \code{x_map}, \code{y_map}, \code{observed_r}, \code{p_spin}.
#' @export
spin_results_table <- function(result, upper_only = FALSE, drop_diag = TRUE) {
  if (!is.list(result) || is.null(result$observed) || is.null(result$p_value)) {
    .spin_stop("E_SPIN_RESULT", "`result` must contain `observed` and `p_value`.")
  }
  obs <- result$observed
  psp <- result$p_value
  if (!is.matrix(obs) || !is.matrix(psp)) {
    .spin_stop("E_SPIN_RESULT", "`observed` and `p_value` must be matrices.")
  }
  if (!all(dim(obs) == dim(psp))) {
    .spin_stop("E_SPIN_RESULT", "Dimension mismatch between `observed` and `p_value`.")
  }

  rn <- rownames(obs) %||% paste0("x_", seq_len(nrow(obs)))
  cn <- colnames(obs) %||% paste0("y_", seq_len(ncol(obs)))
  rownames(obs) <- rn
  colnames(obs) <- cn
  rownames(psp) <- rn
  colnames(psp) <- cn

  idx <- expand.grid(i = seq_len(nrow(obs)), j = seq_len(ncol(obs)), KEEP.OUT.ATTRS = FALSE)
  out <- data.frame(
    x_map = rn[idx$i],
    y_map = cn[idx$j],
    observed_r = as.numeric(obs[cbind(idx$i, idx$j)]),
    p_spin = as.numeric(psp[cbind(idx$i, idx$j)]),
    stringsAsFactors = FALSE
  )

  is_square <- nrow(obs) == ncol(obs) && identical(rn, cn)
  if (is_square && isTRUE(upper_only)) {
    keep <- idx$i <= idx$j
    out <- out[keep, , drop = FALSE]
    idx <- idx[keep, , drop = FALSE]
  }
  if (is_square && isTRUE(drop_diag)) {
    keep <- idx$i != idx$j
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Save Spin Correlation Outputs
#'
#' @description
#' Optional helper to persist spin test outputs. Side effects are opt-in.
#'
#' @param result A \code{spin_correlation_result} object.
#' @param out_dir Output directory.
#' @param prefix File prefix.
#' @param formats Character vector: any of \code{"csv"}, \code{"rds"}.
#' @param include_null Logical. Include null arrays in a separate RDS.
#'
#' @return Invisible list of written files.
#' @export
save_spin_results <- function(result, out_dir, prefix = "spin",
                              formats = c("csv", "rds"),
                              include_null = FALSE) {
  if (!inherits(result, "spin_correlation_result")) {
    .spin_stop("E_SPIN_RESULT", "`result` must be a `spin_correlation_result`.")
  }
  if (missing(out_dir) || is.null(out_dir) || !nzchar(out_dir)) {
    .spin_stop("E_SPIN_IO", "`out_dir` is required.")
  }

  formats <- unique(as.character(formats))
  ok <- c("csv", "rds")
  bad <- setdiff(formats, ok)
  if (length(bad) > 0) {
    .spin_stop("E_SPIN_IO", sprintf("Unsupported format(s): %s.", paste(bad, collapse = ", ")))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list()
  tbl <- spin_results_table(result, upper_only = FALSE, drop_diag = FALSE)

  if ("csv" %in% formats) {
    f_csv <- file.path(out_dir, paste0(prefix, "_summary.csv"))
    utils::write.csv(tbl, f_csv, row.names = FALSE)
    files$csv <- f_csv
  }
  if ("rds" %in% formats) {
    f_rds <- file.path(out_dir, paste0(prefix, "_result.rds"))
    saveRDS(result, f_rds)
    files$rds <- f_rds
  }
  if (isTRUE(include_null) && !is.null(result$null_xy) && !is.null(result$null_yx)) {
    f_null <- file.path(out_dir, paste0(prefix, "_null.rds"))
    saveRDS(list(null_xy = result$null_xy, null_yx = result$null_yx), f_null)
    files$null <- f_null
  }
  invisible(files)
}
