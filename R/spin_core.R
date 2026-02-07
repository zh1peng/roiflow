# Spin-correlation core utilities.

.spin_stop <- function(code, msg) {
  stop(sprintf("%s: %s", code, msg), call. = FALSE)
}

.as_numeric_matrix <- function(x, arg = "x") {
  if (is.vector(x) && !is.list(x)) {
    x <- matrix(as.numeric(x), ncol = 1)
    colnames(x) <- arg
    return(x)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    .spin_stop("E_SPIN_INPUT", sprintf("`%s` must be a numeric vector, matrix, or data.frame.", arg))
  }
  storage.mode(x) <- "double"
  if (!is.numeric(x)) {
    .spin_stop("E_SPIN_INPUT", sprintf("`%s` must be numeric.", arg))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0(arg, "_", seq_len(ncol(x)))
  }
  x
}

.load_default_centroids <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("lh_centroid_data", package = "roiflow", envir = env)
  utils::data("rh_centroid_data", package = "roiflow", envir = env)
  if (!exists("lh_centroid_data", envir = env, inherits = FALSE) ||
      !exists("rh_centroid_data", envir = env, inherits = FALSE)) {
    .spin_stop(
      "E_SPIN_CENTROID_MISSING",
      "Centroid data not found. Provide `lh_centroid_data` and `rh_centroid_data`."
    )
  }
  list(
    lh = get("lh_centroid_data", envir = env, inherits = FALSE),
    rh = get("rh_centroid_data", envir = env, inherits = FALSE)
  )
}

.extract_centroid_matrix <- function(df, hemi) {
  if (!is.data.frame(df)) {
    .spin_stop("E_SPIN_CENTROID", sprintf("`%s_centroid_data` must be a data.frame.", hemi))
  }
  if (all(c("centroid1", "centroid2", "centroid3") %in% names(df))) {
    out <- as.matrix(df[, c("centroid1", "centroid2", "centroid3"), drop = FALSE])
  } else {
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) < 3) {
      .spin_stop("E_SPIN_CENTROID", sprintf(
        "`%s_centroid_data` needs centroid1/2/3 (or at least three numeric columns).", hemi
      ))
    }
    out <- as.matrix(df[, num_cols[seq_len(3)], drop = FALSE])
  }
  if (ncol(out) != 3) {
    .spin_stop("E_SPIN_CENTROID", sprintf("`%s_centroid_data` must provide 3D coordinates.", hemi))
  }
  out
}

.euclid_dist_matrix <- function(a, b) {
  aa <- rowSums(a^2)
  bb <- rowSums(b^2)
  d2 <- outer(aa, bb, "+") - 2 * tcrossprod(a, b)
  d2[d2 < 0] <- 0
  sqrt(d2)
}

.assign_rotated_indices <- function(coord, coord_rot) {
  n <- nrow(coord)
  dist_mat <- .euclid_dist_matrix(coord, coord_rot)
  tmp <- dist_mat

  ref <- integer(n)
  rot <- integer(n)

  for (i in seq_len(n)) {
    row_min <- apply(tmp, 1, min, na.rm = TRUE)
    row_min[!is.finite(row_min)] <- -Inf
    ref_ix <- which(row_min == max(row_min))[1]
    row_vals <- tmp[ref_ix, ]
    rot_ix <- which(row_vals == min(row_vals, na.rm = TRUE))[1]

    ref[i] <- ref_ix
    rot[i] <- rot_ix

    tmp[, rot_ix] <- NA_real_
    tmp[ref_ix, ] <- 0
  }

  list(ref = ref, rot = rot)
}

.rotate_parcellation <- function(coord.l, coord.r, nrot = 10000L, verbose = FALSE) {
  if (!all(dim(coord.l)[2] == 3, dim(coord.r)[2] == 3)) {
    if (all(dim(coord.l)[1] == 3, dim(coord.r)[1] == 3)) {
      coord.l <- t(coord.l)
      coord.r <- t(coord.r)
    } else {
      .spin_stop("E_SPIN_COORD", "Coordinate matrices must have 3 columns (x/y/z).")
    }
  }

  nrot <- as.integer(nrot)
  if (!is.finite(nrot) || nrot < 1) {
    .spin_stop("E_SPIN_NPERM", "`n_perm` / `nrot` must be >= 1.")
  }

  nroi.l <- nrow(coord.l)
  nroi.r <- nrow(coord.r)
  nroi <- nroi.l + nroi.r

  perm.id <- matrix(0L, nrow = nroi, ncol = nrot)
  I1 <- diag(3)
  I1[1, 1] <- -1

  r <- 0L
  while (r < nrot) {
    A <- matrix(stats::rnorm(9), nrow = 3, ncol = 3)
    qrdec <- qr(A)
    TL <- qr.Q(qrdec)
    temp <- qr.R(qrdec)
    TL <- TL %*% diag(sign(diag(temp)))
    if (det(TL) < 0) {
      TL[, 1] <- -TL[, 1]
    }

    TR <- I1 %*% TL %*% I1
    coord.l.rot <- coord.l %*% TL
    coord.r.rot <- coord.r %*% TR

    idx_l <- .assign_rotated_indices(coord.l, coord.l.rot)
    idx_r <- .assign_rotated_indices(coord.r, coord.r.rot)

    ref.lr <- c(idx_l$ref, nroi.l + idx_r$ref)
    rot.lr <- c(idx_l$rot, nroi.l + idx_r$rot)
    rot.lr.sort <- rot.lr[order(ref.lr)]

    if (!all(sort(rot.lr.sort) == seq_len(nroi))) {
      .spin_stop("E_SPIN_PERM", "Generated permutation was invalid.")
    }

    if (!all(rot.lr.sort == seq_len(nroi))) {
      r <- r + 1L
      perm.id[, r] <- as.integer(rot.lr.sort)
      if (verbose && r %% 100 == 0) {
        message(sprintf("[spin] generated %d / %d permutations", r, nrot))
      }
    }
  }

  perm.id
}

#' Generate Spin Permutations from Hemisphere Centroids
#'
#' @description
#' Generates permutation indices for spherical spin tests using left/right
#' hemisphere centroids. This is the core permutation generator used by
#' \code{\link{spin_correlation}} when \code{perm_id} is not supplied.
#'
#' @param n_perm Integer. Number of permutations.
#' @param seed Optional integer seed for reproducibility.
#' @param lh_centroid_data Optional data frame of left hemisphere centroids.
#' @param rh_centroid_data Optional data frame of right hemisphere centroids.
#' @param verbose Logical. Print progress messages.
#'
#' @return Integer matrix of size \code{n_roi x n_perm}, where each column is a
#' permutation of ROI indices.
#' @export
spin_generate_permutations <- function(n_perm = 5000L,
                                       seed = NULL,
                                       lh_centroid_data = NULL,
                                       rh_centroid_data = NULL,
                                       verbose = FALSE) {
  n_perm <- as.integer(n_perm)
  if (!is.finite(n_perm) || n_perm < 1) {
    .spin_stop("E_SPIN_NPERM", "`n_perm` must be >= 1.")
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  if (is.null(lh_centroid_data) || is.null(rh_centroid_data)) {
    cen <- .load_default_centroids()
    if (is.null(lh_centroid_data)) lh_centroid_data <- cen$lh
    if (is.null(rh_centroid_data)) rh_centroid_data <- cen$rh
  }

  coord.l <- .extract_centroid_matrix(lh_centroid_data, "lh")
  coord.r <- .extract_centroid_matrix(rh_centroid_data, "rh")

  .rotate_parcellation(coord.l, coord.r, nrot = n_perm, verbose = verbose)
}

.spin_pair_p <- function(obs, null_xy, null_yx) {
  if (is.na(obs)) return(NA_real_)
  if (obs >= 0) {
    pxy <- mean(null_xy >= obs, na.rm = TRUE)
    pyx <- mean(null_yx >= obs, na.rm = TRUE)
  } else {
    pxy <- mean(null_xy <= obs, na.rm = TRUE)
    pyx <- mean(null_yx <= obs, na.rm = TRUE)
  }
  (pxy + pyx) / 2
}

.spin_p_matrix <- function(observed, null_xy, null_yx) {
  out <- matrix(NA_real_, nrow = nrow(observed), ncol = ncol(observed))
  for (i in seq_len(nrow(observed))) {
    for (j in seq_len(ncol(observed))) {
      out[i, j] <- .spin_pair_p(observed[i, j], null_xy[, i, j], null_yx[, i, j])
    }
  }
  dimnames(out) <- dimnames(observed)
  out
}

.spin_null_at_perm <- function(k, x_mat, y_mat, perm_id, method, use_cor) {
  idx <- perm_id[, k]
  xy <- stats::cor(x_mat[idx, , drop = FALSE], y_mat, method = method, use = use_cor)
  yx <- stats::cor(x_mat, y_mat[idx, , drop = FALSE], method = method, use = use_cor)
  list(k = k, xy = xy, yx = yx)
}

.spin_compute_null <- function(x_mat, y_mat, perm_id, method, use_cor,
                               parallel = FALSE, n_cores = 2L) {
  n_perm <- ncol(perm_id)
  n_x <- ncol(x_mat)
  n_y <- ncol(y_mat)

  null_xy <- array(NA_real_, dim = c(n_perm, n_x, n_y))
  null_yx <- array(NA_real_, dim = c(n_perm, n_x, n_y))

  idx <- seq_len(n_perm)
  if (!isTRUE(parallel)) {
    for (k in idx) {
      tmp <- .spin_null_at_perm(k, x_mat, y_mat, perm_id, method, use_cor)
      null_xy[tmp$k, , ] <- tmp$xy
      null_yx[tmp$k, , ] <- tmp$yx
    }
    return(list(null_xy = null_xy, null_yx = null_yx))
  }

  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1) n_cores <- 1L

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parLapply(
      cl, idx, .spin_null_at_perm,
      x_mat = x_mat, y_mat = y_mat, perm_id = perm_id,
      method = method, use_cor = use_cor
    )
  } else {
    res <- parallel::mclapply(
      idx, .spin_null_at_perm,
      x_mat = x_mat, y_mat = y_mat, perm_id = perm_id,
      method = method, use_cor = use_cor,
      mc.cores = n_cores
    )
  }

  for (tmp in res) {
    null_xy[tmp$k, , ] <- tmp$xy
    null_yx[tmp$k, , ] <- tmp$yx
  }
  list(null_xy = null_xy, null_yx = null_yx)
}

#' Spin Correlation Test for Two Maps or Multi-Map Matrices
#'
#' @description
#' Runs a spatial permutation (spin) correlation test for:
#' \itemize{
#'   \item Two maps: \code{vector vs vector}
#'   \item Multi-map: \code{matrix/data.frame vs matrix/data.frame}
#' }
#'
#' The function is side-effect free and returns a structured object. Optional
#' saving and plotting are handled by separate helpers.
#'
#' @param x Numeric vector/matrix/data.frame. Rows must represent aligned ROIs.
#' @param y Optional numeric vector/matrix/data.frame. If \code{NULL}, performs
#' pairwise testing within \code{x}.
#' @param perm_id Optional permutation matrix (\code{n_roi x n_perm}). If
#' \code{NULL}, permutations are generated from centroid data.
#' @param n_perm Integer. Number of permutations (used when \code{perm_id} is
#' \code{NULL}).
#' @param method Correlation method: \code{"pearson"} or \code{"spearman"}.
#' @param na_action Missing-data handling: \code{"error"} or \code{"pairwise"}.
#' @param seed Optional integer seed.
#' @param parallel Logical. Compute null distributions in parallel.
#' @param n_cores Integer. Number of cores when \code{parallel = TRUE}.
#' @param lh_centroid_data Optional LH centroid data frame (if generating
#' permutations).
#' @param rh_centroid_data Optional RH centroid data frame (if generating
#' permutations).
#' @param atlas Character metadata tag (e.g. \code{"desikan_68"}).
#' @param return_null Logical. Include null distributions in output.
#' @param return_perm_id Logical. Include permutation matrix in output.
#'
#' @return An object of class \code{"spin_correlation_result"} with fields:
#' \itemize{
#'   \item \code{observed}: observed correlation matrix
#'   \item \code{p_value}: spin p-value matrix
#'   \item \code{null_xy}, \code{null_yx}: null arrays (optional)
#'   \item \code{summary}: tidy long table
#'   \item \code{meta}: run metadata (method, n_perm, seed, atlas, timestamp, ...)
#' }
#' @export
spin_correlation <- function(x, y = NULL,
                             perm_id = NULL,
                             n_perm = 10000L,
                             method = c("pearson", "spearman"),
                             na_action = c("error", "pairwise"),
                             seed = NULL,
                             parallel = FALSE,
                             n_cores = 2L,
                             lh_centroid_data = NULL,
                             rh_centroid_data = NULL,
                             atlas = "desikan_68",
                             return_null = TRUE,
                             return_perm_id = FALSE) {
  method <- match.arg(method)
  na_action <- match.arg(na_action)

  x_mat <- .as_numeric_matrix(x, "x")
  y_was_null <- is.null(y)
  y_mat <- if (y_was_null) x_mat else .as_numeric_matrix(y, "y")

  if (nrow(x_mat) != nrow(y_mat)) {
    .spin_stop(
      "E_SPIN_SHAPE",
      sprintf("Row mismatch: nrow(x)=%d vs nrow(y)=%d.", nrow(x_mat), nrow(y_mat))
    )
  }
  n_roi <- nrow(x_mat)

  if (na_action == "error" && (anyNA(x_mat) || anyNA(y_mat))) {
    .spin_stop("E_SPIN_NA", "Input contains NA values; set `na_action = 'pairwise'` to allow NA.")
  }
  use_cor <- if (na_action == "pairwise") "pairwise.complete.obs" else "everything"

  perm_generated <- FALSE
  if (is.null(perm_id)) {
    perm_generated <- TRUE
    perm_id <- spin_generate_permutations(
      n_perm = n_perm,
      seed = seed,
      lh_centroid_data = lh_centroid_data,
      rh_centroid_data = rh_centroid_data,
      verbose = FALSE
    )
  } else {
    if (!is.matrix(perm_id)) {
      .spin_stop("E_SPIN_PERM", "`perm_id` must be a matrix.")
    }
    perm_id <- apply(perm_id, 2, as.integer)
    if (nrow(perm_id) != n_roi) {
      .spin_stop("E_SPIN_PERM", sprintf(
        "`perm_id` row count (%d) must equal ROI count (%d).", nrow(perm_id), n_roi
      ))
    }
    if (any(perm_id < 1L) || any(perm_id > n_roi)) {
      .spin_stop("E_SPIN_PERM", "`perm_id` contains indices out of range.")
    }
    if (!is.null(seed)) set.seed(as.integer(seed))
  }

  observed <- stats::cor(x_mat, y_mat, method = method, use = use_cor)
  if (!is.matrix(observed)) observed <- matrix(observed, nrow = ncol(x_mat), ncol = ncol(y_mat))
  rownames(observed) <- colnames(x_mat)
  colnames(observed) <- colnames(y_mat)

  nulls <- .spin_compute_null(
    x_mat = x_mat,
    y_mat = y_mat,
    perm_id = perm_id,
    method = method,
    use_cor = use_cor,
    parallel = parallel,
    n_cores = n_cores
  )
  p_value <- .spin_p_matrix(observed, nulls$null_xy, nulls$null_yx)

  if (y_was_null && ncol(x_mat) == ncol(y_mat) && identical(colnames(x_mat), colnames(y_mat))) {
    diag(p_value) <- NA_real_
    diag(observed) <- 1
  }

  tbl <- spin_results_table(list(observed = observed, p_value = p_value))

  meta <- list(
    method = method,
    n_perm = ncol(perm_id),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    atlas = atlas,
    n_roi = n_roi,
    n_x_maps = ncol(x_mat),
    n_y_maps = ncol(y_mat),
    symmetric = y_was_null,
    na_action = na_action,
    parallel = isTRUE(parallel),
    n_cores = as.integer(n_cores),
    permutations_generated = perm_generated,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  out <- list(
    observed = observed,
    p_value = p_value,
    null_xy = if (isTRUE(return_null)) nulls$null_xy else NULL,
    null_yx = if (isTRUE(return_null)) nulls$null_yx else NULL,
    summary = tbl,
    meta = meta,
    perm_id = if (isTRUE(return_perm_id)) perm_id else NULL
  )
  class(out) <- c("spin_correlation_result", "list")
  out
}

#' @export
print.spin_correlation_result <- function(x, ...) {
  cat("<spin_correlation_result>\n")
  cat(sprintf("  method: %s\n", x$meta$method))
  cat(sprintf("  n_perm: %d\n", x$meta$n_perm))
  cat(sprintf("  n_roi: %d\n", x$meta$n_roi))
  cat(sprintf("  maps: %d x %d\n", x$meta$n_x_maps, x$meta$n_y_maps))
  cat(sprintf("  symmetric: %s\n", ifelse(isTRUE(x$meta$symmetric), "yes", "no")))
  invisible(x)
}

#' Legacy Wrapper: Generate 68-Region Spin Permutations
#'
#' @description
#' Backward-compatible wrapper for legacy code. Prefer
#' \code{\link{spin_generate_permutations}}.
#'
#' @param nrot Number of rotations/permutations.
#' @param seed Optional seed.
#' @return Integer permutation matrix.
#' @export
get_68perm.id <- function(nrot = 10000L, seed = 666L) {
  warning("`get_68perm.id()` is legacy; use `spin_generate_permutations()`.", call. = FALSE)
  spin_generate_permutations(n_perm = nrot, seed = seed)
}

#' Legacy Wrapper: Spin P Matrix for a Correlation Table
#'
#' @description
#' Backward-compatible wrapper for legacy code. Prefer
#' \code{\link{spin_correlation}} and use \code{$p_value}.
#'
#' @param df2cor Numeric matrix/data.frame with maps in columns.
#' @param perm.id Permutation matrix.
#' @param corr.type Correlation method.
#' @param parallel Logical.
#' @param n_cores Integer.
#' @return Matrix of spin p-values.
#' @export
caculate_spin_p <- function(df2cor, perm.id, corr.type = "pearson",
                            parallel = FALSE, n_cores = 2L) {
  warning("`caculate_spin_p()` is legacy; use `spin_correlation()`.", call. = FALSE)
  res <- spin_correlation(
    x = df2cor,
    y = NULL,
    perm_id = perm.id,
    method = corr.type,
    parallel = parallel,
    n_cores = n_cores,
    return_null = FALSE
  )
  res$p_value
}

#' Legacy Wrapper: Two-Map Spin P-Value
#'
#' @description
#' Backward-compatible wrapper for legacy code. Prefer
#' \code{\link{spin_correlation}} for two-map inputs.
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param perm.id Permutation matrix.
#' @param corr.type Correlation method.
#' @return Numeric p-value.
#' @export
perm.sphere.p <- function(x, y, perm.id, corr.type = "spearman") {
  warning("`perm.sphere.p()` is legacy; use `spin_correlation()`.", call. = FALSE)
  res <- spin_correlation(
    x = x,
    y = y,
    perm_id = perm.id,
    method = corr.type,
    return_null = FALSE
  )
  as.numeric(res$p_value[1, 1])
}
