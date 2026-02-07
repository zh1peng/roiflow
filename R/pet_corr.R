# PET map correlation utilities built on spin permutations.

.pet_stop <- function(code, msg) {
  stop(sprintf("%s: %s", code, msg), call. = FALSE)
}

.load_pet_maps <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("pet_maps_dk68", package = "roiflow", envir = env)
  utils::data("pet_maps_meta", package = "roiflow", envir = env)
  if (!exists("pet_maps_dk68", envir = env, inherits = FALSE)) {
    .pet_stop("E_PET_DATA_MISSING", "Packaged PET maps (`pet_maps_dk68`) are not available.")
  }
  if (!exists("pet_maps_meta", envir = env, inherits = FALSE)) {
    .pet_stop("E_PET_DATA_MISSING", "Packaged PET metadata (`pet_maps_meta`) are not available.")
  }
  list(
    maps = get("pet_maps_dk68", envir = env, inherits = FALSE),
    meta = get("pet_maps_meta", envir = env, inherits = FALSE)
  )
}

.pet_resolve_centroids <- function(centroids = NULL) {
  if (is.null(centroids)) {
    return(.load_default_centroids())
  }
  if (!is.list(centroids)) {
    .pet_stop("E_PET_CENTROIDS", "`centroids` must be NULL or a list with `lh` and `rh`.")
  }
  if (is.null(names(centroids)) && length(centroids) == 2) {
    names(centroids) <- c("lh", "rh")
  }
  if (is.null(centroids$lh) || is.null(centroids$rh)) {
    .pet_stop("E_PET_CENTROIDS", "`centroids` must provide both `lh` and `rh` data.")
  }
  list(lh = centroids$lh, rh = centroids$rh)
}

.pet_as_map_matrix <- function(x, arg) {
  if (is.vector(x) && !is.list(x)) {
    out <- matrix(as.numeric(x), ncol = 1)
    colnames(out) <- arg
    nm <- names(x)
    if (!is.null(nm) && length(nm) == length(x)) {
      rownames(out) <- as.character(nm)
    }
    return(out)
  }
  .as_numeric_matrix(x, arg = arg)
}

.pet_centroid_labels <- function(df, fallback_prefix) {
  if (!is.data.frame(df)) {
    .pet_stop("E_PET_CENTROIDS", "Centroid entries must be data.frames.")
  }
  if ("Row" %in% names(df)) {
    return(as.character(df$Row))
  }
  if ("row" %in% names(df)) {
    return(as.character(df$row))
  }
  paste0(fallback_prefix, "_", seq_len(nrow(df)))
}

.pet_expected_roi <- function(hemi, lh_centroid_data, rh_centroid_data) {
  lh_names <- .pet_centroid_labels(lh_centroid_data, "lh")
  rh_names <- .pet_centroid_labels(rh_centroid_data, "rh")
  if (hemi == "both") return(c(lh_names, rh_names))
  if (hemi == "left") return(lh_names)
  rh_names
}

.pet_align_rows <- function(mat, expected_roi, arg) {
  rn <- rownames(mat)
  if (is.null(rn)) {
    if (nrow(mat) != length(expected_roi)) {
      .pet_stop(
        "E_PET_ROI",
        sprintf(
          "`%s` has %d rows but %d are expected for selected centroids/hemi.",
          arg, nrow(mat), length(expected_roi)
        )
      )
    }
    return(mat)
  }

  missing <- setdiff(expected_roi, rn)
  if (length(missing) > 0) {
    .pet_stop(
      "E_PET_ROI",
      sprintf("`%s` row names are missing expected ROI(s): %s", arg, paste(missing, collapse = ", "))
    )
  }
  mat[expected_roi, , drop = FALSE]
}

.pet_rotate_single_hemi <- function(coord, nrot = 10000L, verbose = FALSE) {
  if (ncol(coord) != 3) {
    if (nrow(coord) == 3) {
      coord <- t(coord)
    } else {
      .pet_stop("E_PET_CENTROIDS", "Hemisphere centroid coordinates must have 3 columns (x/y/z).")
    }
  }

  nrot <- as.integer(nrot)
  if (!is.finite(nrot) || nrot < 1) {
    .pet_stop("E_PET_NPERM", "`n_perm` must be >= 1.")
  }

  n_roi <- nrow(coord)
  perm_id <- matrix(0L, nrow = n_roi, ncol = nrot)
  r <- 0L
  while (r < nrot) {
    A <- matrix(stats::rnorm(9), nrow = 3, ncol = 3)
    qrdec <- qr(A)
    Tm <- qr.Q(qrdec)
    temp <- qr.R(qrdec)
    Tm <- Tm %*% diag(sign(diag(temp)))
    if (det(Tm) < 0) {
      Tm[, 1] <- -Tm[, 1]
    }

    coord_rot <- coord %*% Tm
    idx <- .assign_rotated_indices(coord, coord_rot)
    perm <- idx$rot[order(idx$ref)]

    if (!all(sort(perm) == seq_len(n_roi))) {
      .pet_stop("E_PET_PERM", "Generated hemisphere permutation was invalid.")
    }
    if (!all(perm == seq_len(n_roi))) {
      r <- r + 1L
      perm_id[, r] <- as.integer(perm)
      if (verbose && r %% 100 == 0) {
        message(sprintf("[pet-spin] generated %d / %d permutations", r, nrot))
      }
    }
  }
  perm_id
}

.pet_generate_permutations <- function(n_perm, seed, hemi, lh_centroid_data, rh_centroid_data) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  if (hemi == "both") {
    return(spin_generate_permutations(
      n_perm = n_perm,
      seed = NULL,
      lh_centroid_data = lh_centroid_data,
      rh_centroid_data = rh_centroid_data,
      verbose = FALSE
    ))
  }

  if (hemi == "left") {
    coord <- .extract_centroid_matrix(lh_centroid_data, "lh")
  } else {
    coord <- .extract_centroid_matrix(rh_centroid_data, "rh")
  }
  .pet_rotate_single_hemi(coord, nrot = n_perm, verbose = FALSE)
}

.pet_select_maps <- function(pet, pet_select) {
  pkg <- .load_pet_maps()
  pet_from_package <- is.null(pet)
  pet_mat <- if (pet_from_package) pkg$maps else .pet_as_map_matrix(pet, arg = "pet")

  if (!is.null(pet_select)) {
    if (!is.character(pet_select)) {
      .pet_stop("E_PET_SELECT", "`pet_select` must be NULL or a character vector.")
    }
    pet_select <- unique(as.character(pet_select))
    invalid <- setdiff(pet_select, colnames(pet_mat))
    if (length(invalid) > 0) {
      .pet_stop(
        "E_PET_SELECT",
        sprintf("Unknown PET map name(s): %s", paste(invalid, collapse = ", "))
      )
    }
    pet_mat <- pet_mat[, pet_select, drop = FALSE]
  }

  meta <- pkg$meta
  meta <- meta[meta$pet_map %in% colnames(pet_mat), , drop = FALSE]
  if (nrow(meta) > 0) {
    meta <- meta[match(colnames(pet_mat), meta$pet_map), , drop = FALSE]
  }

  list(
    pet = pet_mat,
    meta = meta,
    from_package = pet_from_package
  )
}

.pet_pair_stat_row <- function(x, y, method) {
  keep <- stats::complete.cases(x, y)
  n_obs <- sum(keep)
  if (n_obs < 3) {
    return(list(
      n_obs = n_obs,
      p_param = NA_real_,
      stat = NA_real_,
      stat_name = NA_character_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  ci <- ct$conf.int %||% c(NA_real_, NA_real_)
  list(
    n_obs = n_obs,
    p_param = ct$p.value %||% NA_real_,
    stat = if (is.null(ct$statistic)) NA_real_ else as.numeric(ct$statistic[[1]]),
    stat_name = if (is.null(ct$statistic)) NA_character_ else names(ct$statistic)[1],
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

.pet_result_row <- function(k, summary_tbl, brain_mat, pet_mat, null_xy, null_yx, method) {
  brain_name <- summary_tbl$x_map[k]
  pet_name <- summary_tbl$y_map[k]

  i <- match(brain_name, colnames(brain_mat))
  j <- match(pet_name, colnames(pet_mat))
  if (is.na(i) || is.na(j)) {
    .pet_stop("E_PET_INTERNAL", "Could not match map names to input matrices.")
  }

  x <- brain_mat[, i]
  y <- pet_mat[, j]
  pair_stat <- .pet_pair_stat_row(x, y, method = method)

  null_xy_vec <- null_xy[, i, j]
  null_yx_vec <- null_yx[, i, j]
  null_combined <- c(null_xy_vec, null_yx_vec)
  q <- stats::quantile(null_combined, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)

  list(
    result = data.frame(
      brain_map = brain_name,
      pet_map = pet_name,
      n_obs = pair_stat$n_obs,
      r_obs = summary_tbl$observed_r[k],
      r_sq = summary_tbl$observed_r[k]^2,
      fisher_z = if (abs(summary_tbl$observed_r[k]) < 1) atanh(summary_tbl$observed_r[k]) else NA_real_,
      p_spin = summary_tbl$p_spin[k],
      p_param = pair_stat$p_param,
      stat = pair_stat$stat,
      stat_name = pair_stat$stat_name,
      ci_low = pair_stat$ci_low,
      ci_high = pair_stat$ci_high,
      null_q025 = q[1],
      null_q50 = q[2],
      null_q975 = q[3],
      null_mean = mean(null_combined, na.rm = TRUE),
      null_sd = stats::sd(null_combined, na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    null = list(
      brain_map = brain_name,
      pet_map = pet_name,
      null_xy = null_xy_vec,
      null_yx = null_yx_vec,
      null_combined = null_combined
    )
  )
}

#' List Available Packaged PET Maps
#'
#' @description
#' Returns a metadata table for PET maps bundled with \pkg{roiflow}, including
#' map names, target labels, modality, and parcellation assumptions.
#'
#' @return A data frame with one row per packaged PET map.
#' @export
#'
#' @examples
#' head(pet_available())
pet_available <- function() {
  obj <- .load_pet_maps()
  out <- obj$meta
  out$n_roi <- nrow(obj$maps)
  out
}

#' Correlate Brain Map(s) with PET Map(s) Using Spin Tests
#'
#' @description
#' Runs observed + spin-test correlations between one or more query brain maps
#' and one or more PET maps.
#'
#' PET maps can be supplied directly via \code{pet}, or selected from packaged
#' PET data via \code{pet_select}. By default, \code{pet = NULL} uses all
#' packaged PET maps.
#'
#' @param brain Numeric vector, matrix, or data.frame. Columns are brain maps.
#' @param pet Optional numeric vector/matrix/data.frame of PET maps.
#' @param pet_select Optional character vector of PET map names to select.
#' @param ... Optional named arguments. Supported: \code{perm_id}.
#' @param n_perm Integer number of spin permutations.
#' @param method Correlation method: \code{"pearson"} or \code{"spearman"}.
#' @param seed Optional integer seed.
#' @param centroids Optional list with \code{lh} and \code{rh} centroid data.
#' @param hemi Hemisphere selection: \code{"both"}, \code{"left"}, \code{"right"}.
#' @param na_action Missing data policy: \code{"omit"} or \code{"error"}.
#' @param return_null Logical. If \code{TRUE}, return null vectors.
#' @param p_adjust Multiple-comparison adjustment for \code{p_spin}:
#' \code{"none"} or \code{"fdr"}.
#' @param p_adjust_scope Adjustment scope when \code{p_adjust = "fdr"}:
#' \code{"global"} or \code{"brain_map"}.
#' @param parallel Logical. Parallel null computation.
#' @param n_cores Integer. Number of cores for parallel computation.
#' @param save Logical. If \code{TRUE}, write csv/rds output files.
#' @param out_dir Output directory used when \code{save = TRUE}.
#' @param out_prefix Output file prefix used when \code{save = TRUE}.
#'
#' @return Object of class \code{"pet_corr_result"} with:
#' \itemize{
#'   \item \code{results}: tidy results table
#'   \item \code{nulls}: optional null vectors per map pair
#'   \item \code{meta}: metadata, settings, and provenance
#'   \item \code{inputs}: aligned matrices used for correlation
#' }
#' @export
#'
#' @examples
#' data(pet_maps_dk68)
#' brain <- pet_maps_dk68[, "D1"] + stats::rnorm(nrow(pet_maps_dk68), sd = 0.05)
#' perm_id <- replicate(25, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
#' res <- pet_corr(
#'   brain = brain,
#'   pet_select = c("D1", "D2"),
#'   n_perm = 25,
#'   seed = 123,
#'   perm_id = perm_id
#' )
#' head(res$results)
pet_corr <- function(brain,
                     pet = NULL,
                     pet_select = NULL,
                     ...,
                     n_perm = 5000L,
                     method = c("pearson", "spearman"),
                     seed = NULL,
                     centroids = NULL,
                     hemi = c("both", "left", "right"),
                     na_action = c("omit", "error"),
                     return_null = FALSE,
                     p_adjust = c("none", "fdr"),
                     p_adjust_scope = c("global", "brain_map"),
                     parallel = FALSE,
                     n_cores = 2L,
                     save = FALSE,
                     out_dir = ".",
                     out_prefix = "pet_corr_results") {
  method <- match.arg(method)
  hemi <- match.arg(hemi)
  na_action <- match.arg(na_action)
  p_adjust <- match.arg(p_adjust)
  p_adjust_scope <- match.arg(p_adjust_scope)

  dots <- list(...)
  if (length(dots) > 0) {
    if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
      .pet_stop("E_PET_DOTS", "All `...` arguments must be named.")
    }
    bad <- setdiff(names(dots), "perm_id")
    if (length(bad) > 0) {
      .pet_stop("E_PET_DOTS", sprintf("Unsupported `...` argument(s): %s", paste(bad, collapse = ", ")))
    }
  }

  n_perm <- as.integer(n_perm)
  if (!is.finite(n_perm) || n_perm < 1) {
    .pet_stop("E_PET_NPERM", "`n_perm` must be >= 1.")
  }
  if (!is.null(seed)) seed <- as.integer(seed)

  centroid_data <- .pet_resolve_centroids(centroids = centroids)
  expected_roi <- .pet_expected_roi(hemi, centroid_data$lh, centroid_data$rh)

  brain_mat <- .pet_as_map_matrix(brain, arg = "brain")
  pet_info <- .pet_select_maps(pet = pet, pet_select = pet_select)
  pet_mat <- pet_info$pet

  brain_mat <- .pet_align_rows(brain_mat, expected_roi = expected_roi, arg = "brain")
  pet_mat <- .pet_align_rows(pet_mat, expected_roi = expected_roi, arg = "pet")
  if (nrow(brain_mat) != nrow(pet_mat)) {
    .pet_stop(
      "E_PET_SHAPE",
      sprintf("Row mismatch after alignment: brain=%d, pet=%d.", nrow(brain_mat), nrow(pet_mat))
    )
  }

  if (na_action == "error" && (anyNA(brain_mat) || anyNA(pet_mat))) {
    .pet_stop("E_PET_NA", "NA values detected. Set `na_action = 'omit'` to allow pairwise NA omission.")
  }

  perm_id <- dots$perm_id %||% NULL
  if (!is.null(perm_id)) {
    if (!is.matrix(perm_id)) {
      .pet_stop("E_PET_PERM", "`perm_id` must be a matrix.")
    }
    if (nrow(perm_id) != nrow(brain_mat)) {
      .pet_stop(
        "E_PET_PERM",
        sprintf("`perm_id` has %d rows; expected %d.", nrow(perm_id), nrow(brain_mat))
      )
    }
    perm_id <- apply(perm_id, 2, as.integer)
  } else {
    perm_id <- .pet_generate_permutations(
      n_perm = n_perm,
      seed = seed,
      hemi = hemi,
      lh_centroid_data = centroid_data$lh,
      rh_centroid_data = centroid_data$rh
    )
  }

  spin_res <- spin_correlation(
    x = brain_mat,
    y = pet_mat,
    perm_id = perm_id,
    method = method,
    na_action = if (na_action == "omit") "pairwise" else "error",
    seed = seed,
    parallel = parallel,
    n_cores = n_cores,
    return_null = TRUE,
    atlas = sprintf("pet_%s", hemi)
  )

  summary_tbl <- spin_results_table(spin_res, upper_only = FALSE, drop_diag = FALSE)
  n_pairs <- nrow(summary_tbl)
  rows <- vector("list", n_pairs)
  null_rows <- vector("list", n_pairs)
  for (k in seq_len(n_pairs)) {
    tmp <- .pet_result_row(
      k = k,
      summary_tbl = summary_tbl,
      brain_mat = brain_mat,
      pet_mat = pet_mat,
      null_xy = spin_res$null_xy,
      null_yx = spin_res$null_yx,
      method = method
    )
    rows[[k]] <- tmp$result
    null_rows[[k]] <- tmp$null
  }

  results <- do.call(rbind, rows)
  results$n_perm <- ncol(perm_id)
  results$method <- method
  results$seed <- if (is.null(seed)) NA_integer_ else seed
  results$hemi <- hemi

  if (p_adjust == "fdr") {
    if (p_adjust_scope == "global") {
      results$p_fdr <- stats::p.adjust(results$p_spin, method = "fdr")
    } else {
      results$p_fdr <- NA_real_
      by_map <- split(seq_len(nrow(results)), results$brain_map)
      for (idx in by_map) {
        results$p_fdr[idx] <- stats::p.adjust(results$p_spin[idx], method = "fdr")
      }
    }
  } else {
    results$p_fdr <- NA_real_
  }

  null_tbl <- NULL
  if (isTRUE(return_null)) {
    null_tbl <- data.frame(
      brain_map = vapply(null_rows, `[[`, character(1), "brain_map"),
      pet_map = vapply(null_rows, `[[`, character(1), "pet_map"),
      stringsAsFactors = FALSE
    )
    null_tbl$null_xy <- I(lapply(null_rows, `[[`, "null_xy"))
    null_tbl$null_yx <- I(lapply(null_rows, `[[`, "null_yx"))
    null_tbl$null_combined <- I(lapply(null_rows, `[[`, "null_combined"))
  }

  meta <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package_version = as.character(utils::packageVersion("roiflow")),
    n_roi = nrow(brain_mat),
    n_brain_maps = ncol(brain_mat),
    n_pet_maps = ncol(pet_mat),
    brain_map_names = colnames(brain_mat),
    pet_map_names = colnames(pet_mat),
    method = method,
    n_perm = ncol(perm_id),
    seed = if (is.null(seed)) NA_integer_ else seed,
    hemi = hemi,
    na_action = na_action,
    p_adjust = p_adjust,
    p_adjust_scope = p_adjust_scope,
    centroids = list(
      source = if (is.null(centroids)) "package_default" else "user_supplied",
      lh_n = nrow(centroid_data$lh),
      rh_n = nrow(centroid_data$rh)
    ),
    pet_source = if (isTRUE(pet_info$from_package)) "package_data" else "user_input",
    pet_metadata = pet_info$meta
  )

  out <- list(
    results = results,
    nulls = null_tbl,
    meta = meta,
    inputs = list(
      brain = brain_mat,
      pet = pet_mat
    ),
    call = match.call()
  )
  class(out) <- c("pet_corr_result", "list")

  if (isTRUE(save)) {
    out$meta$saved_files <- list(
      csv = file.path(out_dir, paste0(out_prefix, ".csv")),
      rds = file.path(out_dir, paste0(out_prefix, ".rds"))
    )
    files <- pet_corr_save(
      pet_corr_obj = out,
      out_dir = out_dir,
      out_prefix = out_prefix,
      include_null = isTRUE(return_null)
    )
    out$meta$saved_files <- files
  }
  out
}

#' Save PET Correlation Outputs
#'
#' @description
#' Saves \code{pet_corr()} outputs to disk (CSV + RDS). Side effects are
#' opt-in via this helper or \code{pet_corr(save = TRUE)}.
#'
#' @param pet_corr_obj A \code{pet_corr_result} object.
#' @param out_dir Output directory.
#' @param out_prefix File prefix. Defaults to \code{"pet_corr_results"}.
#' @param include_null Logical. Save null vectors in RDS when available.
#'
#' @return Invisible named list of written file paths.
#' @export
#'
#' @examples
#' data(pet_maps_dk68)
#' perm_id <- replicate(10, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
#' res <- pet_corr(
#'   brain = pet_maps_dk68[, "D1"],
#'   pet_select = "D2",
#'   perm_id = perm_id,
#'   n_perm = 10
#' )
#' tmp <- tempdir()
#' files <- pet_corr_save(res, out_dir = tmp)
pet_corr_save <- function(pet_corr_obj,
                          out_dir = ".",
                          out_prefix = "pet_corr_results",
                          include_null = TRUE) {
  if (!inherits(pet_corr_obj, "pet_corr_result")) {
    .pet_stop("E_PET_RESULT", "`pet_corr_obj` must be a `pet_corr_result`.")
  }
  if (is.null(out_dir) || !nzchar(out_dir)) {
    .pet_stop("E_PET_IO", "`out_dir` must be a non-empty path.")
  }
  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    .pet_stop("E_PET_IO", "`out_prefix` must be a non-empty string.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(out_dir, paste0(out_prefix, ".csv"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".rds"))

  utils::write.csv(pet_corr_obj$results, csv_path, row.names = FALSE)

  obj_to_save <- pet_corr_obj
  if (!isTRUE(include_null)) {
    obj_to_save$nulls <- NULL
  }
  saveRDS(obj_to_save, rds_path)

  invisible(list(csv = csv_path, rds = rds_path))
}

#' @export
print.pet_corr_result <- function(x, ...) {
  cat("<pet_corr_result>\n")
  cat(sprintf("  map pairs: %d\n", nrow(x$results)))
  cat(sprintf("  brain maps: %d\n", x$meta$n_brain_maps))
  cat(sprintf("  pet maps: %d\n", x$meta$n_pet_maps))
  cat(sprintf("  n_perm: %d\n", x$meta$n_perm))
  cat(sprintf("  method: %s\n", x$meta$method))
  invisible(x)
}
