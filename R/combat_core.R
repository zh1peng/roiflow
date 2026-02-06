################################################################################
#                                                                              #
# combat_fit and combat_apply functions                                        #
# =====================================                                        #
#                                                                              #
# Authors:                                                                     #
#                                                                              #
#  1) The original ComBat function was in the sva package that can be found at #
#     https://bioconductor.org/packages/release/bioc/html/sva.html             #
#                                                                              #
#  2) First modification by Jean-Philippe Fortin for the harmonization of MRI  #
#     data. Please cite https://10.1016/j.neuroimage.2017.08.047               #
#                                                                              #
#  3) Second modification by Joaquim Radua to separate functions for fitting   #
#     and applying the harmonization, allow missings and constant rows and mi- #
#     nor changes in the arguments of the functions to facilitate their use.   #
#     Please cite "Increased power by harmonizing structural MRI site diffe-   #
#     rences with the ComBat batch adjustment method in ENIGMA".                #
#                                                                              #
# The original and present code is under the Artistic License 2.0. If using    #
# this code, make sure you agree and accept this license.                      #
#                                                                              #
################################################################################


# Calculate some characteristics on the batches (patched: store batch/mod)
.combat_tmp1 <- function(dat, batch, levels_batch, mod) {
  batchmod <- model.matrix(~ -1 + batch)

  n.batch <- nlevels(batch)
  batches <- lapply(seq_len(n.batch), function(i) which(batch == levels_batch[i]))
  n.batches <- vapply(batches, length, integer(1))
  n.array <- sum(n.batches)

  design <- cbind(batchmod, mod)

  # Drop intercept columns if present (all ones)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check, drop = FALSE])

  batch_cols <- colnames(batchmod)
  design_cols <- colnames(design)
  batch_col_idx <- match(batch_cols, design_cols)
  if (any(is.na(batch_col_idx))) {
    stop("Internal error: batch columns dropped from design; check batch levels.")
  }
  batch.design <- design[, batch_col_idx, drop = FALSE]

  list(
    dat = dat,
    batch = batch,  # NEW
    mod = mod,      # NEW
    batchmod = batchmod,
    n.batch = n.batch,
    batches = batches,
    n.batches = n.batches,
    n.array = n.array,
    design = design,
    design_cols = design_cols,
    batch.col.idx = batch_col_idx,
    batch.design = batch.design
  )
}

# Estimate B.hat, grand.mean and var.pooled
.combat_tmp2 <- function (tmp1, verbose = TRUE) {
  # Number of covariates or covariate levels
  if (verbose) {
    cat(
      "[combat] Adjusting for",
      ncol(tmp1$design) - ncol(tmp1$batch.design),
      "covariate(s) or covariate level(s)\n"
    )
  }
  # Check if the design is confounded
  n_batch_cols <- ncol(tmp1$batch.design)
  n_design <- ncol(tmp1$design)
  n_cov_cols <- n_design - n_batch_cols
  if (qr(tmp1$design)$rank < n_design) {
    if (n_cov_cols == 1) {
      stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if (n_cov_cols > 1) {
      design_cov <- as.matrix(tmp1$design[, -tmp1$batch.col.idx, drop = FALSE])
      if ((qr(design_cov)$rank < ncol(design_cov))) {
        stop("The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.")
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
    stop("Design matrix is rank deficient. Please check batch and covariates.")
  }
  # Standardize data across features
  B.hat <- solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% t(as.matrix(tmp1$dat))
  # Standarization Model
  grand.mean <- t(tmp1$n.batches / tmp1$n.array) %*% B.hat[tmp1$batch.col.idx, , drop = FALSE]
  var.pooled <- ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% rep(1 / tmp1$n.array, tmp1$n.array)
  return(list(
    B.hat = B.hat,
    grand.mean = grand.mean,
    var.pooled = var.pooled
  ))
}
# Standardize data
.combat_tmp3 <- function (dat, tmp1, tmp2, verbose = TRUE) {
  if (verbose) {
    cat("[combat] Standardizing data across features\n")
  }
  stand.mean <- t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
  if (!is.null(tmp1$design)) {
    tmp <- tmp1$design
    tmp[, tmp1$batch.col.idx] <- 0
    stand.mean <- stand.mean + t(tmp %*% tmp2$B.hat)
  } 
  s.data <- (dat-stand.mean) / (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))
  return(list(
    stand.mean = stand.mean,
    s.data = s.data
  ))
}
# Fit L/S model 
.combat_tmp4 <- function (tmp1, tmp2, tmp3, eb = TRUE, verbose = TRUE) {
  # Get regression batch effect parameters
  if (eb) {
    if (verbose) {
      cat("[combat] Fitting L/S model and finding priors\n")
    }
  } else {
    if (verbose) {
      cat("[combat] Fitting L/S model\n")
    }
  }
  gamma.hat <- solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
  delta.hat <- NULL
  for (i in tmp1$batches) {
    delta.hat <- rbind(delta.hat, apply(tmp3$s.data[,i], 1, var, na.rm = T))
  }
  if (!is.null(delta.hat)) {
    delta.hat <- as.matrix(delta.hat)
    delta_bad <- !is.finite(delta.hat) | delta.hat <= 0
    if (any(delta_bad)) {
      delta_floor <- min(delta.hat[delta.hat > 0 & is.finite(delta.hat)], na.rm = TRUE)
      if (!is.finite(delta_floor)) {
        delta_floor <- 1e-6
      }
      delta.hat[delta_bad] <- delta_floor
      warning(
        sprintf("W_DELTA_HAT_FLOOR: Replaced %d non-positive/NA batch variances with %g.",
                sum(delta_bad), delta_floor),
        call. = FALSE
      )
    }
  }
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb) {
    # Find Priors
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, .aprior)
    b.prior <- apply(delta.hat, 1, .bprior)
    # Find EB batch adjustments
    if (verbose) {
      cat("[combat] Finding parametric adjustments\n")
    }
    for (i in 1:tmp1$n.batch) {
      temp <- .it.sol(tmp3$s.data[,tmp1$batches[[i]]], gamma.hat[i,], delta.hat[i,], gamma.bar[i], t2[i], a.prior[i], b.prior[i])
      gamma.star <- rbind(gamma.star, temp[1,])
      delta.star <- rbind(delta.star, temp[2,])
    }
  } 
  return(list(
    gamma.hat = gamma.hat,
    delta.hat = delta.hat, 
    gamma.star = gamma.star,
    delta.star = delta.star, 
    gamma.bar = gamma.bar,
    t2 = t2,
    a.prior = a.prior,
    b.prior = b.prior
  ))
}
# Adjust the data
.combat_tmp5 <- function (tmp1, tmp2, tmp3, tmp4, eb = TRUE, verbose = TRUE) {
  # Normalize the data
  if (verbose) {
    cat("[combat] Adjusting the data\n")
  }
  bayesdata <- tmp3$s.data
  j <- 1
  for (i in tmp1$batches) {
    if (eb) {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.star)) / (sqrt(tmp4$delta.star[j,]) %*% t(rep(1, tmp1$n.batches[j])))
    } else {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.hat )) / (sqrt(tmp4$delta.hat[j,] ) %*% t(rep(1, tmp1$n.batches[j])))
    }
    j <- j + 1
  }
  return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + tmp3$stand.mean)
}
#' Fit a ComBat Harmonization Model
#'
#' @description
#' Fits a ComBat batch harmonization model to a numeric matrix, with optional
#' covariates. This implementation is adapted from the \code{sva} package with
#' MRI-focused modifications and robustness patches. It supports missing values
#' (imputation during fit), constant-feature handling, and optional empirical
#' Bayes (EB) shrinkage.
#'
#' @param dat Numeric matrix or data frame. Either rows or columns must match
#'   the length of \code{batch}. If subjects are rows, the data are transposed
#'   internally. Features are assumed to be rows after any transpose.
#' @param batch Factor indicating batch/site per subject. Must be a factor; unused
#'   levels are dropped with a warning.
#' @param mod Optional numeric design matrix of covariates (no intercept required).
#'   Must have the same number of rows as \code{length(batch)} when provided.
#' @param eb Logical. If \code{TRUE}, use empirical Bayes shrinkage (default).
#' @param verbose Logical. Print progress messages.
#' @param single_batch Character. What to do if there is only one batch level:
#'   \code{"error"} (default) stops; \code{"return"} returns a no-op model that
#'   will pass data through unchanged in \code{\link{combat_apply}}.
#'
#' @return A list containing the fitted model components, including:
#' \itemize{
#'   \item \code{levels_batch}: factor levels used during fit
#'   \item \code{transpose}: whether the input was transposed
#'   \item \code{not_constant}: indices of non-constant features
#'   \item \code{eb}: whether EB was used
#'   \item \code{tmp2}, \code{tmp4}: internal fitted statistics
#'   \item \code{mod_info}: metadata for validating apply-time covariates
#' }
#'
#' @details
#' Missing values are imputed during fit. If \code{mod} is \code{NULL}, missing
#' values are imputed using within-batch means with a global-mean fallback; a
#' warning is emitted to make this behavior explicit. If \code{mod} is provided,
#' a regression-based imputation is used (mirroring the original ComBat logic).
#'
#' Features with zero (or undefined) variance are filtered out prior to fitting.
#'
#' @seealso \code{\link{combat_apply}} for applying a fitted model.
#' @export
combat_fit <- function (dat, batch, mod = NULL, eb = TRUE, verbose = TRUE,
                        single_batch = c("error", "return")) {
  # Modified by Joaquim Radua
  single_batch <- match.arg(single_batch)
  if (!is.factor(batch)) {
    stop("batch must be a factor")
  }
  dat <- as.matrix(dat)
  if (ncol(dat) == length(batch)) {
    transpose <- FALSE
    if (verbose) {
      cat("[combat] Subjects are COLUMNS\n")
    }
  } else if (nrow(dat) == length(batch)) {
    transpose <- TRUE
    dat <- t(dat)
    if (verbose) {
      cat("[combat] Subjects are ROWS\n")
    }
  } else {
    stop("dat must have the same number of columns or rows than the length of batch")
  }
  if (is.data.frame((mod))) {
    mod <- as.matrix(mod)
  } else if (!(is.matrix(mod) || is.null(mod))) {
    stop("mod must be a matrix or NULL")
  }
  if (!is.null(mod) && is.null(colnames(mod))) {
    colnames(mod) <- paste0("V", seq_len(ncol(mod)))
  }
  if (!is.null(mod) && nrow(mod) != length(batch)) {
    stop("mod must have the same number of rows as the length of batch")
  }
  if (any(table(batch) == 0)) {
    warning("W_BATCH_EMPTY_DROPPED: Dropping unused batch levels.", call. = FALSE)
    batch <- droplevels(batch)
  }
  if (nlevels(batch) < 2) {
    if (single_batch == "return") {
      warning("W_SINGLE_BATCH_NOOP: ComBat requires >= 2 batches; returning input unchanged.", call. = FALSE)
      levels_batch <- levels(batch)
      return(list(
        single_batch = TRUE,
        levels_batch = levels_batch,
        transpose = transpose,
        not_constant = seq_len(nrow(dat)),
        eb = eb,
        tmp2 = NULL,
        tmp4 = NULL,
        mod_info = list(
          has_mod = !is.null(mod),
          ncol = if (is.null(mod)) 0L else ncol(mod),
          colnames = colnames(mod),
          design_cols = NULL,
          design_ncol = 0L
        )
      ))
    }
    stop("ComBat requires >= 2 batches; got 1.")
  }
  # Imputation of missing values
  if (any(is.na(dat))) {
    if (verbose) {
      cat("[combat] Imputing missing data (only for fit)\n")
    }
    if (is.null(mod)) {
      warning(
        "W_IMPUTE_MEAN_NO_MOD: Missing values imputed with within-batch mean (global mean fallback) because mod is NULL.",
        call. = FALSE
      )
      for (batch_i in sort(unique(batch))) {
        i <- which(batch == batch_i)
        if (length(i) == 0) next
        dat_i <- dat[,i, drop = FALSE]
        if (any(is.na(dat_i))) {
          means <- rowMeans(dat_i, na.rm = TRUE)
          na_mat <- is.na(dat_i)
          if (any(na_mat)) {
            dat_i[na_mat] <- means[row(na_mat)]
            dat[,i] <- dat_i
          }
        }
      }
      if (any(is.na(dat))) {
        global_means <- rowMeans(dat, na.rm = TRUE)
        if (any(!is.finite(global_means))) {
          warning(
            "W_IMPUTE_ALL_NA_FEATURE: Some features are all-NA; imputing with 0 for global mean fallback.",
            call. = FALSE
          )
          global_means[!is.finite(global_means)] <- 0
        }
        na_mat <- is.na(dat)
        if (any(na_mat)) {
          dat[na_mat] <- global_means[row(na_mat)]
        }
      }
    } else {
      for (batch_i in sort(unique(batch))) {
        i <- which(batch == batch_i)
        dat_i <- dat[,i]
        if (any(is.na(dat_i))) {
          for (j in 1:nrow(dat)) {
            dat_ji <- dat_i[j,]
            is_na <- which(is.na(dat_ji))
            # Some missing, impute from other subjects of the site
            if (length(is_na) > 0 && length(is_na) < length(i)) {
              if (length(is_na) == 1) {
                mod_i_is_na <- matrix(mod[i[is_na],], nrow = 1)
              } else {
                mod_i_is_na <- mod[i[is_na],]
              }
              beta <- matrix(coef(lm(dat_ji ~ mod[i,])))
              beta[which(is.na(beta))] <- 0
              dat[j,i[is_na]] <- cbind(1, mod_i_is_na) %*% beta
            } else {
              dat[j,i[is_na]] <- mean(dat_ji, na.rm = TRUE)
            }
          }
        }
      }
      # If there are still missing data, they are because a site has all missing
      # data for a ROI
      for (batch_i in sort(unique(batch))) {
        i <- which(batch == batch_i)
        if (any(is.na(dat[,i]))) {
          for (j in 1:nrow(dat)) {
            dat_j <- dat[j,]
            if (is.na(dat_j[i[1]])) {
              if (!is.null(mod)) {
                beta <- matrix(coef(lm(dat_j ~ mod)))
                beta[which(is.na(beta))] <- 0
                dat[j,i] <- cbind(1, mod[i,]) %*% beta
              } else {
                dat[j,i] <- mean(dat_j, na.rm = TRUE)
              }
            }
          }
        }
      }
    }
  }
  not_constant <- which(apply(dat, 1, function (x) {var(x, na.rm = TRUE) > 0}))
  dat <- dat[not_constant, , drop = FALSE]
  if (eb) {
    if (verbose) {
      cat("[combat] Performing ComBat with empirical Bayes\n")
    }
  } else {
    if (verbose) {
      cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
    }
  }
  if (verbose) {
    cat("[combat] Found", nlevels(batch), "batches (e.g., sites)\n")
  }
  levels_batch <- levels(batch)
  tmp1 <- .combat_tmp1(dat, batch, levels_batch, mod)
  tmp2 <- .combat_tmp2(tmp1, verbose)
  tmp3 <- .combat_tmp3(dat, tmp1, tmp2, verbose)
  tmp4 <- .combat_tmp4(tmp1, tmp2, tmp3, eb, verbose)
  return(list(
    levels_batch = levels_batch,
    transpose = transpose,
    not_constant = not_constant,
    eb = eb,
    tmp2 = tmp2,
    tmp4 = tmp4,
    mod_info = list(
      has_mod = !is.null(mod),
      ncol = if (is.null(mod)) 0L else ncol(mod),
      colnames = colnames(mod),
      design_cols = tmp1$design_cols,
      design_ncol = ncol(tmp1$design)
    )
  ))
}
#' Apply a ComBat Harmonization Model
#'
#' @description
#' Applies a fitted ComBat model (from \code{\link{combat_fit}}) to new data with
#' matching batch levels and covariate design. Supports no-op models returned
#' by \code{single_batch = "return"}.
#'
#' @param tmp Fitted model list returned by \code{\link{combat_fit}}.
#' @param dat Numeric matrix or data frame. Either rows or columns must match
#'   the length of \code{batch}. If the fitted model stored a transpose flag,
#'   the same orientation is enforced.
#' @param batch Factor indicating batch/site per subject. Levels must match the
#'   levels used during fit.
#' @param mod Optional numeric design matrix of covariates. Must match the
#'   fit-time design (column count and names); missing fit-time columns are
#'   zero-padded, but new columns are rejected.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with harmonized data and fitted parameters used in adjustment:
#' \itemize{
#'   \item \code{dat.combat}: adjusted data matrix in the original orientation
#'   \item batch/EB parameter estimates and priors for auditing
#' }
#'
#' @export
combat_apply <- function (tmp, dat, batch, mod = NULL, verbose = TRUE) {
  # Modified by Joaquim Radua
  if (!is.factor(batch)) {
    stop("batch must be a factor")
  }
  batch <- factor(batch, levels = tmp$levels_batch)
  if (any(is.na(batch))) {
    stop("batch contains levels not seen during fit")
  }
  if (any(table(batch) == 0)) {
    warning("W_BATCH_MISSING_IN_APPLY: Some batch levels from fit are absent in apply data; proceeding.", call. = FALSE)
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  } else if (!is.matrix(dat)) {
    stop("dat must be a matrix")
  }
  if (tmp$transpose) {
    if (nrow(dat) != length(batch)) {
      stop("dat must have the same number of rows than the length of batch")
    }
    dat <- t(dat)
  } else {
    if (ncol(dat) != length(batch)) {
      stop("dat must have the same number of columns than the length of batch")
    }
  }
  if (is.null(dim(dat))) {
    dat <- matrix(dat, nrow = 1)
  }
  dat <- as.matrix(dat)
  if (isTRUE(tmp$single_batch)) {
    if (verbose) {
      cat("[combat] Single-batch model; returning input unchanged\n")
    }
    warning("W_SINGLE_BATCH_NOOP: ComBat requires >= 2 batches; returning input unchanged.", call. = FALSE)
    dat.combat <- dat
    if (tmp$transpose) {
      dat.combat <- t(dat.combat)
    }
    return(list(
      dat.combat = dat.combat, 
      gamma.hat = NULL,
      delta.hat = NULL, 
      gamma.star = NULL,
      delta.star = NULL, 
      gamma.bar = NULL,
      t2 = NULL,
      a.prior = NULL,
      b.prior = NULL,
      batch = batch,
      mod = mod, 
      stand.mean = NULL,
      stand.sd = NULL
    ))
  }
  dat.combat <- as.matrix(dat)
  dat <- dat[tmp$not_constant, , drop = FALSE]
  if (is.data.frame((mod))) {
    mod <- as.matrix(mod)
  } else if (!(is.matrix(mod) || is.null(mod))) {
    stop("mod must be a matrix or NULL")
  }
  if (!is.null(mod) && is.null(colnames(mod))) {
    colnames(mod) <- paste0("V", seq_len(ncol(mod)))
  }
  if (!is.null(mod) && nrow(mod) != length(batch)) {
    stop("mod must have the same number of rows as the length of batch")
  }
  if (!is.null(tmp$mod_info)) {
    if (isTRUE(tmp$mod_info$has_mod)) {
      if (is.null(mod)) {
        stop("mod must be provided at apply time because it was used during fit")
      }
    } else if (!is.null(mod)) {
      stop("mod must be NULL at apply time because it was NULL during fit")
    }
  }
  tmp1 <- .combat_tmp1(dat, batch, tmp$levels_batch, mod)
  if (!is.null(tmp$mod_info)) {
    if (!is.null(tmp$mod_info$design_cols)) {
      fit_cols <- tmp$mod_info$design_cols
      apply_cols <- colnames(tmp1$design)
      extra_cols <- setdiff(apply_cols, fit_cols)
      if (length(extra_cols) > 0) {
        stop("Design matrix at apply time has unexpected columns not seen in fit.")
      }
      missing_cols <- setdiff(fit_cols, apply_cols)
      if (length(missing_cols) > 0) {
        add <- matrix(0, nrow = nrow(tmp1$design), ncol = length(missing_cols),
                      dimnames = list(NULL, missing_cols))
        tmp1$design <- cbind(tmp1$design, add)
      }
      tmp1$design <- tmp1$design[, fit_cols, drop = FALSE]
      tmp1$design_cols <- fit_cols
      tmp1$batch.col.idx <- match(colnames(tmp1$batchmod), tmp1$design_cols)
      if (any(is.na(tmp1$batch.col.idx))) {
        stop("Batch columns missing after design alignment.")
      }
      tmp1$batch.design <- tmp1$design[, tmp1$batch.col.idx, drop = FALSE]
    }
  }
  tmp3 <- .combat_tmp3(dat, tmp1, tmp$tmp2, verbose)
  tmp5 <- .combat_tmp5(tmp1, tmp$tmp2, tmp3, tmp$tmp4, tmp$eb, verbose)
  dat.combat[tmp$not_constant, ] <- tmp5
  if (tmp$transpose) {
    dat.combat <- t(dat.combat)
  }
  return(list(
    dat.combat = dat.combat, 
    gamma.hat = tmp$tmp4$gamma.hat,
    delta.hat = tmp$tmp4$delta.hat, 
    gamma.star = tmp$tmp4$gamma.star,
    delta.star = tmp$tmp4$delta.star, 
    gamma.bar = tmp$tmp4$gamma.bar,
    t2 = tmp$tmp4$t2,
    a.prior = tmp$tmp4$a.prior,
    b.prior = tmp$tmp4$b.prior,
    batch = tmp1$batch,
    mod = tmp1$mod, 
    stand.mean = tmp3$stand.mean,
    stand.sd = sqrt(tmp$tmp2$var.pooled)[,1]
  ))
}


################################################################################
#                                                                              #
# This is a copy of the original code from the standard version of the sva     #
# package that can be found at                                                 #
# https://bioconductor.org/packages/release/bioc/html/sva.html                 #
# The original and present code is under the Artistic License 2.0.             #
# If using this code, make sure you agree and accept this license.             #
#                                                                              #
################################################################################


# Following four find empirical hyper-prior values
.aprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2 * s2 + m^2) / s2
}
.bprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m * s2 + m^3) / s2
}
.postmean <- function (g.hat, g.bar, n, d.star, t2) {
  (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
}
.postvar <- function (sum2, n, a, b) {
  (0.5 * sum2 + b) / (n / 2 + a - 1)
}
# Pass in entire data set, the design matrix for the entire data, the batch
# means, the batch variances, priors (m, t2, a, b), columns of the data 
# matrix for the batch. Uses the EM to find the parametric batch adjustments
.it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) {
  sdat <- as.matrix(sdat)
  n <- apply(!is.na(sdat), 1, sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- .postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 1, sum, na.rm = T)
    d.new <- .postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old) / (abs(g.old) + 1e-10),
                  abs(d.new - d.old) / (abs(d.old) + 1e-10))
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  # cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}
