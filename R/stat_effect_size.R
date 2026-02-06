# Effect size utilities

.hedges_g <- function(x_ref, x_target) {
  n1 <- sum(!is.na(x_ref)); n2 <- sum(!is.na(x_target))
  m1 <- mean(x_ref, na.rm = TRUE); m2 <- mean(x_target, na.rm = TRUE)
  s1 <- stats::sd(x_ref, na.rm = TRUE); s2 <- stats::sd(x_target, na.rm = TRUE)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  d <- (m2 - m1) / sp
  j <- 1 - (3 / (4 * (n1 + n2) - 9))
  g <- d * j
  list(es = g, es_type = "hedges_g_pooled")
}

.cohen_d <- function(x_ref, x_target) {
  n1 <- sum(!is.na(x_ref)); n2 <- sum(!is.na(x_target))
  m1 <- mean(x_ref, na.rm = TRUE); m2 <- mean(x_target, na.rm = TRUE)
  s1 <- stats::sd(x_ref, na.rm = TRUE); s2 <- stats::sd(x_target, na.rm = TRUE)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  d <- (m2 - m1) / sp
  list(es = d, es_type = "cohen_d_pooled")
}

.glass_delta <- function(x_ref, x_target) {
  m1 <- mean(x_ref, na.rm = TRUE); m2 <- mean(x_target, na.rm = TRUE)
  s1 <- stats::sd(x_ref, na.rm = TRUE)
  d <- (m2 - m1) / s1
  list(es = d, es_type = "glass_delta")
}

.effect_size_two_group <- function(x, g, ref, target, type = "hedges_g_pooled",
                                   ci = FALSE) {
  g <- droplevels(factor(g))
  if (nlevels(g) < 2) {
    return(list(es = NA_real_, es_type = type, es_notes = "not_two_group",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }
  if (is.null(ref) || is.null(target)) {
    return(list(es = NA_real_, es_type = type, es_notes = "missing_ref_or_target",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }
  if (!ref %in% levels(g) || !target %in% levels(g)) {
    return(list(es = NA_real_, es_type = type, es_notes = "ref_or_target_not_found",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }
  x_ref <- x[g == ref]
  x_target <- x[g == target]

  res <- switch(
    type,
    hedges_g_pooled = .hedges_g(x_ref, x_target),
    cohen_d_pooled = .cohen_d(x_ref, x_target),
    glass_delta = .glass_delta(x_ref, x_target),
    list(es = NA_real_, es_type = "none")
  )
  res$es_notes <- NULL

  if (isTRUE(ci) && is.finite(res$es)) {
    n1 <- sum(!is.na(x_ref)); n2 <- sum(!is.na(x_target))
    se <- sqrt((n1 + n2) / (n1 * n2) + (res$es^2) / (2 * (n1 + n2 - 2)))
    res$es_ci_low <- res$es - 1.96 * se
    res$es_ci_high <- res$es + 1.96 * se
  } else {
    res$es_ci_low <- NA_real_
    res$es_ci_high <- NA_real_
  }
  res
}

.effect_size_two_group_model <- function(beta, model, g, ref, target,
                                         type = "hedges_g_pooled", ci = FALSE) {
  g <- droplevels(factor(g))
  if (nlevels(g) < 2) {
    return(list(es = NA_real_, es_type = type, es_notes = "not_two_group",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }
  if (is.null(ref) || is.null(target) || !ref %in% levels(g) || !target %in% levels(g)) {
    return(list(es = NA_real_, es_type = type, es_notes = "bad_ref_or_target",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }
  if (is.null(model) || !is.finite(beta)) {
    return(list(es = NA_real_, es_type = type, es_notes = "missing_model_or_beta",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }

  # Use within-group residual SD as denominator (covariate-adjusted scale).
  r <- tryCatch(stats::residuals(model), error = function(e) NULL)
  if (is.null(r)) {
    return(list(es = NA_real_, es_type = type, es_notes = "no_residuals",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }

  r_ref <- r[g == ref]
  r_target <- r[g == target]
  n1 <- sum(!is.na(r_ref)); n2 <- sum(!is.na(r_target))
  if (n1 < 2 || n2 < 2) {
    return(list(es = NA_real_, es_type = type, es_notes = "too_few_per_group",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }

  s1 <- stats::sd(r_ref, na.rm = TRUE)
  s2 <- stats::sd(r_target, na.rm = TRUE)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp <= 0) {
    return(list(es = NA_real_, es_type = type, es_notes = "bad_pooled_sd",
                es_ci_low = NA_real_, es_ci_high = NA_real_))
  }

  d <- beta / sp
  res <- switch(
    type,
    hedges_g_pooled = {
      j <- 1 - (3 / (4 * (n1 + n2) - 9))
      list(es = d * j, es_type = "hedges_g_pooled")
    },
    cohen_d_pooled = list(es = d, es_type = "cohen_d_pooled"),
    glass_delta = list(es = beta / stats::sd(r_ref, na.rm = TRUE), es_type = "glass_delta"),
    list(es = NA_real_, es_type = "none")
  )
  res$es_notes <- "model_residual_sd"

  if (isTRUE(ci) && is.finite(res$es)) {
    se <- sqrt((n1 + n2) / (n1 * n2) + (res$es^2) / (2 * (n1 + n2 - 2)))
    res$es_ci_low <- res$es - 1.96 * se
    res$es_ci_high <- res$es + 1.96 * se
  } else {
    res$es_ci_low <- NA_real_
    res$es_ci_high <- NA_real_
  }

  res
}

.standardized_beta <- function(model, data, term) {
  if (is.null(model)) {
    return(list(es = NA_real_, es_type = "standardized_beta"))
  }
  if (!term %in% names(stats::coef(model))) {
    return(list(es = NA_real_, es_type = "standardized_beta"))
  }
  y <- stats::model.response(stats::model.frame(model))
  x <- data[[term]]
  if (is.null(x)) return(list(es = NA_real_, es_type = "standardized_beta"))
  b <- stats::coef(model)[[term]]
  es <- b * stats::sd(x, na.rm = TRUE) / stats::sd(y, na.rm = TRUE)
  list(es = es, es_type = "standardized_beta")
}
