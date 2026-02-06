# Diagnostics helpers

.diagnostics_lm <- function(model, data, term = NA_character_, level = "light") {
  n <- stats::nobs(model)
  res <- list(
    term = term,
    n = n,
    r2 = summary(model)$r.squared,
    r2_adj = summary(model)$adj.r.squared,
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    sigma = summary(model)$sigma,
    shapiro_p = NA_real_,
    bp_p = NA_real_,
    cooks_max = max(stats::cooks.distance(model), na.rm = TRUE),
    cooks_n_high = sum(stats::cooks.distance(model) > (4 / n), na.rm = TRUE),
    vif_max = NA_real_
  )

  if (level %in% c("light", "full")) {
    if (n <= 5000) {
      res$shapiro_p <- tryCatch(stats::shapiro.test(stats::resid(model))$p.value, error = function(e) NA_real_)
    }
    if (requireNamespace("lmtest", quietly = TRUE)) {
      res$bp_p <- tryCatch(lmtest::bptest(model)$p.value, error = function(e) NA_real_)
    }
    if (requireNamespace("car", quietly = TRUE)) {
      v <- tryCatch(car::vif(model), error = function(e) NULL)
      if (!is.null(v)) res$vif_max <- max(v, na.rm = TRUE)
    }
  }
  as.data.frame(res)
}

.diagnostics_glm <- function(model, data, term = NA_character_, level = "light") {
  n <- stats::nobs(model)
  res <- list(
    term = term,
    n = n,
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    deviance = model$deviance,
    df_residual = model$df.residual,
    dispersion = summary(model)$dispersion,
    pseudo_r2 = NA_real_,
    convergence = isTRUE(model$converged)
  )
  if (level %in% c("light", "full")) {
    if (requireNamespace("pscl", quietly = TRUE)) {
      pr2 <- tryCatch(pscl::pR2(model), error = function(e) NULL)
      if (!is.null(pr2)) res$pseudo_r2 <- pr2["McFadden"]
    }
  }
  as.data.frame(res)
}

.diagnostics_lmer <- function(model, data, term = NA_character_, level = "light") {
  .check_pkg("lme4", "mixed model diagnostics")
  n <- stats::nobs(model)
  res <- list(
    term = term,
    n = n,
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)),
    singular = lme4::isSingular(model, tol = 1e-4),
    convergence = is.null(model@optinfo$conv$lme4$messages)
  )
  if (level %in% c("light", "full") && requireNamespace("performance", quietly = TRUE)) {
    r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
    if (!is.null(r2)) {
      res$r2_marginal <- r2$R2_marginal
      res$r2_conditional <- r2$R2_conditional
    } else {
      res$r2_marginal <- NA_real_
      res$r2_conditional <- NA_real_
    }
  }
  as.data.frame(res)
}

.diagnostics_glmer <- function(model, data, term = NA_character_, level = "light") {
  .check_pkg("lme4", "mixed model diagnostics")
  n <- stats::nobs(model)
  res <- list(
    term = term,
    n = n,
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)),
    singular = lme4::isSingular(model, tol = 1e-4),
    convergence = is.null(model@optinfo$conv$lme4$messages)
  )
  if (level %in% c("light", "full") && requireNamespace("performance", quietly = TRUE)) {
    r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
    if (!is.null(r2)) {
      res$r2_marginal <- r2$R2_marginal
      res$r2_conditional <- r2$R2_conditional
    } else {
      res$r2_marginal <- NA_real_
      res$r2_conditional <- NA_real_
    }
  }
  as.data.frame(res)
}

