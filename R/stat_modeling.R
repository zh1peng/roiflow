# Core modeling utilities and public analysis functions

.normalize_family <- function(family) {
  if (is.null(family)) return(NULL)
  if (inherits(family, "family")) return(family)
  if (is.character(family)) {
    fam_fun <- get(family, mode = "function", envir = asNamespace("stats"))
    return(fam_fun())
  }
  stop("family must be NULL, a family object, or a character name.", call. = FALSE)
}

.choose_engine <- function(spec) {
  if (!identical(spec$model_engine, "auto")) return(spec$model_engine)
  if (!is.null(spec$random_effects)) {
    if (!is.null(spec$family) && spec$family$family == "binomial") return("glmer")
    return("lmer")
  }
  if (!is.null(spec$family)) return("glm")
  "lm"
}

.build_formula <- function(y, group_col = NULL, covariates = NULL, random_effects = NULL) {
  terms <- c()
  if (!is.null(group_col)) terms <- c(terms, .quote_var(group_col))
  if (!is.null(covariates) && length(covariates) > 0) {
    terms <- c(terms, vapply(covariates, .quote_var, character(1)))
  }
  if (!is.null(random_effects)) {
    terms <- c(terms, random_effects)
  }
  rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
  stats::as.formula(paste(.quote_var(y), "~", rhs))
}

.fit_model <- function(formula, data, spec) {
  engine <- .choose_engine(spec)
  family <- .normalize_family(spec$family)

  fit_call <- switch(
    engine,
    lm = quote(stats::lm(formula, data = data)),
    glm = quote(stats::glm(formula, data = data, family = family)),
    lmer = {
      .check_pkg("lme4", "lmer models")
      quote(lme4::lmer(formula, data = data))
    },
    glmer = {
      .check_pkg("lme4", "glmer models")
      quote(lme4::glmer(formula, data = data, family = family))
    },
    stop("Unsupported model_engine.", call. = FALSE)
  )

  cap <- .capture_warnings(eval(fit_call))
  list(model = cap$value, engine = engine, warnings = cap$warnings, family = family)
}

.tidy_model <- function(model, engine, conf_int = FALSE, conf_level = 0.95) {
  if (engine %in% c("lm", "glm")) {
    .check_pkg("broom", "tidy outputs")
    td <- broom::tidy(model, conf.int = conf_int, conf.level = conf_level)
  } else {
    .check_pkg("broom.mixed", "tidy mixed-model outputs")
    td <- broom.mixed::tidy(model, conf.int = conf_int, conf.level = conf_level, effects = "fixed")
  }
  td
}

.standardize_tidy <- function(td) {
  names(td) <- gsub("std.error", "std_error", names(td), fixed = TRUE)
  names(td) <- gsub("p.value", "p_value", names(td), fixed = TRUE)
  names(td) <- gsub("conf.low", "conf_low", names(td), fixed = TRUE)
  names(td) <- gsub("conf.high", "conf_high", names(td), fixed = TRUE)
  if (!"std_error" %in% names(td)) td$std_error <- NA_real_
  if (!"statistic" %in% names(td)) td$statistic <- NA_real_
  if (!"p_value" %in% names(td)) td$p_value <- NA_real_
  if (!"df" %in% names(td)) td$df <- NA_real_
  if (!"conf_low" %in% names(td)) td$conf_low <- NA_real_
  if (!"conf_high" %in% names(td)) td$conf_high <- NA_real_
  td
}

.term_filter <- function(td, pattern = NULL, terms = NULL) {
  if (!is.null(terms)) {
    return(td[td$term %in% terms, , drop = FALSE])
  }
  if (!is.null(pattern)) {
    return(td[grepl(pattern, td$term), , drop = FALSE])
  }
  td
}

.term_to_level <- function(term, group_col) {
  if (!is.character(term) || length(term) == 0) return(NA_character_)
  if (startsWith(term, group_col)) {
    sub(paste0("^", group_col), "", term)
  } else {
    NA_character_
  }
}

.resolve_group_ref <- function(spec, levels_vec) {
  ref <- spec$contrast$ref %||% spec$ref_levels[[spec$group_col]]
  if (is.null(ref)) {
    ref <- .infer_ref_level(levels_vec)
  }
  if (is.null(ref)) {
    stop("E_REF_NOT_SET: No reference level specified for group column.", call. = FALSE)
  }
  ref
}

.resolve_group_targets <- function(spec, levels_vec, ref) {
  target <- spec$contrast$target %||% NULL
  if (!is.null(target)) {
    if (!target %in% levels_vec) {
      stop("E_TARGET_NOT_FOUND: Target level not found in group levels.", call. = FALSE)
    }
    return(target)
  }
  setdiff(levels_vec, ref)
}

.check_sign_alignment <- function(beta, es, var, term) {
  if (is.na(beta) || is.na(es)) return(invisible(TRUE))
  if (beta == 0 || es == 0) return(invisible(TRUE))
  if (sign(beta) != sign(es)) {
    stop(sprintf("Effect size direction mismatch: expected target-ref but computed ref-target (var=%s, term=%s).",
                 var, term), call. = FALSE)
  }
  invisible(TRUE)
}

#' Compare Groups Across Variables
#'
#' @description
#' Fits group comparison models across many variables (e.g., ROIs) and returns
#' a standardized results table.
#'
#' @param df Data frame with outcomes, group column, and covariates.
#' @param vars Character vector of outcome variables.
#' @param spec Analysis specification from \code{\link{analysis_spec}}.
#' @param return Character. \code{"result"} (default) returns a structured result
#'   object; \code{"data"} returns only the results table.
#'
#' @return A \code{profile_result} object or results data frame.
#' @export
compare_groups <- function(df, vars, spec = analysis_spec(), return = c("result", "data")) {
  return <- match.arg(return)
  .ensure_cols(df, c(vars, spec$group_col))
  if (!is.null(spec$contrast$direction) &&
      !identical(spec$contrast$direction, "target_minus_ref")) {
    stop("Only contrast direction 'target_minus_ref' is supported.", call. = FALSE)
  }

  if (!is.null(spec$site_col) && spec$site_effect != "none") {
    if (spec$site_effect == "fixed") {
      spec$covariates <- unique(c(spec$covariates, spec$site_col))
    } else if (spec$site_effect == "random") {
      spec$random_effects <- sprintf("(1|%s)", spec$site_col)
    }
  }

  res_list <- list()
  diag_list <- list()
  log <- list()

  for (v in vars) {
    cols <- unique(c(v, spec$group_col, spec$covariates))
    data_v <- df[, intersect(cols, names(df)), drop = FALSE]
    keep <- .complete_cases(data_v)
    n_dropped <- sum(!keep)
    if (n_dropped > 0) {
      log <- .log_add(log, "W_DROPPED_ROWS_COMPLETE_CASE",
                      "Rows dropped due to missing data.", list(var = v, n_dropped = n_dropped))
    }
    data_v <- data_v[keep, , drop = FALSE]
    if (nrow(data_v) == 0) next

    data_v <- .apply_ref_levels(data_v, spec$ref_levels)

    group_raw <- data_v[[spec$group_col]]
    group_levels_raw <- levels(factor(group_raw))
    ref_level <- .resolve_group_ref(spec, group_levels_raw)
    data_v[[spec$group_col]] <- .set_ref_factor(group_raw, ref_level, spec$group_col)
    group <- data_v[[spec$group_col]]
    group_levels <- levels(group)
    group_counts <- table(group)
    target_levels <- .resolve_group_targets(spec, group_levels, ref_level)

    f <- .build_formula(v, spec$group_col, spec$covariates, spec$random_effects)
    fit <- .fit_model(f, data_v, spec)

    if (length(fit$warnings) > 0) {
      for (w in fit$warnings) {
        log <- .log_add(log, .classify_warning(w), w, list(var = v))
      }
    }

    td <- .standardize_tidy(.tidy_model(fit$model, fit$engine, spec$conf_int, spec$conf_level))
    td <- .term_filter(td, pattern = paste0("^", spec$group_col))
    if (nrow(td) == 0) next

    td$var <- v
    td$model_engine <- fit$engine
    td$family <- if (is.null(fit$family)) NA_character_ else fit$family$family
    td$formula_text <- paste(deparse(f), collapse = "")
    td$n <- nrow(data_v)
    td$n_groups <- length(group_levels)
    td$group_levels <- paste(group_levels, collapse = "|")
    td$group_counts <- paste(paste(group_levels, as.integer(group_counts), sep = ":"), collapse = "|")
    td$ref_level <- ref_level
    td$target_level <- vapply(td$term, .term_to_level, character(1), group_col = spec$group_col)
    td$target_level[is.na(td$target_level)] <- if (length(target_levels) == 1) target_levels else NA_character_
    td$contrast_label <- ifelse(!is.na(td$target_level),
                                paste(td$target_level, "vs", td$ref_level),
                                NA_character_)

    es_ok <- spec$effect_size_scope != "none" &&
      spec$effect_size %in% c("hedges_g_pooled", "cohen_d_pooled", "glass_delta")
    if (spec$effect_size_scope == "two_group_only" && length(group_levels) != 2) {
      es_ok <- FALSE
    }
    if (es_ok) {
      covars <- spec$covariates %||% character(0)
      from <- spec$effect_size_from %||% "raw_data"
      es_vals <- lapply(seq_len(nrow(td)), function(i) {
        tg <- td$target_level[i]
        if (is.na(tg)) return(list(es = NA_real_, es_type = spec$effect_size,
                                   es_notes = "target_missing", es_ci_low = NA_real_, es_ci_high = NA_real_))
        if (identical(from, "model")) {
          .effect_size_two_group_model(beta = td$estimate[i], model = fit$model,
                                       g = group, ref = ref_level, target = tg,
                                       type = spec$effect_size, ci = (spec$effect_size_ci != "none"))
        } else {
          .effect_size_two_group(data_v[[v]], group, ref_level, tg,
                                 type = spec$effect_size, ci = (spec$effect_size_ci != "none"))
        }
      })
      td$es_value <- vapply(es_vals, function(x) x$es, numeric(1))
      td$es_method <- vapply(es_vals, function(x) x$es_type, character(1))
      td$es_notes <- vapply(es_vals, function(x) x$es_notes %||% NA_character_, character(1))
      td$es_ci_low <- vapply(es_vals, function(x) x$es_ci_low, numeric(1))
      td$es_ci_high <- vapply(es_vals, function(x) x$es_ci_high, numeric(1))

      if (identical(from, "raw_data") && length(covars) > 0) {
        td$es_notes[is.na(td$es_notes)] <- "ignores_covariates"
      }
    } else {
      td$es_value <- NA_real_
      td$es_method <- spec$effect_size
      td$es_notes <- NA_character_
      td$es_ci_low <- NA_real_
      td$es_ci_high <- NA_real_
    }
    td$es <- td$es_value
    td$es_type <- td$es_method

    for (i in seq_len(nrow(td))) {
      if (is.na(td$es_value[i]) || td$es_value[i] == 0) next
      from <- spec$effect_size_from %||% "raw_data"
      if (identical(from, "model")) {
        .check_sign_alignment(td$estimate[i], td$es_value[i], v, td$term[i])
      } else {
        tg <- td$target_level[i]
        if (is.na(tg)) next
        x_ref <- data_v[[v]][group == ref_level]
        x_tg <- data_v[[v]][group == tg]
        diff_raw <- mean(x_tg, na.rm = TRUE) - mean(x_ref, na.rm = TRUE)
        .check_sign_alignment(diff_raw, td$es_value[i], v, td$term[i])

        covars <- spec$covariates %||% character(0)
        if (length(covars) > 0 &&
            is.finite(td$estimate[i]) && td$estimate[i] != 0 &&
            sign(td$estimate[i]) != sign(td$es_value[i])) {
          log <- .log_add(
            log,
            "W_ES_SIGN_DIFF_COVAR",
            "Effect size sign differs from covariate-adjusted beta (effect_size_from='raw_data').",
            list(var = v, term = td$term[i], ref = ref_level, target = tg)
          )
        }
      }
    }

    if (spec$diagnostics != "none") {
      diag_fun <- switch(
        fit$engine,
        lm = .diagnostics_lm,
        glm = .diagnostics_glm,
        lmer = .diagnostics_lmer,
        glmer = .diagnostics_glmer
      )
      diag_list[[v]] <- diag_fun(fit$model, data_v, term = spec$group_col, level = spec$diagnostics)
      diag_list[[v]]$var <- v
    }

    res_list[[v]] <- td
  }

  res <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()
  res <- .apply_p_adjust(res, p_col = "p_value", spec = spec)
  summaries <- tryCatch(.summarize_groups(df, vars, spec$group_col), error = function(e) NULL)
  diagnostics <- if (length(diag_list) > 0) do.call(rbind, diag_list) else NULL

  out <- .result_object(
    data = res,
    diagnostics = diagnostics,
    summaries = summaries,
    meta = list(spec = spec, call = match.call()),
    log = log
  )
  if (return == "data") return(out$data)
  out
}

#' Pairwise Group Comparisons
#'
#' @description
#' Performs pairwise group comparisons across all group levels using the same
#' backend and schema as \code{\link{compare_groups}}.
#'
#' @param df Data frame with outcomes and group column.
#' @param vars Character vector of outcome variables.
#' @param spec Analysis specification from \code{\link{analysis_spec}}.
#' @param return Character. \code{"result"} (default) or \code{"data"}.
#'
#' @return A \code{profile_result} object or results data frame.
#' @export
pairwise_groups <- function(df, vars, spec = analysis_spec(), return = c("result", "data")) {
  return <- match.arg(return)
  .ensure_cols(df, c(vars, spec$group_col))
  group <- droplevels(factor(df[[spec$group_col]]))
  levels_g <- levels(group)
  if (length(levels_g) < 2) stop("Need at least two groups for pairwise comparisons.", call. = FALSE)

  comb <- utils::combn(levels_g, 2, simplify = FALSE)
  res_list <- list()
  diag_list <- list()
  log <- list()

  for (pair in comb) {
    df_pair <- df[group %in% pair, , drop = FALSE]
    df_pair[[spec$group_col]] <- droplevels(factor(df_pair[[spec$group_col]]))
    spec_pair <- spec
    spec_pair$contrast$ref <- pair[1]
    spec_pair$contrast$target <- pair[2]
    spec_pair$ref_levels[[spec$group_col]] <- pair[1]
    res_pair <- compare_groups(df_pair, vars, spec_pair, return = "result")
    if (nrow(res_pair$data) == 0) next
    res_pair$data$comparison <- paste(pair[1], "vs", pair[2])
    res_list[[paste(pair, collapse = "_")]] <- res_pair$data
    if (!is.null(res_pair$diagnostics)) {
      res_pair$diagnostics$comparison <- paste(pair[1], "vs", pair[2])
      diag_list[[paste(pair, collapse = "_")]] <- res_pair$diagnostics
    }
    log <- c(log, res_pair$log)
  }

  res <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()
  res <- .apply_p_adjust(res, p_col = "p_value", spec = spec)
  diagnostics <- if (length(diag_list) > 0) do.call(rbind, diag_list) else NULL

  out <- .result_object(
    data = res,
    diagnostics = diagnostics,
    summaries = NULL,
    meta = list(spec = spec, call = match.call()),
    log = log
  )
  if (return == "data") return(out$data)
  out
}

#' Association Scan Across Outcomes
#'
#' @description
#' Tests association between a predictor \code{x} and many outcomes \code{y_vars},
#' with optional covariates and interaction terms.
#'
#' @param df Data frame.
#' @param y_vars Character vector of outcome variables.
#' @param x Character. Predictor variable.
#' @param spec Analysis specification.
#' @param by_group Logical. If \code{TRUE}, run separately within each group.
#' @param interaction Logical. If \code{TRUE}, include \code{x * group_col}.
#' @param return Character. \code{"result"} or \code{"data"}.
#'
#' @return A \code{profile_result} object or results data frame.
#' @export
assoc_scan <- function(df, y_vars, x, spec = analysis_spec(),
                       by_group = FALSE, interaction = FALSE,
                       return = c("result", "data")) {
  return <- match.arg(return)
  .ensure_cols(df, c(y_vars, x))
  if (!is.null(spec$group_col) && (by_group || interaction)) {
    .ensure_cols(df, spec$group_col)
  }

  if (by_group && is.null(spec$group_col)) {
    stop("by_group requires spec$group_col.", call. = FALSE)
  }

  res_list <- list()
  diag_list <- list()
  log <- list()

  if (by_group) {
    for (g in unique(df[[spec$group_col]])) {
      df_g <- df[df[[spec$group_col]] == g, , drop = FALSE]
      res_g <- assoc_scan(df_g, y_vars, x, spec, by_group = FALSE, interaction = interaction, return = "result")
      res_g$data$group <- g
      res_list[[as.character(g)]] <- res_g$data
      if (!is.null(res_g$diagnostics)) {
        res_g$diagnostics$group <- g
        diag_list[[as.character(g)]] <- res_g$diagnostics
      }
      log <- c(log, res_g$log)
    }
  } else {
    for (v in y_vars) {
      cols <- unique(c(v, x, spec$covariates, spec$group_col))
      data_v <- df[, intersect(cols, names(df)), drop = FALSE]
      keep <- .complete_cases(data_v)
      n_dropped <- sum(!keep)
      if (n_dropped > 0) {
        log <- .log_add(log, "W_DROPPED_ROWS_COMPLETE_CASE",
                        "Rows dropped due to missing data.", list(var = v, n_dropped = n_dropped))
      }
      data_v <- data_v[keep, , drop = FALSE]
      if (nrow(data_v) == 0) next

      data_v <- .apply_ref_levels(data_v, spec$ref_levels)

      covars <- spec$covariates %||% character(0)
      # enforce reference level for group if present
      if (!is.null(spec$group_col) && spec$group_col %in% names(data_v)) {
        grp_levels <- levels(factor(data_v[[spec$group_col]]))
        ref_level <- .resolve_group_ref(spec, grp_levels)
        data_v[[spec$group_col]] <- .set_ref_factor(data_v[[spec$group_col]], ref_level, spec$group_col)
      }

      rhs_terms <- c(.quote_var(x), vapply(covars, .quote_var, character(1)))
      if (interaction) {
        rhs_terms <- c(rhs_terms, .quote_var(spec$group_col),
                       paste0(.quote_var(x), ":", .quote_var(spec$group_col)))
      }
      if (!is.null(spec$random_effects)) rhs_terms <- c(rhs_terms, spec$random_effects)
      f <- stats::as.formula(paste(.quote_var(v), "~", paste(rhs_terms, collapse = " + ")))

      fit <- .fit_model(f, data_v, spec)
      if (length(fit$warnings) > 0) {
        for (w in fit$warnings) log <- .log_add(log, .classify_warning(w), w, list(var = v))
      }

      td <- .standardize_tidy(.tidy_model(fit$model, fit$engine, spec$conf_int, spec$conf_level))
      term_pattern <- if (interaction) {
        paste0("^", x, "$|", x, ":")
      } else {
        paste0("^", x, "$")
      }
      td <- .term_filter(td, pattern = term_pattern)
      if (nrow(td) == 0) next

      td$y_var <- v
      td$x_var <- x
      td$model_engine <- fit$engine
      td$family <- if (is.null(fit$family)) NA_character_ else fit$family$family
      td$formula_text <- paste(deparse(f), collapse = "")
      td$n <- nrow(data_v)

      if (spec$effect_size == "standardized_beta" && spec$effect_size_from == "model") {
        es <- .standardized_beta(fit$model, data_v, x)
        td$es <- es$es
        td$es_type <- es$es_type
        td$es_notes <- NULL
      } else {
        td$es <- NA_real_
        td$es_type <- spec$effect_size
        td$es_notes <- NA_character_
      }

      if (spec$diagnostics != "none") {
        diag_fun <- switch(
          fit$engine,
          lm = .diagnostics_lm,
          glm = .diagnostics_glm,
          lmer = .diagnostics_lmer,
          glmer = .diagnostics_glmer
        )
        diag_list[[v]] <- diag_fun(fit$model, data_v, term = x, level = spec$diagnostics)
        diag_list[[v]]$y_var <- v
      }

      res_list[[v]] <- td
    }
  }

  res <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()
  res <- .apply_p_adjust(res, p_col = "p_value", spec = spec)
  diagnostics <- if (length(diag_list) > 0) do.call(rbind, diag_list) else NULL

  out <- .result_object(
    data = res,
    diagnostics = diagnostics,
    summaries = NULL,
    meta = list(spec = spec, call = match.call()),
    log = log
  )
  if (return == "data") return(out$data)
  out
}

#' Correlation Scan
#'
#' @description
#' Computes correlations between multiple variables and a target variable.
#' Supports Pearson/Spearman and partial correlation via residualization.
#'
#' @param df Data frame.
#' @param vars Character vector of variables to correlate with \code{x}.
#' @param x Character. Target variable.
#' @param method Character. \code{"pearson"} or \code{"spearman"}.
#' @param covariates Optional covariates for partial correlation.
#' @param return Character. \code{"result"} or \code{"data"}.
#'
#' @return A \code{profile_result} object or results data frame.
#' @export
corr_scan <- function(df, vars, x, method = c("pearson", "spearman"),
                      covariates = NULL, return = c("result", "data")) {
  return <- match.arg(return)
  method <- match.arg(method)
  .ensure_cols(df, c(vars, x))

  res_list <- list()
  log <- list()

  for (v in vars) {
    cols <- unique(c(v, x, covariates))
    data_v <- df[, intersect(cols, names(df)), drop = FALSE]
    keep <- .complete_cases(data_v)
    if (sum(!keep) > 0) {
      log <- .log_add(log, "W_DROPPED_ROWS_COMPLETE_CASE",
                      "Rows dropped due to missing data.", list(var = v))
    }
    data_v <- data_v[keep, , drop = FALSE]
    if (nrow(data_v) == 0) next

    data_v <- .apply_ref_levels(data_v, spec$ref_levels)

    xv <- data_v[[x]]
    yv <- data_v[[v]]

    if (!is.null(covariates) && length(covariates) > 0) {
      f_x <- .build_formula(x, NULL, covariates, NULL)
      f_y <- .build_formula(v, NULL, covariates, NULL)
      xv <- stats::resid(stats::lm(f_x, data = data_v))
      yv <- stats::resid(stats::lm(f_y, data = data_v))
    }

    ct <- stats::cor.test(xv, yv, method = method)
    res_list[[v]] <- data.frame(
      var = v,
      x_var = x,
      estimate = unname(ct$estimate),
      statistic = unname(ct$statistic),
      p_value = ct$p.value,
      n = length(xv),
      method = method,
      stringsAsFactors = FALSE
    )
  }

  res <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()
  out <- .result_object(data = res, meta = list(call = match.call()), log = log)
  if (return == "data") return(out$data)
  out
}
