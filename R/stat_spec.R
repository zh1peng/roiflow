# Specification helpers for profiling analyses

#' Create Analysis Specification
#'
#' @description
#' Builds a specification object controlling model engine, covariates, effect sizes,
#' diagnostics, and multiple-testing behavior for profiling analyses.
#'
#' @param group_col Character. Grouping column for group comparisons.
#' @param covariates Character vector. Covariate column names.
#' @param model_engine Character. One of \code{"auto","lm","glm","lmer","glmer"}.
#' @param family Optional family for GLM/GLMER. Either a family object or a string
#'   such as \code{"gaussian"} or \code{"binomial"}.
#' @param random_effects Character. Random effects term, e.g. \code{"(1|Site)"}.
#' @param site_col Optional site column name.
#' @param site_effect Character. How to use \code{site_col}: \code{"none"},
#'   \code{"fixed"}, or \code{"random"}.
#' @param diagnostics Character. One of \code{"none","light","full"}.
#' @param effect_size Character. Effect size type. Default is Hedges' g for two groups.
#' @param effect_size_ci Character. Effect size CI method: \code{"wald"} (default),
#'   \code{"boot"}, or \code{"none"}.
#' @param effect_size_from Character. \code{"raw_data"} or \code{"model"}.
#' @param effect_size_scope Character. \code{"two_group_only"}, \code{"pairwise"},
#'   or \code{"none"}.
#' @param p_adjust List. Multiple-testing specification with fields:
#'   \code{method}, \code{scope}, \code{family_id}, \code{family_desc}.
#' @param conf_int Logical. Whether to compute confidence intervals for model terms.
#' @param conf_level Numeric. Confidence level.
#' @param ref_levels Named list. Reference levels for categorical variables
#'   (e.g., \code{list(Group = "Control")}).
#' @param contrast List. Contrast specification with fields \code{var}, \code{ref},
#'   \code{target}, \code{direction}. Direction must be \code{"target_minus_ref"}.
#'
#' @return A list used by high-level analysis functions.
#' @export
analysis_spec <- function(
  group_col = "groups",
  covariates = c("Age", "Sex", "ICV"),
  model_engine = c("auto", "lm", "glm", "lmer", "glmer"),
  family = NULL,
  random_effects = NULL,
  site_col = NULL,
  site_effect = c("none", "fixed", "random"),
  diagnostics = c("none", "light", "full"),
  effect_size = c("hedges_g_pooled", "cohen_d_pooled", "glass_delta",
                  "standardized_beta", "none"),
  effect_size_ci = c("wald", "boot", "none"),
  effect_size_from = c("raw_data", "model"),
  effect_size_scope = c("two_group_only", "pairwise", "none"),
  p_adjust = list(method = "fdr", scope = "within_call",
                  family_id = NA_character_, family_desc = NA_character_),
  conf_int = FALSE,
  conf_level = 0.95,
  ref_levels = list(),
  contrast = list(var = NULL, ref = NULL, target = NULL, direction = "target_minus_ref")
) {
  if (is.null(contrast$var)) contrast$var <- group_col
  list(
    group_col = group_col,
    covariates = covariates,
    model_engine = match.arg(model_engine),
    family = family,
    random_effects = random_effects,
    site_col = site_col,
    site_effect = match.arg(site_effect),
    diagnostics = match.arg(diagnostics),
    effect_size = match.arg(effect_size),
    effect_size_ci = match.arg(effect_size_ci),
    effect_size_from = match.arg(effect_size_from),
    effect_size_scope = match.arg(effect_size_scope),
    p_adjust = p_adjust,
    conf_int = conf_int,
    conf_level = conf_level,
    ref_levels = ref_levels,
    contrast = contrast
  )
}
