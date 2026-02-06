# Orchestrator: config-first pipeline runner

.is_abs_path <- function(p) {
  if (!is.character(p) || length(p) != 1) return(FALSE)
  grepl("^([A-Za-z]:[\\\\/]|/|\\\\\\\\)", p)
}

.ensure_dirs <- function(out_root) {
  dirs <- list(
    out_root = out_root,
    data = file.path(out_root, "data"),
    tables = file.path(out_root, "tables"),
    figures = file.path(out_root, "figures"),
    logs = file.path(out_root, "logs")
  )
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  dirs
}

.roi_to_ggseg_label <- function(roi, drop_suffix = "_thickavg") {
  lab <- roi
  lab <- sub("^L_", "lh_", lab)
  lab <- sub("^R_", "rh_", lab)
  if (!is.null(drop_suffix) && nzchar(drop_suffix)) {
    lab <- sub(paste0(drop_suffix, "$"), "", lab)
  }
  lab
}

.normalize_group_column <- function(x) {
  # deterministic helper for common case-control encodings
  if (all(unique(x) %in% c(0, 1, "0", "1"))) {
    return(factor(as.character(x), levels = c("0", "1"), labels = c("Control", "Case")))
  }
  factor(as.character(x))
}

.step_time <- function(t0, t1) {
  as.numeric(difftime(t1, t0, units = "secs"))
}

.ctx_log_add <- function(ctx, code, message, context = list()) {
  ctx$log <- .log_add(ctx$log %||% list(), code, message, context)
  ctx
}

.build_formula_from_covars <- function(covars) {
  if (is.null(covars) || length(covars) == 0) return(NULL)
  stats::as.formula(paste("~", paste(vapply(covars, .quote_var, character(1)), collapse = " + ")))
}

.export_plots <- function(plots, out_dir, formats = c("png"), overwrite = FALSE) {
  if (is.null(plots) || length(plots) == 0) return(invisible(list()))
  .check_pkg("ggplot2", "plot saving")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  written <- list()

  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (!inherits(p, "ggplot")) next
    # lightweight defaults
    width <- if (grepl("brain", nm, ignore.case = TRUE)) 8 else 6
    height <- if (grepl("brain", nm, ignore.case = TRUE)) 4 else 4

    for (fmt in formats) {
      fn <- file.path(out_dir, paste0(nm, ".", fmt))
      if (!overwrite && file.exists(fn)) {
        stop(sprintf("E_EXISTS: Plot exists (set overwrite=TRUE): %s", fn), call. = FALSE)
      }
      ggplot2::ggsave(fn, plot = p, width = width, height = height, bg = "white")
      written[[nm]] <- c(written[[nm]] %||% character(0), fn)
    }
  }
  invisible(written)
}

.write_json <- function(x, path, overwrite = FALSE) {
  if (!overwrite && file.exists(path)) {
    stop(sprintf("E_EXISTS: File exists (set overwrite=TRUE): %s", path), call. = FALSE)
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  .check_pkg("jsonlite", "writing JSON")
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

.pipeline_prep <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data
  group_col <- cfg$resolved$group_col

  # basic column validation early
  required <- unique(c(cfg$columns$site, group_col, cfg$columns$age, cfg$columns$sex, cfg$columns$icv))
  .ensure_cols(df, required, context = "prep")

  # normalize common categorical columns before prep() (so ref logic is stable)
  df[[cfg$columns$site]] <- factor(as.character(df[[cfg$columns$site]]))
  df[[cfg$columns$sex]] <- factor(as.character(df[[cfg$columns$sex]]))
  df[[group_col]] <- .normalize_group_column(df[[group_col]])

  # determine group ref for consistent downstream contrast
  group_levels <- levels(df[[group_col]])
  group_ref <- params$group_ref %||% cfg$analysis$group_ref %||% .infer_ref_level(group_levels)

  factor_levels <- params$factor_levels %||% list()
  if (length(factor_levels) == 0) {
    factor_levels <- list()
    # keep Sex stable if present
    if (cfg$columns$sex %in% names(df)) {
      sex_levels <- levels(df[[cfg$columns$sex]])
      sex_ref <- params$sex_ref %||% cfg$analysis$sex_ref %||% .infer_ref_level(sex_levels, preferred = c("F", "0"))
      if (!is.null(sex_ref)) factor_levels[[cfg$columns$sex]] <- sex_ref
    }
    if (!is.null(group_ref)) factor_levels[[group_col]] <- group_ref
  }
  # Ensure explicit ref if provided, even when user passed factor_levels.
  if (!is.null(group_ref) && is.null(factor_levels[[group_col]])) {
    factor_levels[[group_col]] <- group_ref
  }

  prep_spec_obj <- prep_spec(
    roi_regex = cfg$roi$regex,
    site_col = cfg$columns$site,
    factor_levels = factor_levels,
    na_action = params$na_action %||% "drop_required",
    outlier_action = params$outlier_action %||% "none"
  )

  res <- prep(df, spec = prep_spec_obj, return = "result", verbose = isTRUE(params$verbose %||% FALSE))

  ctx$data <- res$data
  ctx$specs$prep <- res$spec
  ctx$results$prep <- res
  ctx$log <- c(ctx$log, res$log)
  ctx$meta$n_rows_after_prep <- nrow(res$data)
  ctx
}

.pipeline_combat <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data

  group_col <- cfg$resolved$group_col
  roi_cols <- cfg$resolved$roi_cols
  .ensure_cols(df, unique(c(cfg$columns$site, group_col, cfg$columns$age, cfg$columns$sex, cfg$columns$icv, roi_cols)),
               context = "combat")

  covars <- params$covariates %||% c(cfg$columns$age, cfg$columns$sex, cfg$columns$icv)
  covars <- covars[!is.na(covars) & nzchar(covars)]
  mod_formula <- .build_formula_from_covars(covars)

  combat_spec_obj <- combat_spec(
    roi_cols = roi_cols,
    roi_regex = cfg$roi$regex,
    batch_col = cfg$columns$site,
    mod_formula = mod_formula,
    missing_fit = params$missing_fit %||% "impute",
    single_batch = params$single_batch %||% "error",
    keep_inputs = TRUE
  )

  fit <- combat_fit_df(df, spec = combat_spec_obj, verbose = isTRUE(params$verbose %||% TRUE))
  res <- combat_apply_df(fit, df, verbose = isTRUE(params$verbose %||% TRUE))

  ctx$data <- res$data
  ctx$specs$combat <- fit$spec
  ctx$results$combat_fit <- fit
  ctx$results$combat_apply <- res
  ctx$log <- c(ctx$log, fit$log, res$log)

  ctx$meta$batch_levels <- fit$model$tmp$levels_batch %||% levels(factor(df[[cfg$columns$site]]))
  ctx
}

.pipeline_compare_groups <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data

  group_col <- cfg$resolved$group_col
  roi_cols <- cfg$resolved$roi_cols
  covars <- params$covariates %||% c(cfg$columns$age, cfg$columns$sex, cfg$columns$icv)
  covars <- covars[!is.na(covars) & nzchar(covars)]

  df[[group_col]] <- factor(as.character(df[[group_col]]))
  levels_g <- levels(df[[group_col]])
  ref <- params$ref %||% .infer_ref_level(levels_g)
  if (is.null(ref)) .cfg_stop("E_REF_NOT_SET", sprintf("No reference level for %s.", group_col))
  target <- params$target %||% setdiff(levels_g, ref)[1]
  if (is.null(target) || !target %in% levels_g) {
    .cfg_stop("E_TARGET_NOT_FOUND", sprintf("Target level '%s' not found in %s.", target, group_col))
  }

  spec <- analysis_spec(
    group_col = group_col,
    covariates = covars,
    model_engine = params$model_engine %||% "lm",
    diagnostics = params$diagnostics %||% "light",
    effect_size = params$effect_size %||% "hedges_g_pooled",
    effect_size_ci = params$effect_size_ci %||% "wald",
    effect_size_scope = "two_group_only",
    effect_size_from = params$effect_size_from %||% "raw_data",
    p_adjust = list(method = params$p_adjust_method %||% "fdr",
                    scope = "within_call",
                    family_id = params$family_id %||% "ct_group_diff",
                    family_desc = params$family_desc %||% "CT case-control"),
    ref_levels = setNames(list(ref), group_col),
    contrast = list(var = group_col, ref = ref, target = target, direction = "target_minus_ref")
  )

  res <- compare_groups(df, vars = roi_cols, spec = spec, return = "result")

  ctx$specs$analysis <- spec
  ctx$results$compare_groups <- res
  ctx$log <- c(ctx$log, res$log)

  # derived group counts for manifest
  ctx$meta$group_counts <- as.list(table(df[[group_col]]))
  ctx
}

.pipeline_pairwise_groups <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data

  group_col <- cfg$resolved$group_col
  roi_cols <- cfg$resolved$roi_cols
  covars <- params$covariates %||% c(cfg$columns$age, cfg$columns$sex, cfg$columns$icv)
  covars <- covars[!is.na(covars) & nzchar(covars)]

  df[[group_col]] <- factor(as.character(df[[group_col]]))

  spec <- analysis_spec(
    group_col = group_col,
    covariates = covars,
    model_engine = params$model_engine %||% "lm",
    diagnostics = params$diagnostics %||% "none",
    effect_size = params$effect_size %||% "hedges_g_pooled",
    effect_size_ci = params$effect_size_ci %||% "wald",
    effect_size_scope = "pairwise",
    effect_size_from = "raw_data",
    p_adjust = list(method = params$p_adjust_method %||% "fdr",
                    scope = "within_call",
                    family_id = params$family_id %||% "pairwise_groups",
                    family_desc = params$family_desc %||% "Pairwise group comparisons"),
    ref_levels = params$ref_levels %||% list(),
    contrast = list(var = group_col, ref = NULL, target = NULL, direction = "target_minus_ref")
  )

  res <- pairwise_groups(df, vars = roi_cols, spec = spec, return = "result")

  ctx$specs$pairwise <- spec
  ctx$results$pairwise_groups <- res
  ctx$log <- c(ctx$log, res$log)
  ctx
}

.pipeline_assoc_scan <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data

  roi_cols <- cfg$resolved$roi_cols
  x <- params$x %||% NULL
  if (is.null(x) || !is.character(x) || length(x) != 1) {
    .cfg_stop("E_CFG_MISSING", "assoc_scan requires params.x (string).")
  }
  covars <- params$covariates %||% c(cfg$columns$age, cfg$columns$sex, cfg$columns$icv)
  covars <- covars[!is.na(covars) & nzchar(covars)]

  spec <- analysis_spec(
    group_col = cfg$resolved$group_col,
    covariates = covars,
    model_engine = params$model_engine %||% "lm",
    diagnostics = params$diagnostics %||% "none",
    effect_size = params$effect_size %||% "none",
    p_adjust = list(method = params$p_adjust_method %||% "fdr",
                    scope = "within_call",
                    family_id = params$family_id %||% "assoc_scan",
                    family_desc = params$family_desc %||% "Association scan"),
    ref_levels = params$ref_levels %||% list()
  )

  res <- assoc_scan(df, y_vars = roi_cols, x = x, spec = spec,
                    by_group = isTRUE(params$by_group %||% FALSE),
                    interaction = isTRUE(params$interaction %||% FALSE),
                    return = "result")

  ctx$specs$assoc <- spec
  ctx$results$assoc_scan <- res
  ctx$log <- c(ctx$log, res$log)
  ctx
}

.pipeline_corr_scan <- function(ctx, params) {
  cfg <- ctx$cfg
  df <- ctx$data

  roi_cols <- cfg$resolved$roi_cols
  x <- params$x %||% NULL
  if (is.null(x) || !is.character(x) || length(x) != 1) {
    .cfg_stop("E_CFG_MISSING", "corr_scan requires params.x (string).")
  }

  method <- params$method %||% "pearson"
  covars <- params$covariates %||% character(0)

  res <- corr_scan(df, vars = roi_cols, x = x, method = method, covariates = covars, return = "result")

  ctx$results$corr_scan <- res
  ctx$log <- c(ctx$log, res$log)
  ctx
}

.pipeline_plots <- function(ctx, params) {
  cfg <- ctx$cfg
  enabled <- params$enabled %||% c("brain_map", "pvalue_hist")
  if (is.character(enabled) && length(enabled) == 1) enabled <- c(enabled)

  allowed <- c("brain_map", "pvalue_hist", "top_rois", "raincloud_roi", "scatter_assoc")
  unknown <- setdiff(enabled, allowed)
  if (length(unknown) > 0) {
    .cfg_stop("E_CFG_UNKNOWN", sprintf(
      "Unknown plot(s) in plots.enabled: %s. Allowed: %s.",
      paste(unknown, collapse = ", "),
      paste(allowed, collapse = ", ")
    ))
  }

  res_group <- ctx$results$compare_groups$data %||% NULL
  res_pair <- ctx$results$pairwise_groups$data %||% NULL
  res_assoc <- ctx$results$assoc_scan$data %||% NULL

  # Default table for generic plots (prefer compare_groups, then assoc_scan, then pairwise)
  res_default <- res_group
  if (!is.data.frame(res_default) || nrow(res_default) == 0) res_default <- res_assoc
  if ((!is.data.frame(res_default) || nrow(res_default) == 0) &&
      is.data.frame(res_pair) && nrow(res_pair) > 0) {
    res_default <- res_pair
  }

  plots <- list()

  if ("brain_map" %in% enabled) {
    if (is.data.frame(res_group) && nrow(res_group) > 0) {
      stat_col <- params$stat_col %||% "es_value"
      p_col <- params$p_type %||% "p_adj"
      p_max <- params$p_max %||% NULL

      legend_title <- params$legend_title %||% if (stat_col %in% c("es", "es_value") && "es_method" %in% names(res_group)) {
        unique(res_group$es_method)[1] %||% "Effect size"
      } else {
        stat_col
      }

      brain_spec <- params$brain_plot_spec %||% list()
      brain_spec$layout <- params$brain_layout %||% brain_spec$layout %||% "dispersed"
      brain_spec$view <- params$brain_view %||% brain_spec$view %||% NULL
      brain_spec$hemisphere <- params$brain_hemisphere %||% brain_spec$hemisphere %||% NULL
      brain_spec$limit <- params$brain_limit %||% brain_spec$limit %||% NULL
      brain_spec$midpoint <- brain_spec$midpoint %||% 0
      brain_spec$p_col <- brain_spec$p_col %||% p_col
      brain_spec$p_max <- brain_spec$p_max %||% p_max
      brain_spec$legend_title <- brain_spec$legend_title %||% legend_title
      brain_spec$title <- brain_spec$title %||% (params$title %||% "Brain Map")
      brain_spec$na_color <- brain_spec$na_color %||% (params$na_color %||% "grey90")
      # Brain maps use ggseg-native themes; do not apply a generic ggplot theme here.
      if (is.null(brain_spec$ggseg_theme) && !is.null(params$ggseg_theme)) brain_spec$ggseg_theme <- params$ggseg_theme

      plots$brain_map <- plot_brain_map_results(
        res_group,
        roi_col = "var",
        value_col = stat_col,
        atlas = params$atlas %||% "dk",
        roi_drop_suffix = cfg$roi$drop_suffix,
        plot_spec = brain_spec
      )
    }
  }

  if ("pvalue_hist" %in% enabled) {
    if (is.data.frame(res_default) && nrow(res_default) > 0) {
      p_col <- params$p_col %||% "p_adj"
      spec <- params$pvalue_plot_spec %||% list()
      if (is.null(spec$theme) && !is.null(params$theme)) spec$theme <- params$theme
      plots$pvalue_hist <- plot_qc_pvalue_hist(
        res_default,
        p_col = p_col,
        title = params$pvalue_title %||% "QC: p-value distribution",
        plot_spec = spec
      )
    }
  }

  if ("top_rois" %in% enabled) {
    if (is.data.frame(res_default) && nrow(res_default) > 0) {
      spec <- params$top_plot_spec %||% list()
      if (is.null(spec$theme) && !is.null(params$theme)) spec$theme <- params$theme

      metric <- params$top_metric %||% {
        if ("es_value" %in% names(res_default)) "es_value" else if ("estimate" %in% names(res_default)) "estimate" else NULL
      }
      if (is.null(metric)) .cfg_stop("E_CFG_MISSING", "top_rois requires top_metric (or a table with es_value/estimate).")
      label_col <- params$top_label %||% {
        if ("var" %in% names(res_default)) "var" else if ("y_var" %in% names(res_default)) "y_var" else "var"
      }
      plots$top_rois <- plot_top_rois(
        res_default,
        metric = metric,
        n = params$top_n %||% 20,
        label_col = label_col,
        title = params$top_title %||% "Top ROIs",
        plot_spec = spec
      )
    }
  }

  if ("raincloud_roi" %in% enabled) {
    rc <- params$raincloud_roi %||% list()
    roi <- rc$roi %||% NULL
    if (is.null(roi) || !is.character(roi) || length(roi) != 1) {
      .cfg_stop("E_CFG_MISSING", "plots.raincloud_roi requires params.raincloud_roi.roi (string).")
    }
    group_col <- rc$group_col %||% cfg$resolved$group_col
    covars <- rc$covariates %||% c(cfg$columns$age, cfg$columns$sex, cfg$columns$icv)
    covars <- covars[!is.na(covars) & nzchar(covars)]
    value <- rc$value %||% "raw"
    orientation <- rc$orientation %||% "vertical"
    style <- rc$style %||% "half_violin"
    pairwise <- rc$pairwise %||% NULL
    spec <- rc$plot_spec %||% list()
    if (is.null(spec$theme) && !is.null(params$theme)) spec$theme <- params$theme

    plots$raincloud_roi <- plot_raincloud_roi(
      ctx$data,
      roi = roi,
      group_col = group_col,
      covariates = covars,
      value = value,
      orientation = orientation,
      style = style,
      pairwise = pairwise,
      plot_spec = spec
    )
  }

  if ("scatter_assoc" %in% enabled) {
    sc <- params$scatter_assoc %||% list()
    roi <- sc$roi %||% NULL
    x <- sc$x %||% NULL
    if (is.null(roi) || !is.character(roi) || length(roi) != 1) {
      .cfg_stop("E_CFG_MISSING", "plots.scatter_assoc requires params.scatter_assoc.roi (string).")
    }
    if (is.null(x) || !is.character(x) || length(x) != 1) {
      .cfg_stop("E_CFG_MISSING", "plots.scatter_assoc requires params.scatter_assoc.x (string).")
    }
    covars <- sc$covariates %||% character(0)
    covars <- covars[!is.na(covars) & nzchar(covars)]
    value <- sc$value %||% "raw"
    color <- sc$color %||% NULL
    spec <- sc$plot_spec %||% list()
    if (is.null(spec$theme) && !is.null(params$theme)) spec$theme <- params$theme

    plots$scatter_assoc <- plot_scatter_assoc(
      ctx$data,
      roi = roi,
      x = x,
      covariates = covars,
      color = color,
      value = value,
      plot_spec = spec
    )
  }

  ctx$plots <- c(ctx$plots %||% list(), plots)
  ctx
}

.pipeline_export <- function(ctx, params) {
  cfg <- ctx$cfg
  dirs <- ctx$paths %||% .ensure_dirs(cfg$io$out_root)
  ctx$paths <- dirs

  overwrite <- isTRUE(cfg$io$overwrite)

  # Data export
  if (isTRUE(cfg$io$save_data)) {
    base <- tools::file_path_sans_ext(basename(cfg$io$input))
    suffix <- if (!is.null(ctx$results$combat_apply)) "_combat" else "_final"
    out_csv <- file.path(dirs$data, paste0(base, suffix, ".csv"))
    if (!overwrite && file.exists(out_csv)) {
      stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s", out_csv), call. = FALSE)
    }
    utils::write.csv(ctx$data, out_csv, row.names = FALSE)
    key <- if (!is.null(ctx$results$combat_apply)) "combat" else "final"
    ctx$exports$data <- ctx$exports$data %||% list()
    ctx$exports$data[[key]] <- out_csv
  }

  # Tables export (known analysis outputs)
  if (isTRUE(cfg$io$save_tables)) {
    if (!is.null(ctx$results$compare_groups$data)) {
      export_results_tables(
        ctx$results$compare_groups$data,
        out_dir = dirs$tables,
        base_name = params$base_name %||% "ct_group_diff",
        formats = params$formats %||% c("csv"),
        style = params$style %||% "publication",
        overwrite = overwrite,
        roi_drop_suffix = cfg$roi$drop_suffix
      )
      if (!is.null(ctx$results$compare_groups$summaries)) {
        sum_path <- file.path(dirs$tables, "ct_group_summary.csv")
        if (!overwrite && file.exists(sum_path)) {
          stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s", sum_path), call. = FALSE)
        }
        .write_table(ctx$results$compare_groups$summaries, sum_path, format = "csv")
      }
      if (!is.null(ctx$results$compare_groups$diagnostics)) {
        diag_path <- file.path(dirs$tables, "ct_diagnostics.csv")
        if (!overwrite && file.exists(diag_path)) {
          stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s", diag_path), call. = FALSE)
        }
        .write_table(ctx$results$compare_groups$diagnostics, diag_path, format = "csv")
      }
    }

    if (!is.null(ctx$results$pairwise_groups$data)) {
      export_results_tables(
        ctx$results$pairwise_groups$data,
        out_dir = dirs$tables,
        base_name = "pairwise_groups",
        formats = params$formats %||% c("csv"),
        style = params$style %||% "publication",
        overwrite = overwrite,
        roi_drop_suffix = cfg$roi$drop_suffix
      )
    }

    if (!is.null(ctx$results$assoc_scan$data)) {
      assoc_path <- file.path(dirs$tables, "assoc_scan.csv")
      if (!overwrite && file.exists(assoc_path)) {
        stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s", assoc_path), call. = FALSE)
      }
      .write_table(ctx$results$assoc_scan$data, assoc_path, format = "csv")
    }

    if (!is.null(ctx$results$corr_scan$data)) {
      corr_path <- file.path(dirs$tables, "corr_scan.csv")
      if (!overwrite && file.exists(corr_path)) {
        stop(sprintf("E_EXISTS: Output exists (set overwrite=TRUE): %s", corr_path), call. = FALSE)
      }
      .write_table(ctx$results$corr_scan$data, corr_path, format = "csv")
    }
  }

  # Plots export
  if (isTRUE(cfg$io$save_plots)) {
    .export_plots(ctx$plots %||% list(), dirs$figures,
                  formats = cfg$io$plot_formats %||% c("png"),
                  overwrite = overwrite)
  }

  # Logs + manifest
  if (isTRUE(cfg$manifest$enabled)) {
    manifest <- build_manifest(ctx)
    ctx$manifest <- manifest
    write_manifest(manifest,
                   file.path(cfg$io$out_root, cfg$manifest$path),
                   overwrite = overwrite)
  }
  if (!is.null(ctx$log)) {
    .write_json(ctx$log, file.path(dirs$logs, "run_log.json"), overwrite = overwrite)
  }

  ctx
}

#' Run a Config-First Analysis Pipeline
#'
#' @description
#' Executes an ordered pipeline described by a JSON config (or an equivalent R
#' list). Steps operate on a shared context (data + metadata + results).
#'
#' @param cfg_or_path Either a config list or a path to a JSON config file.
#'
#' @return A pipeline context list with fields: \code{data}, \code{results},
#'   \code{plots}, \code{specs}, \code{log}, \code{timings}, and \code{manifest}
#'   (if enabled).
#' @export
run_pipeline <- function(cfg_or_path) {
  cfg <- cfg_or_path
  if (is.character(cfg_or_path)) cfg <- read_cfg(cfg_or_path)
  cfg <- validate_cfg(cfg)

  if (!file.exists(cfg$io$input)) .cfg_stop("E_INPUT_NOT_FOUND", sprintf("Input file not found: %s", cfg$io$input))

  # load data early so resolve_cfg can auto-detect columns
  df0 <- utils::read.csv(cfg$io$input, stringsAsFactors = FALSE)
  cfg <- resolve_cfg(cfg, df0)

  ctx <- list(
    cfg = cfg,
    input = list(path = cfg$io$input, data = df0),
    data = df0,
    specs = list(),
    results = list(),
    plots = list(),
    exports = list(),
    manifest = NULL,
    log = list(),
    timings = list(),
    paths = .ensure_dirs(cfg$io$out_root),
    meta = list(n_rows_input = nrow(df0))
  )

  for (st in cfg$pipeline) {
    name <- st$step
    params <- st$params %||% list()
    t0 <- Sys.time()
    ctx <- switch(
      name,
      prep = .pipeline_prep(ctx, params),
      combat = .pipeline_combat(ctx, params),
      compare_groups = .pipeline_compare_groups(ctx, params),
      pairwise_groups = .pipeline_pairwise_groups(ctx, params),
      assoc_scan = .pipeline_assoc_scan(ctx, params),
      corr_scan = .pipeline_corr_scan(ctx, params),
      plots = .pipeline_plots(ctx, params),
      export = .pipeline_export(ctx, params),
      .cfg_stop("E_CFG_STEP", sprintf("Unsupported step: %s", name))
    )
    t1 <- Sys.time()
    ctx$timings[[name]] <- list(start = as.character(t0), end = as.character(t1), seconds = .step_time(t0, t1))
  }

  ctx
}

#' Package-Oriented Alias for \code{run_pipeline}
#'
#' @description
#' Convenience alias: \code{roiflow(cfg)} is identical to
#' \code{\link{run_pipeline}(cfg)}.
#'
#' @param cfg_or_path Config list or JSON path.
#'
#' @return Same as \code{\link{run_pipeline}}.
#' @export
roiflow <- function(cfg_or_path) {
  run_pipeline(cfg_or_path)
}
